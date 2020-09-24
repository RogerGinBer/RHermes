#'@import doSNOW
#'@import foreach
#'@import parallel
IsoCalc <- function(DB, FWHM, intTHR, kTHR, instr = "Orbitrap",
    refm = 200, BiocParallelParam = BiocParallel::SerialParam()) {
    data(isotopes, package = "enviPat")
    isotopecode <- data.frame(name = c("13C", "17O", "18O", "2H",
        "15N", "33S", "34S", "36S", "37Cl", "81Br", "41K", "6Li",
        "10B", "21Ne", "22Ne", "25Mg", "26Mg", "29Si", "30Si",
        "42Ca", "43Ca", "44Ca", "48Ca", "46Ti", "47Ti", "49Ti",
        "50Ti", "50Cr", "53Cr", "54Cr", "54Fe", "57Fe", "58Fe",
        "60Ni", "61Ni", "62Ni", "64Ni", "65Cu", "66Zn", "67Zn",
        "68Zn", "70Zn", "76Se", "77Se", "78Se", "82Se", "84Sr",
        "86Sr", "87Sr", "91Zr", "92Zr", "94Zr", "96Zr"), code = c("M",
        "[17O]", "[18O]", "D", "[15N]", "[33S]", "[34S]", "[36S]",
        "[37Cl]", "[81Br]", "[41K]", "[6Li]", "[10B]", "[21Ne]",
        "[22Ne]", "[25Mg]", "[26Mg]", "[29Si]", "[30Si]", "[42Ca]",
        "[43Ca]", "[44Ca]", "[48Ca]", "[46Ti]", "[47Ti]", "[49Ti]",
        "[50Ti]", "[50Cr]", "[53Cr]", "[54Cr]", "[54Fe]", "[57Fe]",
        "[58Fe]", "[60Ni]", "[61Ni]", "[62Ni]", "[64Ni]", "[65Cu]",
        "[66Zn]", "[67Zn]", "[68Zn]", "[70Zn]", "[76Se]", "[77Se]",
        "[78Se]", "[82Se]", "[84Sr]", "[86Sr]", "[87Sr]", "[91Zr]",
        "[92Zr]", "[94Zr]", "[96Zr]"), stringsAsFactors = FALSE)

    message("Starting Envipat isotope distribution calculation: ")
    DB <- as.data.frame(DB)

    isOrbi <- ifelse(instr == "Orbitrap", TRUE, FALSE)
    if (isOrbi) {
        resol_factor <- FWHM/sqrt(refm)
    } else {
        resol_factor <- FWHM
    }

    testiso <- enviPat::isopattern(isotopes = isotopes, threshold = intTHR,
        chemforms = DB[, "envi"], charge = DB[, "ch"], verbose = TRUE)
    names(testiso) <- DB[, "f"]
    rm(DB) #Free up space for Windows SOCK users
    message("")
    message("Calculating distinguishable isotopes given the instrumental resolution: ")

    nsplit <- BiocParallelParam$workers
    suppressWarnings({
        testres <- bplapply(testiso,
                            RHermes:::isocalc_parallel, kTHR, resol_factor,
                            isotopecode, isOrbi, BPPARAM = BiocParallelParam)
    }) #Suppress warnings to avoid the split() "object not multiple of ..."
    # testres <- unlist(testres, recursive = FALSE)

    b <- lapply(testres, function(x) {
        return(any(is.na(x)))
    })

    #All different isopologues detected in the ionic formula set and their
    #deltaM with respect to M0
    factordf <- do.call(rbind, testres)
    factordf <- factordf[!duplicated(factordf[, 1]), ]

    #To save space in the list (avoiding multiple separate factor instances) we
    #just save the number that is represented by the string in the factor object
    facts <- factor(unlist(factordf[, 1]), levels = unlist(factordf[, 1]),
                    ordered = FALSE)
    d4 <- lapply(testres, function(x) {
        return(as.numeric(factor(unlist(x$ID), levels = levels(facts))))
    })

    #Removing entries with NA to avoid errors downstream
    bad <- which(vapply(d4, function(x) {
        any(is.na(unlist(x)))
    }, logical(1)))
    if (length(bad) != 0) {
        d4 <- d4[-bad]
    }
    return(list(d4, factordf))
}

isocalc_parallel <- function(x, kTHR, resol_factor, isotopecode, isOrbi){
    if (isOrbi) {
        limitfactor <- 2 * kTHR * x[1, 1]^(1/2)/resol_factor
    } else {
        limitfactor <- 2 * kTHR * x[1, 1]/resol_factor
    }
    farenough <- diff(x[, 1]) > limitfactor
    good <- c()
    for (i in seq_len(nrow(x))) {
        if (i == nrow(x)) {
            good <- c(good, i)
            next
        }
        if (farenough[i]) {
            good <- c(good, i)
        } else {
            if (x[i, 2] > x[i + 1, 2]) {
                good <- c(good, i)
            }
        }
    }
    good <- unique(good)
    x <- x[good, ]

    curiso <- isotopecode[which(isotopecode[, 1] %in% colnames(x)),]
    x <- as.data.frame(x)
    x$ID <- ""
    for (i in seq_len(nrow(curiso))) {
        col <- which(colnames(x) == curiso[i, 1])[1]
        num <- x[, col]
        tomodify <- which(num != 0)
        x$ID[tomodify] <- paste0(x$ID[tomodify],
                                 paste0(curiso[i,2], num[tomodify]))
    }
    x$deltam <- x[, 1] - as.numeric(unlist(x[1, 1]))
    x <- x[-1, ]
    return(x[, c("ID", "deltam")])
}
