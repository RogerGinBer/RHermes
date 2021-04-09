#'@title processMS1
#'@description Main function to process the input MS1 mzML files. It
#' accesses the files with mzR and processes them with RHermes internal
#' functions. For each file a PL object is generated inside the
#' object_data_PL slot.
#'
#'@param struct RHermesExp S4 object to update with the processed data.
#' Important: The objects needs to have the experimental parameters already
#' set before calling this function.
#'@param files Character vector of all the paths to the files to be processed.
#'@param labelled Logical, whether to check for all 13C isotopic signals.
#' Defaults to FALSE
#'@return An RHermesExp object with the processed PLs.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))
#'          MS1file <- system.file("extdata", "MS1TestData.mzML",
#'                              package = "RHermes")
#'          }
#'processMS1(struct, MS1file, FALSE) #For regular use
#'processMS1(struct, MS1file, TRUE) #For labelled data
#'

#'@export

setGeneric("processMS1", function(struct, files, labelled = FALSE) {
    standardGeneric("processMS1")
})
#' @rdname processMS1
setMethod("processMS1", c("RHermesExp", "character", "ANY"),
function(struct, files, labelled = FALSE) {
    ppm <- struct@metadata@ExpParam@ppm
    noise <- struct@metadata@ExpParam@nthr

    message(paste0("Preprocessing steps, calculating possible ionic formulas ",
                    "and isotopic distributions"))

    prepro <- preprocessing(struct)
    IF_DB <- prepro[[1]]
    IC <- prepro[[2]]

    message("Starting file processing...")
    struct <- setTime(struct, "Started file processing into PL")
    toAdd <- lapply(seq_along(files), function(i) {
        lf = files[i]
        message(paste0("Now processing: ", lf))
        imported <- import_and_filter(lf, 20, noise)
        ss <- OptScanSearch(DB = IF_DB[[1]],
                            raw = imported[[3]],
                            ppm = ppm,
                            labelled = labelled,
                            IsoList = IC)

        #Construction of S4 Object output
        RHermesPL(peaklist = ss, header = imported[[2]], raw = imported[[1]],
                    labelled = labelled, filename = lf)
    })
    struct@data@PL <- c(struct@data@PL, toAdd)
    struct@metadata@ExpParam@ionF <- IF_DB
    struct@metadata@ExpParam@isoList <- IC
    struct@metadata@filenames <- c(struct@metadata@filenames, files)

    struct <- setTime(struct,
                    paste("Ended file processing. The following files were",
                        "successfully processed:",
                        paste(files, collapse = "  ")
                    ))
    return(struct)
})


preprocessing <- function(struct){
    F_DB <- struct@metadata@ExpParam@DB[,c("MolecularFormula", "EnviPatMass")]
    #Could break if colname isn't exactly "MolecularFormula"
    F_DB <- dplyr::distinct_at(F_DB, "MolecularFormula",
                                .keep_all = T)
    colnames(F_DB) <- c("fms", "m")

    IF_DB <- IonicForm(F_DB, struct@metadata@ExpParam@adlist)

    IC <- IsoCalc(
        IF_DB[[1]], FWHM = struct@metadata@ExpParam@res,
        intTHR = 0.02, kTHR = 1, instr = struct@metadata@ExpParam@instr,
        refm = 200
    )
    return(list(IF_DB, IC))
}

import_and_filter <- function(lf, minpks = 20, noise = 1000) {
    #Opening the connection to a single mzML file
    fileml <- mzR::openMSfile(lf)
    plist <- mzR::peaks(fileml)
    h <- mzR::header(fileml)
    #Filtering header
    h <- h[, -which(vapply(h, function(x) all(x == 0), FUN.VALUE = logical(1)))]
    if (any(h$peaksCount < minpks)) {
        #Removing scans with very very few peaks detected
        plist <- plist[-which(h$peaksCount < minpks)]
        h <- h[-which(h$peaksCount < minpks), ]
    }
    raw <- lapply(seq_along(plist), function(x) {
        #Extracting raw data into a DT
        rpeaks <- plist[[x]]
        rt <- h$retentionTime[x]
        return(data.table(mz = rpeaks[, 1], rtiv = rpeaks[, 2], rt = rt))
    })
    raw <- do.call(rbind, raw)
    filtered <- raw[raw$rtiv > noise, ]
    return(list(raw, h, filtered))
}

IonicForm <- function(F_DB, Ad_DB) {
    suppressWarnings({
        RES <- bplapply(seq_len(nrow(F_DB)), calculate_ionic_forms,
                            BPPARAM = bpparam(), F_DB = F_DB,
                            Ad_DB = Ad_DB)
    })
    db <- do.call(rbind, lapply(RES, function(x) {x[[1]]}))
    db <- db[!duplicated(db[, 1]), ]
    connections <- do.call(rbind, lapply(RES, function(x) {x[[2]]}))
    return(list(db, connections))
}

calculate_ionic_forms <- function(i, F_DB, Ad_DB){
    f <- as.character(F_DB$fms[i])
    j <- apply(Ad_DB, 1, function(x) {
        current_f <- f
        if (x[3] != 1) current_f <- multform(current_f, as.numeric(x[3]))
        if (x[6] != "FALSE") current_f <- sumform(current_f, x[6])
        if (x[7] != "FALSE") current_f <- subform(current_f, x[7])
        if (is.na(current_f)) return(NA)

        ch <- ifelse(x[5] == "positive", "+", "-")
        current_f <- paste0("[", current_f, "]",
                            ifelse(abs(as.numeric(x[2])) == 1,
                                    ch,
                                    ifelse(as.numeric(x[2]) > 0, c(x[2], ch),
                                            c(strsplit(x[2], "-")[[1]][2], ch))
                                    )
                            )
        return(current_f)
    })
    good <- which(!is.na(j))
    j <- j[!is.na(j)]
    if(length(j) == 0){return()}
    envi <- strsplit(j, split = "[", fixed = TRUE) %>%
        vapply(function(x) {x[[2]]}, FUN.VALUE = character(1)) %>%
        strsplit(j, split = "]", fixed = TRUE) %>%
        vapply(function(x) {x[[1]]}, FUN.VALUE = character(1))


    charge <- abs(as.numeric(Ad_DB[good, 2]))
    multiplicity <- as.numeric(Ad_DB[good,3])
    adduct_delta <- as.numeric(Ad_DB[good,4])

    mass <- ((F_DB$m[i] * multiplicity + adduct_delta) / charge) %>%
                round(., digits = 5)
    db <- data.table(f = j, m = mass, ch = charge, envi = envi)
    con <- data.table(f = rep(f, length(j)), ion = j, an = Ad_DB[good, 1])
    return(list(db, con))
}

sumform <- function(f1, f2) {
    f1 <- tryCatch({
        CHNOSZ::count.elements(f1)
    }, error = function(cond){return(NA)})
    if (any(is.na(f1))) return(NA)
    n1 <- names(f1)
    f2 <- CHNOSZ::count.elements(f2)
    n2 <- names(f2)
    interAB <- n1[which(n1 %in% n2)]
    if (length(interAB) != 0) {
        s1 <- vapply(interAB, function(x) {
            f1[[x]] + f2[[x]]
        }, numeric(1))
        names(s1) <- interAB
    } else {
        s1 <- integer()
    }
    uniqueA <- n1[which(!(n1 %in% n2))]
    if (length(uniqueA) != 0) {s1 <- c(s1, f1[uniqueA])}
    uniqueB <- n2[which(!(n2 %in% n1))]
    if (length(uniqueB) != 0) {s1 <- c(s1, f2[uniqueB])}
    alln <- unique(c(n1, n2))
    alln <- alln[!(alln %in% c("C", "H"))]
    ord <- c()
    if ("C" %in% unique(c(n1, n2))) {ord <- c("C")}
    if ("H" %in% unique(c(n1, n2))) {ord <- c(ord, "H")}
    ord <- c(ord, sort(alln))
    s1 <- s1[ord]

    output <- paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(s1), s1), collapse = "")

    return(output)
}

subform <- function(f1, f2) {
    f1 <- tryCatch({
        CHNOSZ::count.elements(f1)
    }, error = function(cond){return(NA)})
    if (any(is.na(f1))) {return(NA)}
    n1 <- names(f1)
    f2 <- CHNOSZ::count.elements(f2)
    n2 <- names(f2)
    interAB <- n1[which(n1 %in% n2)]
    if (length(interAB) < length(n2)) {return(NA)}
    if (length(interAB) != 0) {
        s1 <- vapply(interAB, function(x) {
            res <- ifelse(f1[[x]] - f2[[x]] > 0, f1[[x]] - f2[[x]], NA)
        }, numeric(1))
        if (any(is.na(s1))) {return(NA)}
        if (any(s1 == 0)) {
            n1 <- n1[which(n1 %in% n2)[s1 != 0]]
            n2 <- n2[which(n2 %in% n1)[s1 != 0]]
            s1 <- s1[s1 != 0]
        }
        names(s1) <- interAB
    } else {
        s1 <- integer()
    }
    uniqueA <- n1[which(!(n1 %in% n2))]
    if (length(uniqueA) != 0) {s1 <- c(s1, f1[uniqueA])}
    uniqueB <- n2[which(!(n2 %in% n1))]
    if (length(uniqueB) != 0) {s1 <- c(s1, f2[uniqueB])}
    alln <- unique(c(n1, n2))
    alln <- alln[!(alln %in% c("C", "H"))]
    ord <- c()
    if ("C" %in% unique(c(n1, n2))) {ord <- c("C")}
    if ("H" %in% unique(c(n1, n2))) {ord <- c(ord, "H")}
    ord <- c(ord, sort(alln))
    s1 <- s1[ord]

    output <- paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(s1), s1), collapse = "")

    return(output)
}

multform <- function(f, k) {
    f <- CHNOSZ::count.elements(f)
    f <- f * k
    paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(f), f), collapse = "")
}

IsoCalc <- function(DB, FWHM, intTHR, kTHR, instr = "Orbitrap", refm = 200) {
    isotopes <- NULL #To appease R CMD Check "no visible binding"
    data(isotopes, package = "enviPat", envir = environment())
    isotopecode <- data.frame(
        name = c("13C", "17O", "18O", "2H",
            "15N", "33S", "34S", "36S", "37Cl", "81Br", "41K", "6Li",
            "10B", "21Ne", "22Ne", "25Mg", "26Mg", "29Si", "30Si",
            "42Ca", "43Ca", "44Ca", "48Ca", "46Ti", "47Ti", "49Ti",
            "50Ti", "50Cr", "53Cr", "54Cr", "54Fe", "57Fe", "58Fe",
            "60Ni", "61Ni", "62Ni", "64Ni", "65Cu", "66Zn", "67Zn",
            "68Zn", "70Zn", "76Se", "77Se", "78Se", "82Se", "84Sr",
            "86Sr", "87Sr", "91Zr", "92Zr", "94Zr", "96Zr"),
        code = c("M","[17O]", "[18O]", "D", "[15N]", "[33S]", "[34S]", "[36S]",
            "[37Cl]", "[81Br]", "[41K]", "[6Li]", "[10B]", "[21Ne]",
            "[22Ne]", "[25Mg]", "[26Mg]", "[29Si]", "[30Si]", "[42Ca]",
            "[43Ca]", "[44Ca]", "[48Ca]", "[46Ti]", "[47Ti]", "[49Ti]",
            "[50Ti]", "[50Cr]", "[53Cr]", "[54Cr]", "[54Fe]", "[57Fe]",
            "[58Fe]", "[60Ni]", "[61Ni]", "[62Ni]", "[64Ni]", "[65Cu]",
            "[66Zn]", "[67Zn]", "[68Zn]", "[70Zn]", "[76Se]", "[77Se]",
            "[78Se]", "[82Se]", "[84Sr]", "[86Sr]", "[87Sr]", "[91Zr]",
            "[92Zr]", "[94Zr]", "[96Zr]"),
        stringsAsFactors = FALSE
    )

    message("Starting Envipat isotope distribution calculation: ")
    DB <- as.data.frame(DB)

    isOrbi <- ifelse(instr == "Orbitrap", TRUE, FALSE)
    if (isOrbi) {
        resol_factor <- FWHM/sqrt(refm)
    } else {
        resol_factor <- FWHM
    }

    testiso <- enviPat::isopattern(isotopes = isotopes, threshold = intTHR,
                                chemforms = DB$envi, charge = DB$ch,
                                verbose = TRUE)
    names(testiso) <- DB$f
    rm(DB) #Free up space for Windows SOCK users
    message("")
    message("Calculating isotopes given the instrumental resolution: ")

    suppressWarnings({
        testres <- bplapply(testiso,
                            isocalc_parallel, kTHR, resol_factor,
                            isotopecode, isOrbi, BPPARAM = bpparam())
    }) #Suppress warnings to avoid the split() "object not multiple of ..."

    #All different isopologues detected in the ionic formula set and their
    #deltaM with respect to M0
    factordf <- do.call(rbind, testres)
    factordf <- factordf[!duplicated(factordf[, 1]), ]

    #To save space in the list we  save the number that is
    #represented by the level in the factor object
    facts <- factor(unlist(factordf[, 1]), levels = unlist(factordf[, 1]),
                    ordered = FALSE)
    d4 <- lapply(testres, function(x) {
        return(as.numeric(factor(unlist(x$ID), levels = levels(facts))))
    })

    #Removing entries with NA to avoid errors downstream
    bad <- which(vapply(d4, function(x) {any(is.na(unlist(x)))}, logical(1)))
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
    breaks <- which(diff(x[,1]) > 5*median(diff(x[,1]))) #Find iso clusters
    breaks <- data.frame(st = c(1, breaks + 1), end = c(breaks, nrow(x)))
    clust <- lapply(seq_len(nrow(breaks)), function(i){
        x[seq(breaks$st[i], breaks$end[i]), , drop = FALSE]
    })

    #Split the cluster into subclusters and take the most intense peaks
    x <- do.call(rbind, lapply(clust, function(x){
        farenough <- diff(x[, 1]) > limitfactor
        if (!any(farenough)) {
            return(x[which.max(x[,2]), ]) #No subclusters
        }
        farenough <- data.frame(st = c(1, which(farenough) + 1),
                                end = c(which(farenough), nrow(x)))
        clust <- lapply(seq_len(nrow(farenough)), function(i){
            cur <- x[seq(farenough$st[i], farenough$end[i]), , drop = FALSE]
            return(cur[which.max(cur[,2]), ])
        })
        do.call(rbind, clust)
    }))

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

OptScanSearch <- function(DB, raw, ppm, IsoList, labelled = FALSE,
                            minhit = 1) {
    message(paste0("This process can take quite a bit of time, depending on ",
                    "the processing power and RAM your computer has"))
    setkeyv(raw, c("mz"))
    DB <- as.data.table(DB)  #Making sure its a DT
    setkeyv(DB, c("f"))
    if (labelled) {
        DB$numC <- as.numeric(lapply(DB$envi, function(x){
            elements <- CHNOSZ::count.elements(x)
            if ("C" %in% names(elements)) {return(elements[["C"]])}
            else{return(0)}
        }))
    }

    ncores <- ifelse(is.numeric(bpparam()$workers[[1]]),
                    yes = bpparam()$workers, no = 1)

    #Splitting the formulas into a list (with l = number of workers) to reduce
    #time loss associated with variable loading (in SOCK only)
    flist <- split(DB, f = seq_len(ncores))
    PLresults <- bplapply(seq_along(flist), PLparallelfun, flist, raw, IsoList,
                        labelled, ppm, minhit, BPPARAM = bpparam())
    PLresults <- do.call(rbind, PLresults)

    #Output coherence with PLProcesser input
    return(PLresults[, c(3, 2, 4, 5, 1)])
}


PLparallelfun <- function(gr, flist, raw, IsoList, labelled, ppm, minhit){
    curDB <- flist[[gr]]
    pmz <- as.matrix(raw[, 1])
    curiso <- IsoList[[1]][curDB$f]
    localRES <- lapply(seq_len(nrow(curDB)), function(i) {
        n <- curDB$f[i]
        mass <- curDB$m[i]
        formula <- as.character(n)
        if (!labelled) {
            regularProc(curDB, mass, formula, pmz, curiso, ppm,
                                    IsoList, minhit, i, raw)
        } else {
            numC <- curDB$numC[i]
            labelledProc(curDB, mass, formula, pmz, curiso, ppm,
                                    IsoList, minhit, i, raw, numC)
        }
    })
    return(do.call(rbind, localRES))
}

regularProc <- function(curDB, mass, formula, pmz, curiso, ppm, IsoList, minhit,
                        i, raw){
    if (mass < pmz[1, 1] | mass > pmz[nrow(pmz), 1]) {return()}
    ss <- binarySearch(pmz, mass, ppm)
    if (length(ss) < minhit) {return()}  #Return nothing if no M0 hit

    isofactors <- curiso[[i]]  #If hit, let's find the isotopologues
    isodf <- IsoList[[2]][isofactors, ]
    ch <- curDB[i, 3]  #Charge to normalize isotope deltam's

    scanid <- raw[ss]
    scanid$formv <- formula
    scanid$isov <- "M0"

    isom <- mass + (isodf$deltam/abs(ch[[1]]))
    isoss <- lapply(isom, function(m) {
        if (m < pmz[1, 1] | m > pmz[nrow(pmz), 1]) {return()}
        #Added small multiplicative factor to ppm. We've seen that
        #isotope peaks may have a bit more error than M0
        binarySearch(pmz, m, ppm * 1.5)
    })
    isol <- vapply(isoss, length, FUN.VALUE = numeric(1))

    isoidx <- which(isol != 0)
    if (length(isoidx) != 0) {
        isoentries <- do.call(rbind,
            lapply(isoidx,
                function(x) {
                    isoid <- raw[isoss[[x]]]
                    isoid$formv <- formula
                    isoid$isov <- as.character(isodf$ID[x])
                    return(isoid)
                }
            )
        )
        isoentries <- isoentries[isoentries$rt %in% unique(scanid$rt), ]
        scanid <- rbind(scanid, isoentries)
    }
    return(scanid)
}

labelledProc <- function(curDB, mass, formula, pmz, curiso, ppm,
                            IsoList, minhit, i, raw, numC){
    if (mass < pmz[1, 1] | mass > pmz[nrow(pmz), 1]) {
        ss <- numeric()
    } else {
        ss <- binarySearch(pmz, mass, ppm)
    }
    if (length(ss) == 0) {return()}

    ##Note that we don't need any M0 hit to find the other signals
    isofactors <- curiso[[i]]
    isodf <- IsoList[[2]][isofactors, ]
    if (numC != 0) {
        isodf <- rbind(isodf, data.frame(ID = paste0("M", seq_len(numC)),
                                        deltam = 1.003355*seq_len(numC)))
        isodf <- dplyr::distinct(isodf)
    }
    ch <- curDB[i, 3]  #Charge to normalize isotope deltam's


    scanid <- raw[ss, ]
    if (nrow(scanid) != 0) {
        scanid$formv <- formula
        scanid$isov <- "M0"
    }

    isom <- mass + (isodf$deltam/abs(ch[[1]]))
    isoss <- lapply(isom, function(m) {
        if (m < pmz[1, 1] | m > pmz[nrow(pmz), 1]) {return()}
        #Added small multiplicative factor to ppm. We've seen that
        #isotope peaks may have a bit more error than M0
        binarySearch(pmz, m, ppm * 1.5)
    })
    isol <- vapply(isoss, length, FUN.VALUE = numeric(1))
    isoidx <- which(isol != 0)
    if (length(isoidx) != 0) {
        isoentries <- do.call(rbind,
            lapply(isoidx,
                function(x) {
                    isoid <- raw[isoss[[x]]]
                    isoid$formv <- formula
                    isoid$isov <- as.character(isodf$ID[x])
                    return(isoid)
                }
            )
        )
        if (nrow(scanid) == 0) {return(isoentries)}
        return(rbind(scanid, isoentries))
    }
    if (nrow(scanid) == 0) {return()}
    return(scanid)
}

binarySearch <- function(plist, m, ppm) {
    i <- round((m - plist[1, 1])/(plist[nrow(plist), 1] - plist[1,1]) *
                nrow(plist))[[1]]  #Biased start based on target
    if (i == 0) {i <- 1}
    ub <- nrow(plist)
    lb <- 1
    cycles <- 1
    out <- FALSE
    while (cycles < 50) {
        dist <- (plist[i, 1] - m)/m * 1e+06
        if (i == ub) {
            break
        } else if (dist >= ppm) {
            ub <- i
            i <- round(lb + ((ub - lb)/2))
            cycles <- cycles + 1
        } else if (dist <= -ppm) {
            lb <- i
            i <- round(lb + ((ub - lb)/2))
            cycles <- cycles + 1
        } else {
            out <- TRUE
            top <- i
            bot <- i
            #This 'for' loops find the interval of raw entries that
            #satisfy d < ppm error.
            #The 'jump' strategy used avoids unnecessary calculations
            for (j in c(1000, 100, 10, 1)) {
                repeat {
                    if (top + j > ub) {break}
                    if ((plist[top + j, 1] - m)/m * 1e+06 <= ppm) {
                        top <- top + j
                    } else {
                        break
                    }
                }
            }
            for (j in c(1000, 100, 10, 1)) {
                repeat {
                    if (bot - j < lb) {
                        break
                    }
                    if ((plist[bot - j, 1] - m)/m * 1e+06 >= -ppm) {
                        bot <- bot - j
                    } else {
                        break
                    }
                }
            }
            break
        }
    }
    if (!out) {return()}
    return(bot:top)
}



