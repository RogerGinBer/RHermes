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


#' @importFrom dplyr rename distinct_at
preprocessing <- function(struct){
    ## Fix to keep backwards compatibility with already processed files
    if ("EnviPatMass" %in% colnames(DB(struct)) &
       (!"ExactMass" %in% colnames(DB(struct)))) {
        dplyr::rename(DB(struct), "ExactMass" = "EnviPatMass")
    }
    if (!all(c("MolecularFormula", "ExactMass") %in% colnames(DB(struct)))) {
        stop("DB does not have the columns MolecularFormula and ExactMass, please rerun setDB")
    }
    F_DB <- struct@metadata@ExpParam@DB[,c("MolecularFormula", "ExactMass")]
    F_DB <- dplyr::distinct_at(F_DB, "MolecularFormula", .keep_all = T)
    colnames(F_DB) <- c("fms", "m")

    IF_DB <- IonicForm(F_DB, struct@metadata@ExpParam@adlist)

    IC <- IsoCalc(
        IF_DB[[1]], FWHM = struct@metadata@ExpParam@res,
        intTHR = 0.02, kTHR = 1, instr = struct@metadata@ExpParam@instr,
        refm = 200
    )
    return(list(IF_DB, IC))
}

#'@importFrom data.table data.table
import_and_filter <- function(lf, minpks = 20, noise = 1000) {
    ## Open the connection to a single mzML file
    fileml <- mzR::openMSfile(lf)
    plist <- mzR::peaks(fileml)
    h <- mzR::header(fileml)

    ## Filter empty headers 
    h <- h[, -which(vapply(h, function(x) all(x == 0), FUN.VALUE = logical(1)))]
    if (any(h$peaksCount < minpks)) {
        #Remove scans with very very few peaks detected
        plist <- plist[-which(h$peaksCount < minpks)]
        h <- h[-which(h$peaksCount < minpks), ]
    }
    #Remove all MS>1 scans
    if(any(h$msLevel != 1)){
        plist <- plist[-which(h$msLevel != 1)]
        h <- h[-which(h$msLevel != 1), ]
    }
    raw <- lapply(seq_along(plist), function(x) {
        #Extract raw data into a DT
        rpeaks <- plist[[x]]
        rt <- h$retentionTime[x]
        return(data.table(mz = rpeaks[, 1], rtiv = rpeaks[, 2], rt = rt))
    })
    raw <- do.call(rbind, raw)
    filtered <- raw[raw$rtiv > noise, ]
    return(list(raw, h, filtered))
}

IonicForm <- function(F_DB, Ad_DB) {
    RES <- lapply(seq_len(nrow(F_DB)), calculate_ionic_forms,
                        F_DB = F_DB,
                        Ad_DB = Ad_DB)
    db <- do.call(rbind, lapply(RES, function(x) {x[[1]]}))
    db <- db[!duplicated(db[, 1]), ]
    connections <- do.call(rbind, lapply(RES, function(x) {x[[2]]}))
    return(list(db, connections))
}

#' @importFrom MetaboCoreUtils addElements subtractElements
calculate_ionic_forms <- function(i, F_DB, Ad_DB){
    f <- as.character(F_DB$fms[i])
    j <- apply(Ad_DB, 1, function(x) {
        current_f <- f
        if (x[3] != 1) current_f <- multform(current_f, as.numeric(x[3]))
        if (x[6] != "FALSE") current_f <- addElements(current_f, x[6])
        if (x[7] != "FALSE") current_f <- subtractElements(current_f, x[7])
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
    
    testres <- lapply(testiso,
                        isocalc_parallel, kTHR, resol_factor,
                        isotopecode, isOrbi)

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
    breaks <- which(diff(x[,1]) > 5 * median(diff(x[,1]))) #Find iso clusters
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
    pmz <- as.matrix(raw[, 1])
    curiso <- IsoList[[1]][DB$f]
    PLresults <- lapply(seq_len(nrow(DB)), function(i) {
        n <- DB$f[i]
        mass <- DB$m[i]
        formula <- as.character(n)
        if (!labelled) {
            regularProc(DB, mass, formula, pmz, curiso, ppm,
                        IsoList, minhit, i, raw)
        } else {
            numC <- DB$numC[i]
            labelledProc(DB, mass, formula, pmz, curiso, ppm,
                         IsoList, minhit, i, raw, numC)
        }
    })
    PLresults <- rbindlist(PLresults)
    
    #Output coherence with PLProcesser input
    return(PLresults[, c(3, 2, 4, 5, 1)])
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
    scanid$isov <- "M0"
    isom <- mass + (isodf$deltam/abs(ch[[1]]))

    valid <- isom >= pmz[1, 1] & isom <= pmz[nrow(pmz), 1]
    if(any(valid)){
        isoss <- lapply(isom[valid], function(m) {
            #Added small multiplicative factor to ppm. We've seen that
            #isotope peaks may have a bit more error than M0
            binarySearch(pmz, m, ppm * 1.5)
        })
        isol <- vapply(isoss, function(x) length(x) != 0, FUN.VALUE = T)

        if (any(isol)) {
            isoentries <- do.call(rbind,
                                  lapply(which(isol),
                                         function(x) {
                                             isoid <- raw[isoss[[x]]]
                                             isoid$isov <- as.character(isodf$ID[valid][x])
                                             return(isoid)
                                         }
                                  )
            )
            isoentries <- isoentries[isoentries$rt %in% unique(scanid$rt), ]
            scanid <- rbind(scanid, isoentries)
        }
    }
    scanid$formv <- formula
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
    i <- floor(nrow(plist)/2)
    ub <- nrow(plist)
    lb <- 1
    cycles <- 1
    out <- FALSE
    while (cycles < 30) {
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

