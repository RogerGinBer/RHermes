#' @title exportMSP
#' @family Exports
#' @param struct RHermesExp object
#' @param id Index of the MS2Exp that you want to export
#' @param file Name of the output file, without the .msp termination
#' @param whichSpec Which entries to export. Defaults to all of them.
#' @description  Exports the superspectra of a given MS2Exp ID into an
#' .msp format file
#' @return An .msp file
#' @examples
#'   if(FALSE){
#'        exportMSP(struct, 1, "output")
#'    }
#' @export
exportMSP <- function(struct, id, file, whichSpec = NA) {
    il <- struct@data@MS2Exp[[id]]@IL
    ilannot <- il@annotation
    il <- il@IL
    MS2Features <- struct@data@MS2Exp[[id]]@Ident$MS2Features
    charge <- struct@metadata@ExpParam@ion
    
    if (is.na(whichSpec)) {
        whichSpec <- seq_len(nrow(MS2Features))
    }
    
    sink(file = paste0(file, ".msp"))
    for (i in whichSpec) {
        ilentry <- MS2Features$ILentry[i]
        precs <- ilannot[[ilentry]]$metadata
        precs <- sapply(precs,function(x) x$mass)
        precs <- unique(as.numeric(unlist(precs)))
        curss <- MS2Features[i,]
        for (precMZ in precs){
            ssdata <- curss$ssdata[[1]]
            ssdata$int <- 100*ssdata$int/max(ssdata$int)
            cat("BEGIN IONS\n")
            cat("TITLE=Entry", i,"Range:", curss$start[[1]], 
                "-", curss$end[[1]], "\n", sep = "_")
            cat("RTINSECONDS=", curss$apex[[1]], "\n", 
                sep = "")
            cat("PEPMASS=", precMZ, "\n", sep = "")
            cat("CHARGE=", "1", charge, "\n", sep = "")
            for (j in seq_len(nrow(ssdata))) {
                cat(as.numeric(ssdata[j, ]), "\n")
            }
            cat("END IONS\n")
            # cat("\n")
        }
    }
    sink()
}



#' @title exportMGF
#' @family Exports
#' @inheritParams exportMSP
#' @param file Name of the output file, without the .mgf termination
#' @description  Exports the superspectra of a given MS2Exp into an
#'  .mgf format file
#' @return An .msf file
#' @examples
#'   if(FALSE){
#'        exportMGF(struct, 1, "output")
#'    }
#' @export
exportMGF <- function(struct, id, file, whichSpec = NA) {
    sstable <- struct@data@MS2Exp[[id]]@Ident[[1]]
    charge <- struct@metadata@ExpParam@ion
    if (is.na(whichSpec)) {
        whichSpec <- seq_len(nrow(sstable))
    }
    sink(file = paste0(file, ".mgf"))
    for (i in whichSpec) {
        curss <- sstable[i, ]
        precmz <- curss$precmass
        ssdata <- curss$ssdata[[1]]
        # Write into the .ms file
        cat("BEGIN IONS\n")
        cat("PEPMASS=", precmz, "\n", sep = "")
        cat("CHARGE=", "1", charge, "\n", sep = "")
        cat("RTINSECONDS=", curss$apex[[1]], "\n", sep = "")
        cat("TITLE=Entry ", i, " Apex at ", curss$apex[[1]], " seconds. ",
            "Range: ", curss$start[[1]], "-", curss$end[[1]], "\n", sep = "")
        for (j in seq_len(nrow(ssdata))) {
            cat(as.numeric(ssdata[j, ]), "\n")
        }
        cat("END IONS\n")
        cat("\n")
    }
    sink()
}

#' @title exportmzML
#' @family Exports
#' @author Jordi Capellades
#' @inheritParams exportMSP
#' @param fname Name of the output file, without the .mzML termination
#' @description  Exports the superspectra of a given MS2Exp ID into an
#' mzML format file. Warning! The retention time information is not
#' included in the mzML file.
#' @return An .mzML file
#' @examples
#'   if(FALSE){
#'        exportmzML(struct, 1, "output")
#'    }
#' @export
exportmzML <- function (struct, id, fname, 
                        whichSpec = NA, collisionEnergy=35) 
{
    Ident <- struct@data@MS2Exp[[id]]@Ident
    il <- struct@data@MS2Exp[[id]]@IL
    ilannot <- il@annotation
    il <- il@IL
    MS2Data <- struct@data@MS2Exp[[id]]@MS2Data
    MS2Features <- Ident$MS2Features
    if (is.null(MS2Features)) {
        stop("No MS2 spectra found for that MS2Exp ID index")
    }
    out_file <- paste0(fname, ".mzML")
    if (!is.na(whichSpec)) {
        Ident <- Ident[whichSpec]
        il <- il[whichSpec, ]
        MS2Data <- MS2Data[whichSpec]
        MS2Features <- MS2Features[whichSpec]
        out_file <- paste0(fname, "_", paste(whichSpec[1], whichSpec[length(whichSpec)], 
                                             sep = "_"), ".mzML")
    }
    minMZ <- struct@metadata@ExpParam@minmz
    maxMZ <- struct@metadata@ExpParam@maxmz
    if (struct@metadata@ExpParam@ion == "+") {
        polValue <- 1
    }else {
        polValue <- 0
    }
    pks <- list()
    hdr <- structure(list(seqNum = integer(0), acquisitionNum = integer(0), 
                          msLevel = integer(0), polarity = integer(0), peaksCount = integer(0), 
                          totIonCurrent = numeric(0), retentionTime = numeric(0), 
                          basePeakMZ = numeric(0), basePeakIntensity = numeric(0), 
                          collisionEnergy = numeric(0), ionisationEnergy = numeric(0), 
                          lowMZ = numeric(0), highMZ = numeric(0), precursorScanNum = integer(0), 
                          precursorMZ = numeric(0), precursorCharge = integer(0), 
                          precursorIntensity = numeric(0), mergedScan = integer(0), 
                          mergedResultScanNum = integer(0), mergedResultStartScanNum = integer(0), 
                          mergedResultEndScanNum = integer(0), injectionTime = numeric(0), 
                          filterString = character(0), spectrumId = character(0), 
                          centroided = logical(0), ionMobilityDriftTime = numeric(0), 
                          isolationWindowTargetMZ = numeric(0), isolationWindowLowerOffset = numeric(0), 
                          isolationWindowUpperOffset = numeric(0), scanWindowLowerLimit = numeric(0), 
                          scanWindowUpperLimit = numeric(0)), row.names = integer(0), 
                     class = "data.frame")
    p <- 1
    for (i in seq_len(nrow(MS2Features))) {
        ilentry <- MS2Features$ILentry[i]
        precs <- ilannot[[ilentry]]$metadata
        precs <- sapply(precs,function(x) x$mass)
        precs <- unique(unlist(precs))
        for (precMZ in precs){
            rt <- p * 0.5
            for (l in c(1, 2)) {
                if (l == 1) {
                    precInt <- 10000
                    ssdata <- matrix(c(precMZ, precInt), nrow = 1, 
                                     ncol = 2)
                    for (r in seq(5, 1)) {
                        j <- nrow(hdr) + 1
                        hdr[j, ] <- NA
                        pks <- append(pks, list(ssdata))
                        rt2 <- rt - (0.05 * r)
                        hdr$retentionTime[j] <- rt2
                        hdr$msLevel[j] <- l
                        hdr$polarity[j] <- polValue
                        hdr$peaksCount[j] <- nrow(ssdata)
                        hdr$totIonCurrent[j] <- sum(ssdata[, 2])
                        imax <- which.max(ssdata[, 2])
                        hdr$basePeakMZ[j] <- ssdata[imax, 1]
                        hdr$basePeakIntensity[j] <- ssdata[imax, 2]
                        hdr$ionisationEnergy[j] <- 0
                        hdr$lowMZ[j] <- min(ssdata[, 1])
                        hdr$highMZ[j] <- max(ssdata[, 1])
                        hdr$precursorCharge[j] <- 0
                        hdr$injectionTime[j] <- 100
                        hdr$centroided[j] <- TRUE
                        hdr$scanWindowLowerLimit[j] <- minMZ
                        hdr$scanWindowUpperLimit[j] <- maxMZ
                        if (polValue == 1) {
                            polString <- "FTMS + p ESI Full ms ["
                        }else {
                            polString <- "FTMS - p ESI Full ms ["
                        }
                        hdr$filterString[j] <- paste0(polString, hdr$scanWindowLowerLimit[j], 
                                                      ".0000-", hdr$scanWindowUpperLimit[j], ".0000]")
                    }
                }else {
                    j <- nrow(hdr) + 1
                    hdr[j, ] <- NA
                    ssdata <- MS2Features$ssdata[[i]]
                    ssdata <- ssdata[order(ssdata$mz), ]
                    ssdata <- as.matrix(ssdata)
                    colnames(ssdata) <- NULL
                    pks <- append(pks, list(ssdata))
                    rt2 <- rt + 0.025
                    hdr$retentionTime[j] <- rt2
                    hdr$precursorMZ[j] <- precMZ
                    hdr$isolationWindowTargetMZ[j] <- precMZ
                    hdr$precursorIntensity[j] <- precInt
                    hdr$msLevel[j] <- l
                    hdr$polarity[j] <- polValue
                    hdr$peaksCount[j] <- nrow(ssdata)
                    hdr$totIonCurrent[j] <- sum(ssdata[, 2])
                    imax <- which.max(ssdata[, 2])
                    hdr$basePeakMZ[j] <- ssdata[imax, 1]
                    hdr$basePeakIntensity[j] <- ssdata[imax, 2]
                    hdr$collisionEnergy[j] <- collisionEnergy
                    hdr$ionisationEnergy[j] <- 0
                    hdr$lowMZ[j] <- min(ssdata[, 1])
                    hdr$highMZ[j] <- max(ssdata[, 1])
                    hdr$precursorCharge[j] <- 0
                    hdr$injectionTime[j] <- 50
                    hdr$centroided[j] <- TRUE
                    hdr$isolationWindowLowerOffset <- 0.8
                    hdr$isolationWindowUpperOffset <- 0.8
                    minMZfilter <- round(floor(min(hdr$lowMZ[j])),0) - 10
                    maxMZfilter <- max(c(hdr$highMZ[j], hdr$precursorMZ[j]))
                    maxMZfilter <- round(maxMZfilter, 0) + 20
                    if (polValue == 1) {
                        polString2 <- "FTMS + p ESI d Full ms2 "
                    }else {
                        polString2 <- "FTMS - p ESI d Full ms2 "
                    }
                    hdr$filterString[j] <- paste0(polString2, hdr$precursorMZ[j], 
                                                  "@hcd", hdr$collisionEnergy[j], ".00 [", minMZfilter, 
                                                  ".0000-", maxMZfilter, ".0000]")
                    hdr$scanWindowLowerLimit[j] <- minMZfilter
                    hdr$scanWindowUpperLimit[j] <- maxMZfilter
                    p <- p+1
                }
            }
            j <- nrow(hdr) + 1
            hdr[j, ] <- NA
            ssdata <- matrix(c(precMZ, precInt * 0.5), nrow = 1, 
                             ncol = 2)
            pks <- append(pks, list(ssdata))
            rt2 <- rt + 0.1
            hdr$retentionTime[j] <- rt2
            hdr$msLevel[j] <- 1
            hdr$polarity[j] <- polValue
            hdr$peaksCount[j] <- nrow(ssdata)
            hdr$totIonCurrent[j] <- sum(ssdata[, 2])
            imax <- which.max(ssdata[, 2])
            hdr$basePeakMZ[j] <- ssdata[imax, 1]
            hdr$basePeakIntensity[j] <- ssdata[imax, 2]
            hdr$ionisationEnergy[j] <- 0
            hdr$lowMZ[j] <- min(ssdata[, 1])
            hdr$highMZ[j] <- max(ssdata[, 1])
            hdr$precursorCharge[j] <- 0
            hdr$injectionTime[j] <- 100
            hdr$centroided[j] <- TRUE
            hdr$scanWindowLowerLimit[j] <- minMZ
            hdr$scanWindowUpperLimit[j] <- maxMZ
            hdr$filterString[j] <- paste0(polString, hdr$scanWindowLowerLimit[j], 
                                          ".0000-", hdr$scanWindowUpperLimit[j], ".0000]")
        }
    }
    for (i in seq_len(nrow(hdr))) {
        if (hdr$msLevel[i] == 1) {
            hdr$seqNum[i] <- i
            hdr$acquisitionNum[i] <- i
            j <- i
        }else {
            hdr$seqNum[i] <- i
            hdr$acquisitionNum[i] <- i
            hdr$precursorScanNum[i] <- j
        }
    }
    hdr$spectrumId <- paste0("controllerType=0 controllerNumber=1 scan=", 
                             hdr$seqNum)
    rownames(hdr) <- NULL
    suppressWarnings(file.remove(out_file))
    mzR::writeMSData(object = pks, file = out_file, header = hdr)
    of <- mzR::openMSfile(out_file)
    mzR::close(of)
}

#' @title exportSIRIUS
#' @family Exports
#' @description Exports all processed Hermes spectra into a .ms format suitable
#' for SIRIUS analysis.
#' @inheritParams exportMSP
#' @param file Name of the output file, without the .ms termination
#' @return Returns nothing, but generates an .ms file where specified.
#' @examples
#'   if(FALSE){
#'        exportmzML(struct, 1, "output")
#'    }
#' @export
exportSIRIUS <- function(struct, id, file, whichSpec = NA){
    Ident <- struct@data@MS2Exp[[id]]@Ident
    IL <- struct@data@MS2Exp[[id]]@IL@IL
    MS2Features <- Ident$MS2Features
    if(is.null(MS2Features)){
        stop("No MS2 spectra found for that MS2Exp ID index")
    }
    if(!is.na(whichSpec)){
        Ident <- Ident[whichSpec]
        IL <- IL[whichSpec, ]
        MS2Features <- MS2Features[whichSpec, ]
    }

    sink(file = paste0(file, ".ms"))
    for (i in seq_len(nrow(MS2Features))) {
        curentry <- MS2Features[i, ]
        fs <- curentry$anot %>% unlist()
        msms <- curentry$ssdata[[1]]
        #Write into the .ms file
        for (j in fs) {
            comp <- paste("Spec", i, j, sep = "_")
            cat(">compound", comp, "\n")
            cat(">formula", j, "\n")
            cat(">ms2", "\n")
            for (j in seq_len(nrow(msms))) {
                cat(as.numeric(msms[j, ]), "\n")
            }
            cat("\n")
        }
    }
    sink()
}


#'@title exportIdent
#'@description Export all MS2 identifications (hits) into a csv file.
#'@author Roger Gine
#'@param struct An RHermesExp object
#'@param id The id of the RHermesMS2Exp object.
#'@param file Name assigned to the csv file.
#'@return Returns nothing, but generates a csv file with the exported
#'  identifications
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#' if(FALSE){exportIdent(struct, 1, "identifications.csv")}
#'@export
setGeneric("exportIdent", function(struct, id, file){
    standardGeneric("exportIdent")
})
#' @rdname exportIdent
setMethod("exportIdent", c("RHermesExp", "numeric", "character"),
function(struct, id, file) {
    ms2 <- struct@data@MS2Exp[[id]]@Ident[[1]]
    ms2 <- ms2[vapply(ms2$results, is.data.frame, logical(1)), ]

    #Extracting hits from inner data.frame
    ms2$hits <- lapply(ms2$results, function(hits) {
        if (!is.data.frame(hits)) {return(hits)}
        hits$formula
    })
    ms2$bestscore <- lapply(ms2$results, function(hits) {
        if (!is.data.frame(hits)) {return("N/A")}
        as.character(round(max(hits$cos), digits = 3))
    })
    ms2 <- dplyr::select(ms2, -c("ssdata", "results"))

    #Adjusting format
    ms2$start <- as.numeric(ms2$start)
    ms2$end <- as.numeric(ms2$end)
    ms2$apex <- as.numeric(ms2$apex)
    ms2$precmass <- as.numeric(ms2$precmass)
    ms2$bestscore <- as.numeric(ms2$bestscore)
    ms2$start <- format(round(ms2$start, 2), nsmall = 2)
    ms2$end <- format(round(ms2$end, 2), nsmall = 2)
    ms2$apex <- format(round(ms2$apex, 2), nsmall = 2)
    ms2$bestscore <- format(round(ms2$bestscore, 4), nsmall = 4)
    ms2$precmass <- format(round(ms2$precmass, 4), nsmall = 4)

    #Collapsing multiple annotations into a single string
    ms2$anot <- lapply(ms2$anot, function(x) {
        paste(x, collapse = " ")
    })
    ms2$hits <- lapply(ms2$hits, function(x) {
        paste(x, collapse = " ")
    })
    ms2$hits <- lapply(ms2$hits, function(x) {
        gsub(pattern = "\n", replacement = "", x = x, )
    })
    ms2$hits <- lapply(ms2$hits, function(x) {
        gsub(pattern = "\t", replacement = "", x = x, )
    })
    write.csv(ms2, file = file)
})


#'@export
#' @rdname RHermesMS2Exp-class
#' @param object An RHermesMS2Exp object
setMethod("show", "RHermesMS2Exp", function(object){
    message("Info about this MS2Exp object:")
    show(object@IL)
    hasMS2 <- length(object@Ident) != 0
    if(hasMS2){
        num_entries <- length(unique(object@Ident[["MS2Features"]]$ILentry))
        num_sup <- nrow(object@Ident[["MS2Features"]])
        num_ident <- length(which(vapply(object@Ident[["MS2Features"]]$results,
                                            is.data.frame, logical(1))))
        message(paste0("Number of covered entries: ", num_entries ))
        message(paste0("Number of superspectra: ", num_sup ))
        message(paste0("Identified superspectra: ", num_ident ))
    } else {
        message("No identification was performed")
    }
})

#'@title ssNetwork
#'@description Generates a similarity matrix between different MS2 spectra.
#'@author Roger Gine
#'@inheritParams exportIdent
#'@param ss The MS2 spectra from Ident() to be compared
#'@return A square symmetrical similarity matrix
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'ssNetwork(struct, 1, 1:7)
#'@export
ssNetwork <- function(struct, id, ss) {
    validObject(struct)
    MS2Exp <- struct@data@MS2Exp[[id]]
    allspec <- MS2Exp@Ident$MS2Features$ssdata[ss]
    cos <- lapply(allspec, function(pattern) {
        lapply(allspec, function(query) {
            MSMScosineSim(pattern, query, minhits = 1)
        })
    }) %>% unlist() %>% matrix(nrow = length(allspec), byrow = TRUE)
    cos
}
