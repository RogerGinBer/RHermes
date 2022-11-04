#### findSOI-related ####

#'@title findSOI
#'@description Generates a list of SOI (Scans of Interest) given a peaklist and
#'  appends it to the RHermesExp object.
#'@param struct The RHermesExp object to which the SOI list will be saved.
#'@param params A SOIParam object or a list of SOIParam, contains all the
#'  parameters necessary to generate the SOI. See \link[RHermes]{getSOIpar} for
#'  more info on how to set it. If params is a list, findSOI will use each
#'  object for each SOI list it generates until it has used them all. When that
#'  happens, the function will keep using the same SOIparam object for the rest
#'  of SOIs.
#'@param fileID Numeric vector of the indexes of the peaklists you want to
#'  process into a SOI.
#'@param blankID Optional. Numeric vector of the indexes of peaklists to use as
#'  a 'blank' against the correspoding fileID PL. If you enter 0 (or don't input
#'  anything at all), no blank substraction will be performed. In case you want
#'  to process multiple SOIs and only some of them have blank substraction, you
#'  must space the 0s so that the indexes match. (See examples)
#'
#'@examples
#'if(FALSE){
#'p <- getSOIpar('double') #Set SOI parameters
#'myHermes <- findSOI(myHermes, p, 1) #Basic usage, no blank substraction
#'
#'#Blank substraction usage, 1 without blank and 2 using 1 as blank
#'myHermes <- findSOI(myHermes, p, c(1, 2), c(0, 1))
#'
#'#Even more advanced usage, with multiple parameters and blank substraction
#'#Here the first SOI is calculated using p and without blank sub.
#'#The second uses p2 and 1 as a blank
#'#And the same goes for the third. p2 would be reused and 1 as blank
#'p2 <- getSOIpar('triple')
#'myHermes <- findSOI(myHermes, c(p,p2), c(1, 2, 3), c(0, 1, 1))
#'}
#'
#'@return An updated RHermesExp object with the resulting SOI lists appended.
#'
#'@export
#'
#'@importFrom dplyr distinct filter
#'@import magrittr
#'@import reticulate
setGeneric("findSOI", function(struct, params, fileID, blankID = numeric(0)) {
    standardGeneric("findSOI")
})

#' @rdname findSOI
setMethod("findSOI", c("RHermesExp", "ANY", "ANY", "ANY"),
function(struct, params, fileID, blankID = numeric(0)) {
    if (length(blankID) == 0) {
        blankID <- rep(0, length(fileID))
    } else if (length(blankID) < length(fileID)){
        message("Less blanks than files, reusing the last ID provided.")
        blankID <- c(blankID,
                     rep(blankID[length(blankID)],
                         length(fileID) - length(blankID)))
    }
    maxn <- length(struct@data@PL)
    if (any(c(fileID, blankID) > maxn)) {
        stop("Some indexes are above the number of current PL")
    }
    if (is.list(params)) {
        multParams <- TRUE
        specID <- c(seq_along(params),
                    rep(length(params),
                    times = length(fileID) - length(params)))
    } else {
        multParams <- FALSE
    }
    for (i in seq_along(fileID)) {
        idx <- fileID[i]
        if (idx == blankID[i]) {
            stop(paste("You've tried to substract the blank from",
                        "the same file. This is not allowed."))
        }
        if (multParams) {
            cur <- params[[specID[i]]]
            ifelse(blankID[i] != 0, cur@blanksub <- TRUE,
                cur@blanksub <- FALSE)
        } else {
            cur <- params
            ifelse(blankID[i] != 0, params@blanksub <- TRUE,
                params@blanksub <- FALSE)
        }
        if (blankID[i] != 0) {
            cur@blanksub <- TRUE
            cur@blankname <- struct@data@PL[[blankID[i]]]@filename
            blankPL <- struct@data@PL[[blankID[i]]]@peaklist
        } else {
            cur@blanksub <- FALSE
            cur@blankname <- "None"
            blankPL <- NA
        }
        message(paste("--------", "Processing SOI",
                        i, "out of", length(fileID), "--------\n"))
        struct@data@SOI <- c(struct@data@SOI,
                            PLprocesser(struct@data@PL[[idx]],
                            struct@metadata@ExpParam, cur, blankPL,
                            struct@metadata@filenames[idx]
                            ))
        if (blankID[i] == 0) {
            struct <- setTime(struct,
                            paste("Generated SOI list from",
                            struct@metadata@filenames[idx]))
        } else {
            struct <- setTime(struct,
                            paste("Generated SOI list from",
                            struct@metadata@filenames[idx], "using",
                            struct@metadata@filenames[blankID[i]],
                            "as blank"))
        }
    }
    return(struct)
})

#' @importFrom stats  approx  median  sd
PLprocesser <- function(PL, ExpParam, SOIParam, blankPL = NA, filename) {
    ## Extracting info from S4 objects into local variables
    DataPL <- PL@peaklist
    h <- PL@header
    formulaDB <- ExpParam@ionF[[1]]
    FA_to_ion <- ExpParam@ionF[[2]]
    ppm <- ExpParam@ppm

    params <- SOIParam@specs
    maxlen <- SOIParam@maxlen
    noise <- SOIParam@minint
    useblank <- SOIParam@blanksub
    mode <- tryCatch(SOIParam@mode, error = function(cond){"regular"})

    ## Setting up PeakList
    DataPL <- as.data.table(DataPL)
    #Filtering M0 beforehand to improve indexing performance later on
    DataPL <- DataPL[DataPL$isov == "M0", ]
    DataPL <- DataPL[DataPL$rtiv > noise, ]
    setkeyv(DataPL, c("formv", "rt"))
    setkeyv(formulaDB, c("f"))
    setkeyv(FA_to_ion, c("ion"))

    if(mode == "regular"){
        ## Density filtering
        message("Starting density filtering: ")
        GR <- apply(params, 1, function(x) {
            densityProc(x, DataPL, h)
        })
        ## Grouping different filtered results
        message("Now Grouping:")
        Groups <- GR[[1]]
        setkeyv(Groups, c("formula", "start", "end"))
        if (length(GR) > 1) {
            for (i in 2:length(GR)) {
                Groups <- groupGrouper(GR, i, Groups)
            }
        }
        Groups <- distinct(Groups)
    } else {
        message("Starting Centwave-based SOI detection")
        cwp <- SOIParam@cwp
        Groups <- calculateSOICentwave(filename, formulaDB, cwp, ppm = ppm)
    }

    ## Initial Peak Retrieval
    message("Retrieving datapoints from SOIs:")
    Groups$peaks <- lapply(seq_len(nrow(Groups)), retrievePeaks,
                        Groups, DataPL)


    ## SOI Characterization
    message("Starting SOI characterization:")

    message("Mass calculation:")
    setkeyv(Groups, "formula")
    Groups$mass <- formulaDB[.(Groups$formula),2]

    message("Number of scans:")
    Groups$nscans <- apply(Groups, 1, function(x) {
        d <- x["peaks"][[1]]
        return(dim(d)[1])
    })
    Groups <- Groups[Groups$nscans > 5, ]

    if (any(Groups$length > maxlen) & mode == "regular") {
        Groups <- groupShort(Groups[,1:6], maxlen)
        Groups$nscans <- apply(Groups, 1, function(x) {
            d <- x["peaks"][[1]]
            return(dim(d)[1])
        })
        Groups <- Groups[Groups$nscans > 5, ]
    }

    ## Blank substraction
    if (useblank) {
        Groups <- blankSubstraction(Groups, blankPL)
        if (nrow(Groups) == 0) {return(RHermesSOI())}
    }

    ## Rest of characterization
    message("Width calculation:")
    width <- apply(Groups, 1, function(x) {
        d <- x["peaks"][[1]]
        sidx <- dim(d)[1]
        w <- 0
        if (sidx > 2) {
            w <- log10(max(d[, 2])/min(d[, 2]))
        }
        return(w)
    })
    Groups$width <- width

    message("Maximum intensity calculation:")
    suppressWarnings({
        MaxI <- apply(Groups, 1, function(x) {
        d <- x["peaks"][[1]]
        return(max(d$rtiv))
        })
    })
    Groups$MaxInt <- MaxI
    Groups <- Groups[Groups$MaxInt > 0, ]

    message("Converting from ionic formula to F/A combinations:")
    Groups$anot <- lapply(Groups$formula, function(x) {
        apply(FA_to_ion[.(x), c(1, 3)], 1, function(y) {
            paste(y, collapse = "_")
        })
    })

    message("Calculating chaos:")
    Groups$chaos <- lapply(Groups$peaks, function(x) {
        tryCatch({rho_chaos(x, nlevels = 10, fillGaps = TRUE)},
                    error = function(cond){return(0)})
    }) %>% as.numeric()

    setkeyv(Groups, c("formula"))
    message("Generating peaklist for plotting:")
    plist <- lapply(unique(Groups$formula), preparePlottingDF,
                        Groups)
    plist <- do.call(rbind, plist)
    plist$isov <- rep("M0", nrow(plist))

    # Constructing S4 output object
    output <- RHermesSOI(SOIList = Groups, PlotDF = as.data.table(plist),
                        SOIParam = SOIParam, filename = filename)

    return(output)
}

densityProc <- function(x, DataPL, h){
    rtbin <- x[1]
    scanspercent <- x[2]
    shift <- x[3]
    message("Running Density Filter")
    BinRes <- densityFilter(DataPL, h, rtbin, "M0", shift)
    cutoff <- BinRes[[1]] * scanspercent

    #Correcting CUT < 1 cases to avoid errors (eg. if it was 0, any time region
    #would be included). Added a 1 for robustness.
    cutoff[cutoff < 1] <- median(c(cutoff[cutoff > 1], 1))
    message("Running Density Interpreter")
    uf <- unique(DataPL$formv)
    RES <- lapply(BinRes[-c(1,2)], parallelInterpreter, cutoff)
    RES <- RES[vapply(RES, is.data.frame, logical(1))]
    if(length(RES) == 0){return()}
    RES <- lapply(names(RES), function(anot){
        RES[[anot]]$formula <- anot
        RES[[anot]]
    })
    RES <- do.call(rbind, RES)

    #Changing from index to corresponding RT
    RES[, 1] <- (RES[, 1] * rtbin) - rtbin + floor(min(h$retentionTime)) + shift
    RES[, 2] <- (RES[, 2] * rtbin) - rtbin + floor(min(h$retentionTime)) + shift
    RES %<>% dplyr::mutate(., length = .data$end - .data$start)
    RES <- RES[, c(1, 2, 4, 3)]
    return(as.data.table(RES))
    }

resolveGroup <- function(x, Groups, matched, first_df){
    corresp <- rbind(Groups[unlist(matched[.(x), 1]), ], first_df[x, ])
    if (nrow(corresp) == 1) {return(corresp)}
        corresp[1, 1] <- min(corresp[, 1])  #Min start
        corresp[1, 2] <- max(corresp[, 2])  #Max end
        return(corresp[1, ])
    }

groupGrouper <- function(GR, i, Groups){
    first_df <- GR[[i]]
    setkeyv(first_df, c("formula", "start", "end"))
    message(paste("Merging filter", i, "out of", length(GR),
                "with previous SOI list"))
    ## Overlap between current group element and the next
    # Ordered by formula, start and end time
    links <- foverlaps(Groups, first_df, by.x = c(4, 1, 2),
                    by.y = c(4, 1, 2), type = "any", which = TRUE)

    nomatched <- links[is.na(links$yid), ]$xid
    #Which are found to be overlapping in both dfs
    matched <- links[!is.na(links$yid), ]
    setkeyv(matched, "yid")

    ## Grouping entries and choosing lowest start and highest end time
    res <- lapply(unique(matched$yid), resolveGroup, Groups,
                matched, first_df)
    res <- do.call(rbind, res)

    end <- NULL; start <- NULL #To appease R CMD Check "no visible binding"
    res[, `:=`(length, end - start)]
    res <- distinct(res)

    message("Adding entries unique from the first DF")
    res <- rbind(res, Groups[nomatched, ])

    message("Adding entries unique from the second DF")
    links <- foverlaps(first_df, Groups, by.x = c(4, 1, 2),
                    by.y = c(4, 1, 2), type = "any", which = TRUE)
    other_unmatched <- links[is.na(links$yid), ]$xid
    res <- rbind(res, first_df[other_unmatched, ])

    iter <- 1
    repeat {
        setkeyv(res, c("formula", "start", "end"))
        message(paste("Solving inner conflicts, Round:", iter))
        # Inner overlaps
        links2 <- foverlaps(res, res, by.x = c(4, 1,2), by.y = c(4, 1, 2),
                            type = "any", which = TRUE)
        match2 <- links2[links2$yid %in% which(table(links2$yid) > 1), ]
        setkeyv(match2, "yid")
        if (nrow(match2) == 0) {
            break
        }
        res2 <- lapply(unique(match2$yid), function(x) {
        corresp <- res[unlist(match2[.(x), 1]), ]
        if (nrow(corresp) == 1) {
            return(corresp)
        }
        corresp[1, 1] <- min(corresp[, 1])
        corresp[1, 2] <- max(corresp[, 2])
        return(corresp[1, ])
        })
        res2 <- do.call(rbind, res2)
        res2[, `:=`(length, end - start)]
        res2 <- distinct(res2)
        res <- rbind(res[links2[links2$yid %in%
                                which(table(links2$yid) == 1)]$yid,], res2)
        iter <- iter + 1
        if (iter > 10) {break} #Break for safety
    }
    return(res)
}

#' @importFrom MSnbase readMSData
calculateSOICentwave <- function(filename, DB, CentWaveParam, ppm){
    msdata <- readMSData(filename, mode = "onDisk")
    pks <- findChromPeaks(msdata, CentWaveParam)
    pks <- as.data.frame(chromPeaks(pks))
    pks$anot <- lapply(pks$mz, function(mz){
        range <- mz * c(1 - ppm * 1e-6,
                        1 + ppm * 1e-6)
        DB$f[between(DB$m, range[1], range[2])]
    })
    pks <- pks[sapply(pks$anot, length) > 0, ]
    soi <- lapply(seq_len(nrow(pks)), function(i){
        peak <- pks[i, ]
        data.frame(start = peak$rtmin,
                   end = peak$rtmax,
                   length = peak$rtmax - peak$rtmin,
                   formula = peak$anot[[1]])
    })
    soi <- do.call(rbind, soi)
    soi <- as.data.table(soi)
    return(soi)
}

retrievePeaks <- function(i, Groups, PL, delta = 0){
    x <- Groups[i, ]
    rtmin <- x[1, 1]
    rtmax <- x[1, 2]
    f <- x[1, 4]
    pks <- PL[.(f)] 
    pks <- filter(pks, between(pks$rt, rtmin[[1]] - delta, rtmax[[1]] + delta))
    return(pks[, c(1, 2, 5)])
}

groupShort <- function(Groups, maxlen, BPPARAM = bpparam()){
    message("Shortening and selecting long groups:")
    SG <- filter(Groups, length <= maxlen)
    LG <- filter(Groups, length > maxlen)
    LG <- lapply(seq_len(nrow(LG)), parallelGroupShort, LG, maxlen)
    LG <- do.call(rbind, LG)
    return(rbind(SG, LG))
}


#Experimental Centwave approach to SOI partitioning

parallelGroupShort <- function(i, LG, maxlen){
    ms1data <- LG$peaks[[i]]
    curGR <- LG[i,]
    tryCatch({
        #Suppress the "No peaks found" warning
        suppressWarnings(
            pks <- xcms::peaksWithCentWave(int = ms1data$rtiv, rt = ms1data$rt,
                            peakwidth = c(8,60), prefilter = c(0,100),
                            snthresh = 0, noise = 0, fitgauss = FALSE,
                            firstBaselineCheck = FALSE)
        )
    }, error = function(cond){
        pks <- data.frame() #To avoid problems, act as if no peaks found
    })

    #If XCMS detects any peak
    if (nrow(pks) != 0) {
        pks <- as.data.frame(pks)
        pks <- pks[order(pks[,2]), ]
        starts <- pks[,2]
        ends <- pks[,3]
        known_peak <- rep(TRUE, nrow(pks))
        if (nrow(pks) > 1) {
            for (i in seq_len(nrow(pks) - 1)) {
                if (ends[i] < starts[i + 1]) {
                starts <- c(starts, ends[i])
                ends <- c(ends, starts[i + 1])
                known_peak <- c(known_peak, FALSE)
                }
            }
        }
        if (min(ms1data[,1]) < min(starts)) {
            ends <- c(ends, min(starts))
            starts <- c(starts, min(ms1data[, 1]))
            known_peak <- c(known_peak, FALSE)
            }
        if (max(ms1data[,1]) > max(ends)) {
            starts <- c(starts, max(ends))
            ends <- c(ends, max(ms1data[, 1]))
            known_peak <- c(known_peak, FALSE)
        }
        #Split remaining long traces into smaller ones
        if (any(!known_peak)) {
            too_long <- ends - starts > maxlen
            if (any(too_long & !known_peak)) {
                for (i in which(too_long & !known_peak)) {
                curstart <- starts[i]
                curend <- ends[i]
                times <- seq(from = curstart, to = curend,
                        length.out = ceiling((curend - curstart) / maxlen) + 1)
                starts <- c(starts, times[-length(times)])
                ends <- c(ends, times[-1])
                }
            starts <- starts[-which(too_long & !known_peak)]
            ends <- ends[-which(too_long & !known_peak)]
            }
        }
        NewGR <- data.table(start = starts, end = ends,
                            length = ends - starts, formula = curGR$formula,
                            peaks = curGR$peaks, mass = curGR$mass)
    } else {
        #Divide long traces into equal-sized smaller traces
        times <- seq(from = curGR[1, 1][[1]], to = curGR[1, 2][[1]],
                    length.out = ceiling(curGR[1, 3][[1]] / maxlen) + 1)
        deltat <- times[2] - times[1]
        NewGR <- data.table(start = times[-length(times)], end = times[-1],
                            length = deltat, formula = curGR$formula,
                            peaks = curGR$peaks, mass = curGR$mass)
    }
    #Data point redistribution within each new SOI
    NewGR[, "peaks"] <- apply(NewGR, 1, function(x) {
        pks <- x[5][[1]]
        return(pks[between(pks$rt, x[1][[1]], x[2][[1]]), ])
    })
    return(NewGR)
}

blankSubstraction <- function(Groups, blankPL){
    message("Blank substraction:")
    blankPL <- blankPL[blankPL$isov == "M0", ]
    setkeyv(blankPL, c("formv", "rt"))
    message("First cleaning")
    toKeep <- lapply(seq_len(nrow(Groups)), firstCleaning, Groups, blankPL) %>%
        unlist()
    sure <- Groups[which(toKeep), ]
    if (any(toKeep)) {Groups <- Groups[-which(toKeep),]}
    setkeyv(blankPL, "formv")
    message("Preparing input for cosine similarity")
    RES <- lapply(seq_len(nrow(Groups)), prepareNetInput, Groups, blankPL)
    cosines <- sapply(RES, function(x){
        if(is.na(x)[1]){return(Inf)}
        philentropy::cosine_dist(x[1,] + 1e-12, x[2,] + 1e-12,
                                 testNA = FALSE)
    })
    Groups <- Groups[cosines < 0.8, ]
    Groups <- rbind(sure, Groups)
    return(Groups)
}

#'@importFrom stats IQR quantile
firstCleaning <- function(i, Groups, blankPL){
    peaks <- Groups$peaks[[i]]
    deltat <- 10
    blankpks <- retrievePeaks(i, Groups, blankPL, deltat)
    blankpks <- distinct(blankpks[, c(1, 2)])
    if (nrow(blankpks) < 5) {return(TRUE)} #No blank signals

    sampleCV <- IQR(peaks$rtiv) /
        (quantile(peaks$rtiv, 0.25) + quantile(peaks$rtiv, 0.75))

    blankCV <-  IQR(blankpks$rtiv) /
        (quantile(blankpks$rtiv, 0.25) + quantile(blankpks$rtiv, 0.75))

    if(is.na(sampleCV) | is.na(blankCV)){return(FALSE)}

    sampleMax <- max(peaks$rtiv)
    # blankMax <- max(blankpks$rtiv)
    q90_ratio <- quantile(peaks$rtiv,0.9) / quantile(blankpks$rtiv,0.9)

    #We have to be restrictive with the conditions, otherwise we collect junk
    if (sampleCV/blankCV > 5) {return(TRUE)}
    if (q90_ratio > 3 & sampleMax > 15000) {return(TRUE)}

    return(FALSE) #Can't decide if the SOI is good enough, let the ANN decide
}

prepareNetInput <- function(i, Groups, blankPL){
    cur <- Groups[i, ]
    st <- cur$start
    end <- cur$end
    f <- cur$formula
    deltat <- 10
    peaks <- cur$peaks[[1]]
    Npoints <- 200
    tryCatch({
        if (nrow(peaks) <= 2 | length(unique(peaks$rt)) <= 2) {
            return(NA) #No signal
        }
        #Interpolate to N points
        smooth_pks <- data.frame(approx(x = peaks,
                                        xout = seq(from = st, to = end,
                                                    length.out = Npoints),
                                        rule = 1, ties = min))
        smooth_pks[is.na(smooth_pks[, "y"]), "y"] <- 0
        blankpks <- retrievePeaks(i, Groups, blankPL, deltat)
        blankpks <- distinct(blankpks[, c(1, 2)])
        if (length(unique(blankpks$rt)) < 2) {
            blankpks <- data.frame(rt = c(st, end), rtiv = c(0, 0))
            smooth_blankpks <- data.frame(
                approx(x = blankpks,
                    xout = seq(from = st, to = end, length.out = Npoints),
                    rule = 1, ties = min))
        } else {
            smooth_blankpks <- data.frame(
                approx(x = blankpks,
                    xout = seq(from = st, to = end, length.out = Npoints),
                    rule = 1, ties = min))
            smooth_blankpks[is.na(smooth_blankpks[, "y"]), "y"] <- 0
        }
        return(as.matrix(rbind(smooth_pks[, 2], smooth_blankpks[,2])) /
                max(smooth_pks[, 2]))
        }, error = function(cond) {
            return(NA)
    })
}

preparePlottingDF <- function(i, Groups){
    data <- Groups[.(i), ]
    res <- data.frame(rt = 0.1, rtiv = 0.1, form = "", stringsAsFactors = FALSE)
    for (j in seq_len(nrow(data))) {
        peaks <- data[j, ]$peaks[[1]]
        if (is.null(dim(peaks)[1])) {
            next
        }
        peaks <- peaks[,1:2]
        form <- rep(unlist(data[j, ]$formula), nrow(peaks))
        peaks <- cbind(peaks, form)
        res <- rbind(res, peaks)
    }
    res <- res[-1, ]
    return(res)
}

parallelFilter <- function(anot, ScanResults, bins, timebin){
    data <- ScanResults[.(anot), "rt"][[1]]
    mint <- min(data)
    maxt <- max(data)
    res <- rep(0, length(bins))

    #Shortcut to get in which bin each point is
    l <- table(ceiling((data - bins[1]) / timebin))
    if (length(l) == 0) return(res)
    res[as.numeric(names(l))] <- l
    return(res)
}

densityFilter <- function(ScanResults, h, timebin, iso = "M0", tshift = 0) {
    #Setting time bins
    rtmin <- floor(min(h$retentionTime))
    rtmax <- ceiling(max(h$retentionTime))
    bins <- seq(from = rtmin - tshift, to = rtmax + tshift, by = timebin)

    #Getting how many scans were taken on each bin
    scans <- lapply(bins, function(curbin) {
        return(h %>% filter(., .data$retentionTime > curbin &
                                .data$retentionTime <= curbin + timebin) %>%
                   dim(.) %>% .[1])
    })

    #Counting how many scan entries are on each bin
    setkeyv(ScanResults, c("formv", "rt"))
    RES <- lapply(unique(ScanResults$formv), parallelFilter, ScanResults,
                  bins, timebin)
    names(RES) <- unique(ScanResults$formv)
    RES <- c(list(unlist(scans)), list(bins), RES)
    names(RES)[c(1,2)] <- c("Bins","Scans")
    return(RES)
}

parallelInterpreter <- function(x, cutoff) {
    interlist <- densityInterpreter(x, cutoff)
    if (length(interlist[[1]]) == 0) {return()}
    res <- data.frame(start = interlist[[1]], end = interlist[[2]])
    return(res)
}

densityInterpreter <- function(list, cutoff) {
    ls <- rbind(list, cutoff)
    #Calculate bins have >thr scans (good)
    bool <- apply(ls, 2, function(x) {return(x[1] >= x[2])})
    if (!any(bool)) {return(list(c(), c()))}
    good <- which(bool)
    if (length(good) == 1) {return(list(good, good + 1))}

    #Define stretches of consecutive good bins
    diff <- diff(good)
    start <- good[c(1, which(diff != 1) + 1)]
    end <- good[c(which(diff != 1), length(good))] + 1
    return(list(start, end))
}

rho_chaos <- function(data, nlevels = 20, fillGaps = TRUE){
    data$rtiv <- (data$rtiv - min(data$rtiv)) /
                    (max(data$rtiv) - min(data$rtiv))
    ns <- sapply(seq_len(nlevels), function(n){
        t_n <- n/nlevels
        data$above <- data$rtiv > t_n
        nr <- nrow(data)
        #Filling 1-gaps
        if(fillGaps){
            for(i in seq_len(nr - 2)){
                if(data$above[i] & !data$above[i+1] & data$above[i+2]){
                    data$above[i+1] <- TRUE
                }
            }
        }

        bad <- which(!data$above)
        if(length(bad) == 0){return(1)}
        else{
            add <- 1
            if(length(bad) == 1){return(1+add)}
            return(length(which(diff(bad) > 1)) + add)
        }
    })
    rho <- 1 - sum(ns) / (nrow(data) * nlevels)
    return(rho)
}

#### filterSOI-related ####

#'@title filterSOI
#'@author Roger Gine
#'@description Performs a series of filters and quality checks to a given SOI
#'  list, removing unwanted SOIs in the process.
#'@inheritParams findSOI
#'@param id ID of the SOI list to be filtered/checked.
#'@param minint Minimun SOI intensity. All SOIs below this value will be removed
#'  from the SOI list
#'@param isofidelity Boolean. Whether to perform an isotopic fidelity check.
#'@param minscore Numeric. Minimum value (between 0 and 1) of isofidelity to
#'  retain an entry. Defaults to 0.8.
#'@return A filtered SOI list.
#' @examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#' struct <- filterSOI(struct, id = 1, minint = 10000, isofidelity = TRUE)
#'@export
setGeneric("filterSOI", function(struct, id, minint = 1e4, isofidelity, minscore = 0.8) {
    standardGeneric("filterSOI")
})
#' @rdname filterSOI
setMethod("filterSOI", signature = c("RHermesExp", "numeric", "ANY", "ANY", "ANY"),
    function(struct, id, minint = 1e4, isofidelity, minscore = 0.8) {
        soiobject <- struct@data@SOI[[id]]
        fname <- soiobject@filename
        PLid <- which(struct@metadata@filenames == fname)
        if(length(PLid) > 1){
            PLid <- PLid[1]
            warning("There's > 1 file with the same name, choosing the first")
        }
        PL <- struct@data@PL[[PLid]]
        ppm <- struct@metadata@ExpParam@ppm
        soilist <- soiobject@SOIList

        ##Filter by maximum intensity
        intense_enough <- which(soilist$MaxInt > minint)
        soilist <- soilist[intense_enough, ]
        if(nrow(soilist) == 0){
            warning("No SOIs remaining after filtering by intensity. ",
                    "Aborting filter")
            return(struct)
        }

        ##Filter by isotopic fidelity
        if (isofidelity) {
            # Isotopic elution similarity
            message("Computing isotopic elution similarity:")
            soilist <- isoCos(soilist, PL, isothr = 0.85)
            good <- which(!(soilist$MaxInt > 1e+06 & soilist$isofound == 0))
            soilist <- soilist[good, ]
            with_isos <- intense_enough[good]

            if(nrow(soilist) == 0){
                warning("No SOIs remaining after filtering by isotopic elution",
                        " similarity. Aborting filter")
                return(struct)
            }

            # Isotopic pattern similarity
            message("Calculating isotopic fidelity metrics:")
            isodata <- lapply(with_isos, plotFidelity, struct = struct,
                                id = id, plot = FALSE)

            cos <- vapply(isodata, function(x){x[[3]]}, numeric(1))
            soilist$isofidelity <- cos
            soilist <- soilist[cos > minscore, ]

            if(nrow(soilist) == 0){
                warning("No SOIs remaining after filtering by isotopic ",
                        "fidelity. Aborting filter")
                return(struct)
            }

            rtmargin <- 20
            # Removing confirmed isotopic signals
            soilist <- soilist[order(-soilist$MaxInt), ]
            message("Removing confirmed isotope entries:")
            toRemove <- numeric()
            for (i in seq_len(nrow(soilist))) {
                isomasses <- soilist[i, ]$isodf[[1]][, 2]
                st <- soilist[i, ]$start
                end <- soilist[i, ]$end
                #Entries in window
                idx <- which(between(soilist$start, st - rtmargin, end) &
                            between(soilist$end, st, end + rtmargin))
                entrymass <- soilist[idx, ]$mass
                overlaps <- lapply(isomasses, function(x) {
                    #Multiply ppm range to cover possible isotope mishaps
                    thr <- c(x - 3 * ppm * 1e-06 * x, x + 3 * ppm * 1e-06 * x)
                    if (any(entrymass > thr[1] & entrymass < thr[2])) {
                        return(which(entrymass > thr[1] & entrymass < thr[2]))
                    }
                    return()
                })
                overlaps <- do.call(rbind, overlaps)
                if (length(overlaps) != 0) {
                    idx <- idx[overlaps]
                    toRemove <- c(toRemove, idx)
                }
            }
            if (length(toRemove) != 0) {soilist <- soilist[-unique(toRemove)]}
        }

        ##Recalculate peaklist for plotting
        setkeyv(soilist, c("formula"))
        message("Recalculating peaklist for plotting:")
        plist <- lapply(unique(soilist$formula), recalculateDF, soilist)
        plist <- do.call(rbind, plist)
        plist$isov <- rep("M0", nrow(plist))

        ##Annotate adducts by cosine similarity
        soilist <- adCos(soilist, adthr = 0.8,
                                FATable = struct@metadata@ExpParam@ionF[[2]])
        struct@data@SOI[[id]]@SOIList <- soilist
        struct@data@SOI[[id]]@PlotDF <- as.data.table(plist)
        return(struct)
    })

recalculateDF <- function(i, soilist){
    data <- soilist[.(i), ]
    res <- data.frame(rt = numeric(), rtiv = numeric(),
                    form = character(), stringsAsFactors = FALSE)
    for (j in seq_len(nrow(data))) {
        peaks <- data[j, ]$peaks[[1]]
        if (is.null(dim(peaks)[1])) {next}
        peaks <- peaks[,1:2]
        form <- rep(unlist(data[j, ]$formula), nrow(peaks))
        peaks <- cbind(peaks, form)
        res <- rbind(res, peaks)
    }
    return(res)
}

isoCos <- function(soilist, PL, isothr = 0.99) {
    PL <- PL@peaklist
    setkeyv(PL, c("formv"))
    message("Calculating isotope similarity:")
    clist <- lapply(seq_len(nrow(soilist)), parallelIsoCos, soilist, PL,
                    isothr)
    hits <- lapply(clist, function(x) {return(x[[1]])})
    hitdf <- lapply(clist, function(x) {return(x[[2]])})
    soilist$isofound <- as.numeric(hits)
    soilist$isodf <- hitdf
    return(soilist)
}

adCos <- function(soilist, FATable, adthr = 0.8) {
    message("Calculating adduct similarity:")
    parseanot <- lapply(soilist$anot, function(x) {
        strsplit(x = x, split = "_")
    })
    soilist$f <- lapply(parseanot, function(x) {
        vapply(x, function(y) {y[[1]]}, character(1))
    })
    soilist$ad <- lapply(parseanot, function(x) {
        vapply(x, function(y) {y[[2]]}, character(1))
    })
    setkey(FATable, "f")
    clist <- lapply(seq_len(nrow(soilist)), parallelAdCos, soilist, FATable,
                    adthr)
    soilist$adrows <- clist
    return(soilist)
}

parallelIsoCos <- function(i, soilist, PL, isothr){
    SOI <- soilist[i, ]
    rti <- SOI$start[[1]]
    rtend <- SOI$end[[1]]
    cur <- PL[.(SOI$formula[[1]])]
    cur <- cur[cur$rt >= rti & cur$rt <= rtend, ]
    setkeyv(cur, c("isov", "rt"))
    count <- 0
    hitdf <- data.frame(iso = character(), mass = numeric(),
                    stringsAsFactors = FALSE)
    if (length(which(cur$isov == "M0")) < 5) {
        return(list(count, hitdf))
    }
    for (j in unique(cur$isov)) {
        if (j == "M0") {next}
        cos <- cosineSim(pattern = cur[which(cur$isov == "M0"), ],
                                query = cur[which(cur$isov == j),])
        if (!is.na(cos) & cos > isothr) {
        count <- count + 1
        hitdf <- rbind(hitdf, data.frame(iso = j,
                                        mass = mean(cur$mz[cur$isov == j]),
                                        stringsAsFactors = FALSE))
        }
    }
    return(list(count, hitdf))
}

parallelAdCos <- function(i, soilist, FATable, adthr){
    SOI <- soilist[i, ]
    rti <- SOI$start[[1]]
    rtend <- SOI$end[[1]]
    f <- SOI$f[[1]]
    equiv <- lapply(f, function(x) {FATable[.(x)]})
    equiv <- lapply(seq_along(f), function(x) {
        equiv[[x]][equiv[[x]]$ion != SOI$formula, ]
    })
    ids <- lapply(equiv, function(entries) {
        lapply(entries$ion, function(x) {
            candidates <- which(soilist$formula == x &
                                soilist$start >= (rti - 10) &
                                soilist$end <= (rtend + 10))
            if (length(candidates) == 0) {
                return()
            } else {
                return(lapply(candidates, function(row){
                    score <- cosineSim(pattern = SOI$peaks[[1]],
                                        query = soilist$peaks[row][[1]],
                                        nscans = 5)
                if (score > adthr) {return(row)}
                else {return()}
                }))
            }
        }) %>% unlist()
    })
    names(ids) <- unlist(SOI$ad)
    return(ids)
}


#### removeISF-related ####
#' @rawNamespace import(data.table, except = between)
#' @import magrittr
generateDiffDB <- function(DB, formulas, polarity = 1){
    metaesp <- DB$df_spectra
    metamet <- DB$df_metabolite
    espmet <- DB$df_spectraMetabolite
    list_fragments <- DB$list_fragments

    lapply(formulas, function(f){
        idmet <- metamet[.(f), ]$ID_metabolite
        if (is.na(idmet[1])) return(list())
        RES <- lapply(idmet, function(id) {
            idesp <- espmet[.(id), ]$ID_spectra
            metadata <- metaesp[.(idesp), ]

            #Filtering by low energy
            idesp <- espmet[.(id), ]$ID_spectra[metadata$REFCE %in%
                                                    c("0", "10","20","10eV")]
            metadata <- metadata[metadata$REFCE %in% c("0","10","20","10eV"), ]

            #Filtering by adduct
            idesp <- espmet[.(id), ]$ID_spectra[metadata$REFadduct %in%
                                                    c("M+H", "M-H", "M+")]
            metadata <- metadata[metadata$REFadduct %in% c("M+H","M-H", "M+"), ]
            which_polarity <- which(metadata$REFpolarity == polarity)
            spec <- list_fragments[.(idesp[which_polarity]),
                                "spectra"] %>% unlist(., recursive = FALSE)
            if (length(spec) == 0) {return(spec)}
            energies <- apply(metadata[which_polarity, c("REFprecursor_mz",
                                                        "REFCE", "REFnature")],
                                1, function(x){paste(x, collapse = "_")})
            names(spec) <- energies
            return(spec)
        })
        novalidspec <- vapply(RES, function(x) {length(x) == 0}, logical(1))
        if (any(novalidspec)) {
            RES <- RES[!novalidspec]
            if (length(RES) == 0) {
                return(list())
            }
            names(RES) <- paste(rep(f, times = nrow(metamet[.(f),]) -
                                length(which(novalidspec))),
                                metamet[.(f), ]$ID_metabolite[!novalidspec],
                                unlist(lapply(metamet[.(f), ]$REFname,
                                        function(x) {x[[1]]}))[!novalidspec],
                                metamet[.(f),]$REFsmiles[!novalidspec],
                                sep = "#")
        } else {
        names(RES) <- paste(rep(f, times = nrow(metamet[.(f),])),
                            metamet[.(f), ]$ID_metabolite,
                            lapply(metamet[.(f), ]$REFname,
                                    function(x) {x[[1]]}) %>% unlist(),
                            metamet[.(f), ]$REFsmiles, sep = "#")
        }
        common_deltas <- lapply(RES, function(structure){
            prec_m <- lapply(names(structure),
                            function(x){as.numeric(strsplit(x,
                                                            "_")[[1]][1])}) %>%
                            unlist(.)
            structure <- structure[!is.na(prec_m)]
            prec_m <- prec_m[!is.na(prec_m)]
            observed_mzs <- c()
            for (i in seq_along(structure)) {
                spectrum <- as.matrix(t(structure[[i]]))
                maxi <- max(spectrum[,2])
                spectrum <- spectrum[spectrum[,2] > 0.3*maxi, , drop = FALSE]
                deltas <- spectrum[,1] - prec_m[i]
                observed_mzs <- c(observed_mzs,
                                spectrum[abs(deltas) > 0.01 & deltas < 0, 1])
            }
            if (length(observed_mzs) == 0) {return()}
            return(sort(unique(observed_mzs)))
        })
    common_deltas <- common_deltas[vapply(common_deltas, length,
                            FUN.VALUE = numeric(1)) != 0]
    unlist(common_deltas) %>% unique() %>% sort()

    })
}

#' @rawNamespace import(data.table, except = between)
anotateISF <- function(SL, DB, polarity = 1){
    DB$df_spectra <- as.data.table(DB$df_spectra)
    DB$df_metabolite <- as.data.table(DB$df_metabolite)
    DB$df_spectraMetabolite <- as.data.table(DB$df_spectraMetabolite)
    DB$list_fragments <- data.table(ID_spectra = DB$list_fragments[[1]],
                                    spectra = DB$list_fragments[[2]])
    setkeyv(DB$df_spectra, c("ID_spectra"))
    setkeyv(DB$df_spectraMetabolite, c("ID_metabolite"))
    setkeyv(DB$df_metabolite, c("REFformula"))
    setkeyv(DB$list_fragments, c("ID_spectra"))

    SL$ISF <- lapply(seq_len(nrow(SL)), anotateParallelISF,
                        SL = SL, DB = DB, polarity = polarity)
    return(SL)
}

anotateParallelISF <- function(i, SL, DB, polarity){
    cur <- SL[i, ]
    masses <- generateDiffDB(DB, cur$f[[1]], polarity)
    masses <- masses[vapply(masses, length,
                            FUN.VALUE = numeric(1)) != 0]
    if (length(masses) == 0) {return(integer())}
    ISF <- lapply(masses, function(x){
        diffs <- do.call(cbind, lapply(x, function(mass){
        SL$mass - mass}))
    candidates <- which(apply(diffs, 1, function(m){any(abs(m) < 0.02)}))
    cosines <- vapply(candidates, function(cand){
        cosineSim(cur$peaks[[1]], SL$peaks[[cand]], nscans = 5)
    }, FUN.VALUE = numeric(1))
    candidates[cosines > 0.95 & candidates != i]
    })
    return(unlist(ISF))
}

#' @import igraph
#' @import visNetwork
#' @importFrom grDevices colorRamp rgb
#' @title plotISF
#' @family Plots
#' @author Roger Gine
#' @description Plots all annotated ISF relations between SOIs. NOTE: you have
#' to first run removeISF with the parameter justAnnotate set to TRUE, otherwise
#' the ISF SOIs are removed from the SOI list.
#' @inheritParams plotSOI
#' @examples
#' if(FALSE){
#'     plotISF(struct, 1)
#' }
#' @return A visNetwork plot
#' @export
plotISF <- function(struct, id){
    SOIlist <- SOI(struct, id)@SOIList
    net <- graph_from_adj_list(SOIlist$ISF)
    net <- set.vertex.attribute(net,
        "value",
        value = log10(SOIlist$MaxInt),
        index = which(seq_along(V(net)) %in% unique(SOIlist$originalID))
    )
    net <- set.vertex.attribute(net,
        "title",
        value = paste0(
            "<p>Mz: ", SOIlist$mass,"</p>",
            "<p>Intensity: ",
            round(SOIlist$MaxInt,digits = 2),
            "</p>", "<p>Anot: ",
            lapply(SOIlist$anot, function(an){
            paste(an, collapse = ",")}) %>%
            unlist(), "</p>"),
        index = which(seq_along(V(net)) %in% unique(SOIlist$originalID))
    )
    col <- colorRamp(c(rgb(1,0,0),rgb(0,1,0), rgb(0,0,1)), bias = 0.1)
    defaultcol <- rep("#777777", length(V(net)))
    node_color <- apply(col(SOIlist$mass/max(SOIlist$mass)), 1, function(x){
        x <- x/255
        rgb(x[1],x[2],x[3])
    })
    defaultcol[seq_along(V(net)) %in% unique(SOIlist$originalID)] <- node_color

    net <- set.vertex.attribute(net, "color", value = defaultcol)
    visnet <- toVisNetworkData(net)

    # distances <- apply(visnet[[2]], 1, function(x){
    #   SOIlist$mass[[x[2]]] - SOIlist$mass[[x[1]]]
    # })
    # edge_color <- apply(col(abs(distances)/max(abs(distances))),
    #   1, function(x){
    #   x <- x/255
    #   rgb(x[1],x[2],x[3])})
    #
    # visnet[[2]]$label <-  as.character(round(distances, 5))

    visNetwork(nodes = visnet[[1]], edges = visnet[[2]]) %>%
        visEdges(arrows = list(to = list(enabled = TRUE,
                                        scaleFactor = 0.5))) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(stabilization = FALSE)
}

#'@import igraph
cleanupISF <- function(SL){
    net <- igraph::graph_from_adj_list(SL$ISF)
    cl <- igraph::groups(igraph::cluster_walktrap(net))
    SL$originalID <- seq(nrow(SL))
    do.call(rbind, lapply(cl, function(group){
        subnet <- igraph::subgraph(net, group)
        df <- data.frame(out = degree(subnet, mode = "out"),
                        into = degree(subnet, mode = "in"))
        isf <- which(df$int > df$out)
        if (length(isf) == 0) {
            sois <- SL[group,]
            sois$group <- rep(list(group), nrow(sois))
        return(sois)
        }
        mInt <- max(SL$MaxInt[group[-isf]])
        isf <- isf[SL$MaxInt[group[isf]] < mInt]
        if (length(isf) != 0) {group <- group[-isf]}
        sois <- SL[group,]
        sois$group <- rep(list(group), nrow(sois))
        return(sois)
    }))
}

#' @title removeISF
#' @description Detect and remove ISF signals from a SOI list using low
#' collision energy MS2 data
#' @inheritParams filterSOI
#' @param DBpath Path to a MS2 spectral database with a specific format
#' @param justAnnotate Whether to actually remove the ISF or just keep them
#' annotated. Defaults to FALSE (removes the ISF). If set to TRUE, you can then
#' plot the ISF network with plotISF()
#' @return An updated RHermesExp object where the selected SOI has had its ISF
#' removed. If using justAnnotate = TRUE, the ISF aren't removed.
#' @seealso plotISF
#' @export
#' @examples
#' if(FALSE){
#'     #You need a MS2 Database with the proper format
#'     #(we are working on releasing one to the public soon enough)
#'     struct <- removeISF(struct, id = 1, DBpath = "MS2Database.rds")
#' }
removeISF <- function(struct, id, DBpath = "D:/MS2ID_B2R_20201113_083214.rds",
                        justAnnotate = FALSE){
    DB <- readRDS(DBpath)
    polarity <- ifelse(struct@metadata@ExpParam@ion == "+", 1, 0)

    #Annotate and remove ISF
    SoiObj <- struct@data@SOI[[id]]
    SL <- SoiObj@SOIList
    message("Anotating ISF:")
    SL <- anotateISF(SL, DB, polarity)
    message("Creating ISF network and cleaning:")
    SL <- cleanupISF(SL)
    SL <- as.data.table(SL)
    setkeyv(SL, "formula")

    #Recalculate plotting DF
    plist <- lapply(unique(SL$formula), preparePlottingDF,
                    SL)
    plist <- do.call(rbind, plist)
    plist$isov <- rep("M0", nrow(plist))

    #Update object and return
    SoiObj@SOIList <- SL
    SoiObj@PlotDF <- as.data.table(plist)
    struct@data@SOI[[id]] <- SoiObj
    struct <- setTime(struct, paste("Removed ISF from SOI list", id,
                                            "using database file", "DBpath"))
    return(struct)
}



