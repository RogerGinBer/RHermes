
#'@title SOIfinder
#'@description Generates a list of SOI (Scans of Interest) given a peaklist and appends it to the RHermesExp object.
#'
#'@param struct The RHermesExp object to which the SOI list will be saved.
#'@param params A SOIParam object or a list of SOIParam, contains all the parameters necessary to generate the SOI.
#'See \link[RHermes]{getSOIpar}
#'for more info on how to set it. If params is a list, SOIfinder will use each object for each SOI list it generates
#'until it has used them all. When that happens, the function will keep using the same SOIparam object for the rest of SOIs.
#'@param fileID Numeric vector of the indexes of the peaklists you want to process into a SOI.
#'@param against Optional. Numeric vector of the indexes of peaklists to use as a 'blank' against the correspoding fileID PL.
#'If you enter 0 (or don't input anything at all), no blank substraction will be performed. In case you want to process
#'multiple SOIs and only some of them have blank substraction, you must space the 0s so that the indexes match. (See examples)
#'
#'@examples
#'if(FALSE){
#'p <- getSOIpar('double') #Check getSOIpar documentation page to understand what this does
#'myHermes <- SOIfinder(myHermes, p, 1) #Basic usage, no blank substraction or anything
#'
#'#Blank substraction usage, 1 without blank and 2 using 1 as blank
#'myHermes <- SOIfinder(myHermes, p, c(1, 2), c(0, 1))
#'
#'#Even more advanced usage, with multiple parameters and blank substraction
#'#Here the first SOI is calculated using p and without blank sub.
#'#The second uses p2 and 1 as a blank
#'#And the same goes for the third. p2 would be reused and 1 as blank
#'p2 <- getSOIpar('triple')
#'myHermes <- SOIfinder(myHermes, c(p,p2), c(1, 2, 3), c(0, 1, 1))
#'}
#'
#' @return An updated RHermesExp object with the resulting SOI lists appended.
#' This allows for successive iterations of SOI generation, visualization and refining.
#'
#'@export
#'
#'@importFrom keras load_model_tf predict_classes
#'@importFrom dplyr distinct filter
#'@import magrittr
#'@import tensorflow
#'
#'

setGeneric("SOIfinder", function(struct, params, fileID, against = numeric(0)) {
    standardGeneric("SOIfinder")
})
setMethod("SOIfinder", c("RHermesExp", "ANY", "ANY", "ANY"),
    function(struct, params, fileID, against = numeric(0)) {
        if (length(against) == 0) {
            against <- rep(0, length(fileID))
            noBlSub <- TRUE
        }
        maxn <- length(struct@data@PL)
        if (any(c(fileID, against) > maxn)) {
            stop("Some of the indexes are above the number of current PL")
        }
        if (is.list(params)) {
            multParams <- TRUE
            specID <- c(seq_along(params), rep(length(params),
                times = length(fileID) - length(params)))
        } else {
            multParams <- FALSE
        }
        for (i in seq_along(fileID)) {
            idx <- fileID[i]
            if (idx == against[i]) {
                stop(paste0("You've tried to substract the blank from the same",
                  "file. That's not allowed."))
            }
            if (multParams) {
                cur <- params[[specID[i]]]
                ifelse(against[i] != 0, cur@blanksub <- TRUE,
                  cur@blanksub <- FALSE)
            } else {
                cur <- params
                ifelse(against[i] != 0, params@blanksub <- TRUE,
                  params@blanksub <- FALSE)
            }
            if (against[i] != 0) {
                cur@blanksub <- TRUE
                cur@blankname <- struct@data@PL[[against[i]]]@filename
                blankPL <- struct@data@PL[[against[i]]]@peaklist
            } else {
                cur@blanksub <- FALSE
                cur@blankname <- "None"
                blankPL <- NA
            }
            struct@data@SOI <- c(struct@data@SOI,
                                PLprocesser(struct@data@PL[[idx]],
                                    struct@metadata@ExpParam, cur, blankPL,
                                    struct@metadata@filenames[idx],
                                    struct@metadata@cluster)
                                )
            closeAllConnections()  #Closes parallel connections, just in case
            if (against[i] == 0) {
                struct <- setTime(struct, paste("Generated SOI list from",
                  struct@metadata@filenames[idx]))
            } else {
                struct <- setTime(struct, paste("Generated SOI list from",
                  struct@metadata@filenames[idx], "using",
                  struct@metadata@filenames[against[i]], "as blank"))
            }
        }
        return(struct)
    })



PLprocesser <- function(PL, ExpParam, SOIParam, blankPL = NA, filename,
                        BiocParallelParam = SerialParam()) {
    ## Extracting info from S4 objects into local variables ####
    DataPL <- PL@peaklist
    RAWdata <- PL@raw
    h <- PL@header
    ppm <- ExpParam@ppm
    formulaDB <- ExpParam@ionF[[1]]
    FA_to_ion <- ExpParam@ionF[[2]]

    params <- SOIParam@specs
    maxlen <- SOIParam@maxlen
    noise <- SOIParam@minint
    useblank <- SOIParam@blanksub

    ## Setting up PeakList ## ----------------------------------------------
    DataPL <- as.data.table(DataPL)
    DataPL <- DataPL[DataPL$isov == "M0", ]  #Filtering M0 beforehand to improve indexing performance later on
    DataPL <- DataPL[DataPL$rtiv > noise, ]
    setkeyv(DataPL, c("formv", "rt"))
    setkeyv(formulaDB, c("f"))
    setkeyv(FA_to_ion, c("ion"))

    ## Density filtering ## ------------------------------------------------
    message("Starting density filtering: ")
    GR <- apply(params, 1, function(x) {
        RHermes:::densityProc(x, DataPL, h, BiocParallelParam)
    })
    ## Grouping different filtered results ## -----------------------------
    message("Now Grouping:")
    Groups <- GR[[1]]
    setkeyv(Groups, c("formula", "start", "end"))
    if (length(GR) > 1) {
        for (i in 2:length(GR)) {
            Groups <- RHermes:::groupGrouper(GR, i, Groups, BiocParallelParam)
        }
    }
    Groups <- distinct(Groups)

    ## Initial Peak Retrieval ## ------------------------------------------
    message("Initial peak retrieval:")
    peakscol <- bplapply(seq_len(nrow(Groups)), RHermes:::retrievePeaks,
                         Groups, DataPL, BPPARAM = BiocParallelParam)
    Groups$peaks <- peakscol
    ## Group Characterization ## -------------------------------
    message("")
    message("Starting group characterization:")
    message("Mass calculation:")
    setkey(Groups, formula)
    Groups$mass <- formulaDB[.(Groups$formula),2]

    if (any(Groups$length > maxlen)) {
        Groups <- RHermes:::groupShort(Groups, maxlen, BiocParallelParam)
    }
    ## Blank substraction ## ----------------------------------
    if (useblank) {
        Groups <- RHermes:::blankSubstraction(Groups, blankPL,
                                              BiocParallelParam)
        if(nrow(Groups) == 0){return(RHermesSOI())}
    }
    ## Rest of characterization ## ---------------------------------------
    message("")
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

    message("Number of scans:")
    nscans <- apply(Groups, 1, function(x) {
        d <- x["peaks"][[1]]
        return(dim(d)[1])
    })
    Groups$nscans <- nscans

    message("Converting from ionic formula to F/A combinations:")
    Groups$anot <- lapply(Groups$formula, function(x) {
        apply(FA_to_ion[.(x), c(1, 3)], 1, function(y) {
            paste(y, collapse = "_")
        })
    })

    setkeyv(Groups, c("formula"))
    message("Generating peaklist for plotting:")
    plist <- bplapply(unique(Groups$formula), RHermes:::preparePlottingDF,
                      Groups, BPPARAM = BiocParallelParam)
    plist <- do.call(rbind, plist)
    plist$isov <- rep("M0", nrow(plist))

    # Constructing S4 output object
    output <- RHermesSOI(SoiList = Groups, PlotDF = as.data.table(plist),
        SoiParam = SOIParam, filename = filename)

    return(output)
}

densityProc <- function(x, DataPL, h, BiocParallelParam){
    rtbin <- x[1]
    scanspercent <- x[2]
    shift <- x[3]
    message("Running Density Filter")
    BinRes <- RHermes:::densityFilter(DataPL, h, rtbin, "M0", shift,
                                      BiocParallelParam)
    cutoff <- BinRes[[1]] * scanspercent
    #Correcting CUT < 1 cases to avoid errors (eg. if it was 0, any time region
    #would be included)
    cutoff[cutoff < 1] <- median(c(cutoff[cutoff > 1], 1)) #Added a 1 for robustness
    message("Running Density Interpreter")
    nwork <- BiocParallelParam$workers
    if(!is.numeric(nwork)) nwork <- 1

    uf <- unique(DataPL$formv)
    suppressWarnings({uf <- split(uf, seq_len(nwork))})
    id <- cumsum(vapply(c(0, uf), length, numeric(1))) - 1
    RES <- bplapply(seq_along(uf), RHermes:::parallelInterpreter, uf, cutoff,
                    BinRes, id, BPPARAM = BiocParallelParam)
    RES <- do.call(rbind, RES)
    RES$formula <- unlist(uf)[RES$formula]
    RES[, 1] <- (RES[, 1] * rtbin) - rtbin + floor(min(h$retentionTime)) +
        shift  #Changing from index to corresponding RT
    RES[, 2] <- (RES[, 2] * rtbin) - rtbin + floor(min(h$retentionTime)) +
        shift
    RES %<>% dplyr::mutate(., length = end - start)
    RES <- RES[, c(1, 2, 4, 3)]
    return(as.data.table(RES))
}

resolveGroup <- function(x, Groups, matched, first_df){
    corresp <- rbind(Groups[unlist(matched[.(x), 1]),
    ], first_df[x, ])
    if (nrow(corresp) == 1) {
        return(corresp)
    }
    corresp[1, 1] <- min(corresp[, 1])  #Min start
    corresp[1, 2] <- max(corresp[, 2])  #Max end
    return(corresp[1, ])
}

groupGrouper <- function(GR, i, Groups, BiocParallelParam){
    first_df <- GR[[i]]
    setkeyv(first_df, c("formula", "start", "end"))

    message(paste("Merging filter", i, "out of", length(GR),
                  "with previous SOI list"))
    ## Overlap between current group element and the next
    links <- foverlaps(Groups, first_df, by.x = c(4, 1, 2),
                       by.y = c(4, 1, 2), type = "any", which = TRUE)  #By formula, start and end time

    nomatched <- links[is.na(links$yid), ]$xid
    matched <- links[!is.na(links$yid), ]  #Which are found to be overlapping in both dfs
    setkeyv(matched, "yid")

    ## Grouping entries and choosing lowest start and highest end
    ## time
    res <- bplapply(unique(matched$yid), RHermes:::resolveGroup, Groups,
                    matched, first_df, BPPARAM = BiocParallelParam)
    res <- do.call(rbind, res)
    res[, `:=`(length, end - start)]
    res <- distinct(res)

    message("Adding entries unique from the first DF")
    res <- rbind(res, Groups[nomatched, ])

    message("Adding entries unique from the second DF")
    links <- foverlaps(first_df, Groups, by.x = c(4, 1, 2),
                        by.y = c(4, 1, 2), type = "any", which = TRUE)  #By formula, start and end time
    other_unmatched <- links[is.na(links$yid), ]$xid
    res <- rbind(res, first_df[other_unmatched, ])


    iter <- 1
    repeat {
        setkeyv(res, c("formula", "start", "end"))
        message(paste("Solving inner conflicts, Round:",
                      iter))
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
                                    which(table(links2$yid) ==1)]$yid,],
                     res2)
        iter <- iter + 1
        if(iter > 10){break} #For safety
    }
    return(res)
}

retrievePeaks <- function(i, Groups, PL){
    x <- Groups[i, ]
    rtmin <- x[1, 1]
    rtmax <- x[1, 2]
    f <- x[1, 4]
    pks <- PL[.(f)] %>% filter(., rt > rtmin[[1]] & rt < rtmax[[1]])
    return(pks[, c(1, 2)])
}

groupShort <- function(Groups, maxlen, BiocParallelParam){
    message("Shortening and selecting long groups:")
    SG <- filter(Groups, length <= maxlen)
    LG <- filter(Groups, length > maxlen)
    LG <- bplapply(seq_len(nrow(LG)), parallelGroupShort, LG,
                   maxlen, BPPARAM = BiocParallelParam)
    LG <- do.call(rbind, LG)
    return(rbind(SG, LG))
}

parallelGroupShort <- function(i, LG, maxlen){
    curGR <- LG[i,]
    #Divide long group into equal-sized smaller groups
    times <- seq(from = curGR[1, 1][[1]], to = curGR[1, 2][[1]],
                 length.out = ceiling(curGR[1, 3][[1]] / maxlen) + 1)
    deltat <- times[2]-times[1]

    NewGR <- data.table(start = times[-length(times)], end = times[-1],
                        length = deltat, formula = curGR$formula,
                        peaks = curGR$peaks,
                        mass = curGR$mass)

    #Reestructuration of the groups: data point 'splitting'
    NewGR[, "peaks"] <- apply(NewGR, 1, function(x) {
        pks <- x[5][[1]]
        return(pks[between(pks$rt, x[1][[1]], x[2][[1]]), ])
    })
    return(NewGR)
}

blankSubstraction <- function(Groups, blankPL, BiocParallelParam){
    message("Blank substraction:")
    setkeyv(blankPL, c("formv", "rt"))

    toKeep <- bplapply(seq_len(nrow(Groups)), firstCleaning, Groups, blankPL,
                  BPPARAM = BiocParallelParam) %>% unlist()
    sure <- Groups[which(toKeep), ]
    if(any(toKeep)){Groups <- Groups[-which(toKeep),]}
    reticulate::py_available(initialize = TRUE)
    if(reticulate::py_module_available("keras") &
       reticulate::py_module_available("tensorflow")){
        model <- load_model_tf(system.file("extdata", "model",
                                           package = "RHermes"))  #ANN
        setkeyv(blankPL, "formv")

        RES <- bplapply(seq_len(nrow(Groups)), prepareNetInput, Groups, blankPL,
                        BPPARAM = BiocParallelParam)

        Groups$MLdata <- RES
        NAgroups <- do.call(rbind, lapply(RES, function(x){is.na(x[1])}))

        if (any(NAgroups)) {
            Groups <- Groups[-which(NAgroups), ]
            RES <- RES[-which(NAgroups)]
        }
        organizeddata <- do.call(rbind, lapply(RES,
                                               function(x) {
                                                   return(c(x[1, ], x[2, ]))
                                               }))
        if(nrow(organizeddata) != 0){
            organizeddata <- keras::array_reshape(organizeddata,
                                                  c(nrow(organizeddata), 400),
                                                  order = "C")  #ANN input
            q <- model %>% keras::predict_classes(organizeddata)

            Groups <- Groups[-which(q == 0), ]  #ANN output
            Groups <- Groups[, -c("MLdata")]
        }
    } else {
        warning(paste0("A Keras installation was not found and blank",
                        "substraction was not performed"))
    }
    Groups <- rbind(sure, Groups)
}

firstCleaning <- function(i, Groups, blankPL){
  cur <- Groups[i, ]
  st <- cur$start
  end <- cur$end
  f <- cur$formula
  deltat <- 10
  peaks <- cur$peaks[[1]]
  blankpks <- blankPL[.(f)] %>% filter(., rt >= st - deltat &
                                         rt <= end + deltat &
                                         isov == "M0")
  blankpks <- distinct(blankpks[, c(1, 2)])
  if(nrow(blankpks) < 5){return(TRUE)} #No blank signals
  
  sampleCV <- sd(peaks$rtiv)/mean(peaks$rtiv)
  blankCV <- sd(blankpks$rtiv)/mean(blankpks$rtiv)
  sampleMax <- max(peaks$rtiv)
  blankMax <- max(blankpks$rtiv)
  
  #We have to be a bit stringent with the conditions, otherwise we collect junk
  if(sampleCV/blankCV > 5){return(TRUE)} 
  if(sampleMax/blankMax > 3 & sampleMax > 30000){return(TRUE)}
  
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
        smooth_pks <- data.frame(approx(x = peaks,
                                        xout = seq(from = st, to = end,
                                                   length.out = Npoints),
                                        rule = 1)) #Interpolate to N points
        smooth_pks[is.na(smooth_pks[, "y"]), "y"] <- 0
        blankpks <- blankPL[.(f)] %>% filter(., rt >= st - deltat &
                                                 rt <= end + deltat &
                                                 isov == "M0")
        blankpks <- distinct(blankpks[, c(1, 2)])
        if (length(unique(blankpks$rt)) < 2) {
            blankpks <- data.frame(rt = c(st, end), rtiv = c(0, 0))
            smooth_blankpks <- data.frame(
                approx(x = blankpks,
                        xout = seq(from = st, to = end, length.out = Npoints),
                        rule = 1))
        } else {
            smooth_blankpks <- data.frame(
                approx(x = blankpks,
                        xout = seq(from = st, to = end, length.out = Npoints),
                        rule = 1))
            smooth_blankpks[is.na(smooth_blankpks[, "y"]), "y"] <- 0
        }
        return(as.matrix(rbind(smooth_pks[, 2], smooth_blankpks[,2]))/
                   max(smooth_pks[, 2]))
    }, error = function(cond) {
        return(NA)
    })

}

preparePlottingDF <- function(i, Groups){
    data <- Groups[.(i), ]
    res <- data.frame(rt = 0.1, rtiv = 0.1, form = "",
                      stringsAsFactors = FALSE)
    for (j in seq_len(nrow(data))) {
        peaks <- data[j, ]$peaks[[1]]
        if (is.null(dim(peaks)[1])) {
            next
        }
        form <- rep(unlist(data[j, ]$formula), nrow(peaks))
        peaks <- cbind(peaks, form)
        res <- rbind(res, peaks)
    }
    res <- res[-1, ]
    return(res)
}
