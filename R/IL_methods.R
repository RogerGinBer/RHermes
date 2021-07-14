#'@title filterIL
#'@description Filters out IL entries lower than a specified intensity in a
#' given RT interval.
#'@param struct The RHermesExp object.
#'@param id The IL ID in the RHermesExp object. The IDs are assigned by the
#' order in which the IL are generated.
#'@param rts Time interval to filter (two numeric values - start, end),
#' in seconds
#'@param minint Minimum entry intensity to be retained. All entries
#'<= minint will be removed in the specified rt interval. Defaults to
#'infinity, so all IL entries in the range are removed.
#'@return Nothing. As a side effect, it generates one/multiple .csv
#' files with the inclusion list data
#'@examples
#'if(FALSE){
#' filterIL(myHermes, 1, c(0,200), minint = 1e6)
#'}
#'@export
setGeneric("filterIL", function(struct, id, rts, minint = Inf) {
    standardGeneric("filterIL")
})
#' @rdname filterIL
setMethod("filterIL", c("RHermesExp", "numeric", "numeric", "ANY"),
function(struct, id, rts, minint = Inf) {
    if (length(rts) != 2) {
        stop("Please input just two RT values corresponding to the starting RT
                and ending RT you want to filter")
    }
    IL <- struct@data@MS2Exp[[id]]@IL@IL
    anot <- struct@data@MS2Exp[[id]]@IL@annotation
    which_to_remove <- which(IL$start >= rts[1] & IL$end <= rts[2] &
                                IL$MaxInt < minint)
    if (length(which_to_remove) != 0) {
        struct@data@MS2Exp[[id]]@IL@IL <- IL[-which_to_remove]
        struct@data@MS2Exp[[id]]@IL@annotation <- anot[-which_to_remove]
        struct <- setTime(struct, paste("IL entry", id, "filtered between",
                                        rts[1], "and", rts[2], "taking", minint,
                                        "as minimum intensity"))
    }
    return(struct)
})


#'@title exportIL
#'@md
#'@description Organizes the IL entries into multiple injections taking into
#'  account the user-specified parameters. Outputs a single or multiple csv
#'  files that serve as input for the MS to performed MSMS analysis.
#'
#'@param struct The RHermesExp object.
#'@param id The IL ID in the RHermesExp object. The IDs are assigned by the
#'  order in which the IL are generated.
#'@param file A string containing the folder to save the IL csv/s into and the
#'  basename for the file. By default will be your working directory and the
#'  default name is 'InclusionList' as in './InclusionList'.
#'@param mode Whether to plan set of continuous MS2 entries (PRM), adaptative
#'  injection time scans at the apexes of the entries' peaks or both (which
#'  optimizes instrument free time). Options are: "continuous", "deep" or
#'  "both".
#'@param maxOver Numeric, very important. It's the number of mz-rt segments that
#'  can be monitored at the same time by the MS instrument. Higher numbers lead
#'  to less injections but the number of scans for each IL entry will be reduced
#'  and gives problems when deconvoluting the MS2 spectras. It is ignored if
#'  mode = "deep". Defaults to 5.
#'@param defaultIT Numeric, the default IT in ms for continuous MS2 scans (only
#'  aplicable in Orbitrap instruments and for "continuous" and "both" modes).
#'  Defaults to 100ms.
#'@param maxInjections Numeric, the maximum number of planned injections to
#'  export. Defaults to 9999 to export all of them.
#'@param sepFiles Logical, whether to generate a single csv file or multiple
#'  csvs, each corresponding to each injection/chromatographic run. 
#'@return Nothing. As a side effect, it generates one/multiple .csv files with
#'  the inclusion list data
#'@examples
#'if(FALSE){
#'    exportIL(myHermes, 1, 'C:/SomeFolder', maxOver = 5, sepFiles = FALSE)
#'}
#'@export
setGeneric("exportIL", function(struct, id, file = "./InclusionList",
                                mode = "both", maxOver = 5, defaultIT = 100,
                                sepFiles = TRUE, maxInjections = 9999) {
    standardGeneric("exportIL")
})

#'@rdname exportIL
setMethod("exportIL", c("RHermesExp", "numeric", "ANY", "ANY", "ANY", "ANY"),
function(struct, id, file = "./InclusionList", mode = "both", maxOver = 5,
         defaultIT = 100, sepFiles = TRUE, maxInjections = 9999) {
    validObject(struct)
    if (length(struct@data@MS2Exp) == 0) {
        stop("This object doesn't have any ILs")
    }
    if (!between(id, 1, length(struct@data@MS2Exp))) {
        stop("Please enter a valid IL number")
    }
    IL <- struct@data@MS2Exp[[id]]@IL@IL
    IL$IT <- defaultIT
    plan <- injectionPlanner(IL, maxOver = maxOver, byMaxInt = TRUE,
                             injections = maxInjections,
                             mode = mode)
    if (sepFiles) {
    for (x in seq_along(plan)) {
        p <- plan[[x]]
        p$IT[p$IT > 1500] <- 1500
        p$f <- character(nrow(p))
        p$ad <- character(nrow(p))
        p$z <- 1
        p <- p[, c("f", "ad", "mass", "z", "start", "end",  "IT")]
        #Setting column style for Thermo Xcalibur import
        colnames(p) <- c("Formula", "Adduct", "m/z", "z", "t start (min)",
                         "t stop (min)", "Maximum Injection Time (ms)")
        p[, 5] <- floor(p[, 5]*100/60)/100
        p[, 6] <- ceiling(p[, 6]*100/60)/100
        # JC: need to consider max absolute RT
        #Added as result of Michi's comment
        p <- cbind(data.frame(Compound = seq_len(nrow(p))), p)
        write.csv(p, file = paste0(file, "_Injection_", x, ".csv"),
                  quote=F, row.names = FALSE)
    }
    } else {
        plandf <- do.call(rbind, lapply(seq_along(plan), function(x) {
            p <- plan[[x]]
            p$ID <- x
            return(p[,c("start", "end", "mass", "MaxInt", "ID")])
        }))
        write.csv(plandf, paste0(file, "_complete.csv"), row.names = FALSE)
    }
    return()
})


injectionPlanner <- function(IL, injections, maxOver, byMaxInt = TRUE,
                             mode = "continuous") {
    if (byMaxInt) {IL <- IL[order(-IL$MaxInt), ]}
    idx <- which(is.na(IL$start) | is.na(IL$end))  #NA depuration
    if (length(idx) != 0) {IL <- IL[-idx, ]}
    plan <- list()

    if(mode != "continuous"){
        message("Calculating high IT scans")
        deep_IL <- calculate_deep_IL(IL)
    } else {
        deep_IL <- data.frame()
    }

    mint <- min(IL$start)
    maxt <- max(IL$end)

    if(mode %in% c("continuous", "both")){
        while (nrow(IL) != 0 & injections > 0) {
            timeInt <- seq(mint,
                           maxt,
                           by = 0.005)
            OL <- rep(0, length(timeInt))
            ok_entries <- c()
            for (i in seq_len(nrow(IL))) {
                timeidx <- which(timeInt >= IL$start[i] & timeInt <= IL$end[i])
                if (any(OL[timeidx] >= maxOver)) {next}
                ok_entries <- c(ok_entries, i)
                OL[timeidx] <- OL[timeidx] + 1
            }
            curinj <- IL[ok_entries, ]
            IL <- IL[-ok_entries, ]

            if(mode == "both" & any(OL == 0)){
                deep_ok_entries <- c()
                for (i in seq_len(nrow(deep_IL))) {
                    timeidx <- which(timeInt >= deep_IL$start[i] &
                                         timeInt <= deep_IL$end[i])
                    if (any(OL[timeidx] >= 1)) {next}
                    deep_ok_entries <- c(deep_ok_entries, i)
                    OL[timeidx] <- OL[timeidx] + 1
                }
                if(length(deep_ok_entries) != 0){
                    deep_curinj <- deep_IL[deep_ok_entries, ]
                    deep_IL <- deep_IL[-deep_ok_entries, ]
                    curinj <- rbind(curinj, deep_curinj, fill = TRUE)
                }
            }

            plan <- c(plan, list(curinj))
            injections <- injections - 1
        }
        if(mode == "both" & nrow(deep_IL) != 0){
            message(paste(nrow(deep_IL), "high IT scans could not be planned",
                          "within the continuous MS2 injections.",
                          "Adding them separately after injection",
                          length(plan)))
            while (nrow(deep_IL) != 0 & injections > 0) {
                timeInt <- seq(min(deep_IL$start, na.rm = TRUE),
                               max(deep_IL$end, na.rm = TRUE),
                               by = 0.005)
                OL <- rep(0, length(timeInt))
                deep_ok_entries <- c()
                for (i in seq_len(nrow(deep_IL))) {
                    timeidx <- which(timeInt >= deep_IL$start[i] &
                                         timeInt <= deep_IL$end[i])
                    if (any(OL[timeidx] >= 1)) {next}
                    deep_ok_entries <- c(deep_ok_entries, i)
                    OL[timeidx] <- OL[timeidx] + 1
                }
                deep_curinj <- deep_IL[deep_ok_entries, ]
                deep_IL <- deep_IL[-deep_ok_entries, ]
                plan <- c(plan, list(deep_curinj))
                injections <- injections - 1
            }
        }
    } else {
        while (nrow(deep_IL) != 0 & injections > 0) {
            timeInt <- seq(min(deep_IL$start, na.rm = TRUE),
                           max(deep_IL$end, na.rm = TRUE),
                           by = 0.005)
            OL <- rep(0, length(timeInt))
            deep_ok_entries <- c()
            for (i in seq_len(nrow(deep_IL))) {
                timeidx <- which(timeInt >= deep_IL$start[i] &
                                     timeInt <= deep_IL$end[i])
                if (any(OL[timeidx] >= 1)) {next}
                deep_ok_entries <- c(deep_ok_entries, i)
                OL[timeidx] <- OL[timeidx] + 1
            }
            deep_curinj <- deep_IL[deep_ok_entries, ]
            deep_IL <- deep_IL[-deep_ok_entries, ]
            plan <- c(plan, list(deep_curinj))
            injections <- injections - 1
        }
    }


    if (nrow(IL) != 0) {
        message(paste0(nrow(IL), "ILs haven't been added to the injection plan",
                       " due to lack of space. Try again with more injections,",
                       " more maxOver or maxInjections"))
    }
    # print(plot(cumsum(sapply(plan, nrow))))
    # hist_data <- lapply(seq_along(plan), function(x){
    #     loc_plan <- plan[[x]]
    #     p <- hist(x = loc_plan$MaxInt %>% log10,
    #          breaks = seq(4,10,0.5))$counts
    #     # p <- p / nrow(loc_plan) * 100
    #     data.frame(sample = x, int = seq(4.25,10,0.5), p = p)
    # })
    # hist_data <- do.call(rbind, hist_data)
    # print(ggplotly(ggplot(hist_data) + geom_tile(aes(x=sample, y = int, fill = p))))
    return(plan)
}


#'@export
#'@rdname RHermesIL-class
#' @param object An RHermesIL object
setMethod("show", "RHermesIL", function(object){
    message("Info about the IL:")
    message(paste("\tIL entries:", nrow(object@IL)))
    message(paste("\tSOI index:", object@SOInum))
})

calculate_XIC_estimation <- function(raw, mzs, rts){
    points <- filter(raw, between(.data$mz, mzs[1], mzs[2]))
    points <- filter(points, between(.data$rt, rts[1], rts[2]))
    xic <- sapply(unique(points$rt), function(rt){
        sum(points$rtiv[points$rt == rt])
    })
    xic <- data.frame(rt = unique(points$rt), int = xic)
    return(xic)
}


calculate_deep_IL <- function(IL, intThr = 1e4){
    IL_deep <- IL %>%
        filter(.data$MaxInt > intThr)
    IL_deep <- lapply(seq_len(nrow(IL_deep)),function(entry){
        cur <- IL_deep[entry,]
        scans <- calculate_best_interval(cur$XIC[[1]], objective = 3e4)
        scans <- as.data.frame(scans)
        if(nrow(scans) == 0){return()}
        colnames(scans) <- c("start", "end")
        scans$mass <- cur$mass
        scans$MaxInt <- cur$MaxInt
        scans$entrynames <- cur$entrynames
        scans$XIC <- cur$XIC
        return(scans)
    })
    IL_deep <- do.call(rbind, IL_deep)
    IL_deep$IT <- (IL_deep$end - IL_deep$start) * 1000
    return(IL_deep)
}

calculate_apex <- function(scans){
    scans$rt[which.max(scans$int)]
}

calculate_best_interval <- function(scans, objective = 1e6, maxIT = 2000){
    if(nrow(scans) > 10){
        apexes <- scans$rt[inflect(scans$int, 4)]
        if(length(apexes) == 0){apexes <- calculate_apex(scans)}
    } else {
        apexes <- calculate_apex(scans)
    }
    intervals <- lapply(apexes, function(apex){
        maxt <- min(max(apex - min(scans$rt), max(scans$rt) - apex), maxIT/1000)
        scans <- approx(scans$rt, scans$int, n = 1000) %>% do.call(rbind, .) %>%
            t %>% as.data.frame
        for(t in seq(0.1, maxt, 0.01)){
            integral <- calculate_integral(filter(scans,
                                                  between(.data$x,
                                                          apex-t, apex+t)))
            if(integral > objective){return(c(apex-t, apex+t))}
        }
        return(c(apex-maxt, apex+maxt))
    })
    return(do.call(rbind, intervals))
}

calculate_integral <- function(scans){
    sum(diff(scans[[1]]) * (scans[[2]][-nrow(scans)] + diff(scans[[2]])/2))
}


# Authorship Evan Friedland, Stackoverflow #6836409
inflect <- function(x, threshold = 1){
    up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
    down <-  sapply(-1:-threshold, function(n){
        c(rep(NA, abs(n)),
            x[-seq(length(x), length(x) - abs(n) + 1)])
    })
    a    <- cbind(x,up,down)
    which(apply(a, 1, max) == a[,1])
}




