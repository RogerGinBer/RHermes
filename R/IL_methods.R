injectionPlanner <- function(IL, injections, maxover, byMaxInt = TRUE,
                                returnAll = FALSE) {
    if (returnAll) {injections <- 1e4}
    if (byMaxInt) {IL <- IL[order(-IL$MaxInt), ]}
    idx <- which(is.na(IL$start) | is.na(IL$end))  #NA depuration
    if (length(idx) != 0) {IL <- IL[-idx, ]}
    plan <- list()
    while (nrow(IL) != 0 & injections > 0) {
        TimeInt <- seq(min(IL$start, na.rm = TRUE), max(IL$end,
                                                        na.rm = TRUE), by = 0.5)
        OL <- rep(0, length(TimeInt))
        ok_entries <- c()
        for (i in seq_len(nrow(IL))) {
            timeidx <- which(TimeInt >= IL$start[i] & TimeInt <= IL$end[i])
            if (any(OL[timeidx] >= maxover)) {next}
            ok_entries <- c(ok_entries, i)
            OL[timeidx] <- OL[timeidx] + 1
        }
        curinj <- IL[ok_entries, ]
        IL <- IL[-ok_entries, ]
        plan <- c(plan, list(curinj))
        injections <- injections - 1
    }
    if (nrow(IL) != 0) {
        warning(paste0("Some SOI haven't been added to the injection plan due",
                        " to lack of space. Try again with more injections,",
                        "  more maxover or returnAll = TRUE"))
    }
    return(plan)
}

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
#'@description Organizes the IL entries into multiple injections
#'  taking into account the user-specified parameters. Outputs a
#'  single or multiple csv files that serve as input for the MS
#'  to performed MSMS analysis.
#'
#'@param struct The RHermesExp object.
#'@param id The IL ID in the RHermesExp object. The IDs are
#'  assigned by the order in which the IL are generated.
#'@param folder A string containing the folder to save the IL
#'  csv/s into. By default will be your working directory
#'@param maxOver Numeric, very important. It's the number of
#'  mz-rt segments that can be monitored at the same time by the
#'  MS instrument. Higher numbers lead to less injections but the
#'  number of scans for each IL entry will be reduced and gives
#'  problems when deconvoluting the MS2 spectras.
#'@param sepFiles Logical, whether to generate a single csv file
#'  or multiple csvs, each corresponding to each
#'  injection/chromatographic run. From our experience with an
#'  Orbitrap Fusion, separate csvs will simplify the task.
#'@return Nothing. As a side effect, it generates one/multiple
#'  .csv files with the inclusion list data
#'@examples
#'if(FALSE){
#'    exportIL(myHermes, 1, 'C:/SomeFolder', maxOver = 5, sepFiles = FALSE)
#'}
#'@export
setGeneric("exportIL", function(struct, id, file = "./InclusionList",
                                maxOver = 5, sepFiles = FALSE) {
    standardGeneric("exportIL")
})

#'@rdname exportIL
setMethod("exportIL", c("RHermesExp", "numeric", "ANY", "ANY", "ANY"),
function(struct, id, file = "./InclusionList", maxOver = 5, sepFiles = FALSE) {
    validObject(struct)
    if (length(struct@data@MS2Exp) == 0) {
        stop("This object doesn't have any ILs")
    }
    if (!between(id, 1, length(struct@data@MS2Exp))) {
        stop("Please enter a valid IL number")
    }
    IL <- struct@data@MS2Exp[[id]]@IL@IL
    plan <- injectionPlanner(IL, 10, maxOver, byMaxInt = TRUE, returnAll = TRUE)
    if (sepFiles) {
    for (x in seq_along(plan)) {
        p <- plan[[x]]
        p <- p[, c("start", "end", "mass")]
        #Setting column style for Thermo Xcalibur import
        colnames(p) <- c("t start (min)", "t stop (min)", "m/z")
        p[, 1] <- p[, 1]/60
        p[, 2] <- p[, 2]/60
        #Added as result of Michi's comment
        p <- cbind(data.frame(Compound = seq_len(nrow(p))), p)
        write.csv(p, file = paste0(file, "_Injection_", x, ".csv"),
                    row.names = FALSE)
    }
    } else {
        plandf <- do.call(rbind, lapply(seq_along(plan), function(x) {
            p <- plan[[x]]
            p$ID <- x
            return(p)
        }))
        write.csv(plandf, paste0(file, "_complete.csv"))
    }
    return()
})


#'@export
#'@rdname RHermesIL-class
#' @param object An RHermesIL object
setMethod("show", "RHermesIL", function(object){
    message("Info about the IL:")
    message(paste("\tIL entries:", nrow(object@IL)))
    message(paste("\tSOI index:", object@SOInum))
})
