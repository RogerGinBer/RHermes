#'@export
#' @rdname RHermesPL-class
#' @param object An RHermesPL object
setMethod("show", "RHermesPL", function(object){
    message("Info about this PeakList:")
    message(paste("\tOriginal filename:", object@filename))
    message(paste("\tNumber of raw data points:", nrow(object@raw)))
    message(paste("\tNumber of PL points: ", nrow(object@peaklist)))
    message(paste("\tRedundance quotient: (Number of redundancies per",
                "unique entry) ",
                round(nrow(object@peaklist) / nrow(distinct(object@peaklist)),
                        digits = 3)))
    message(paste("\tNumber of scans", nrow(object@header)))
    rtrange <- range(object@header$retentionTime)
    message(paste0("\tAverage scan time: ",
                format((rtrange[2] - rtrange[1]) / nrow(object@header),
                        digits = 3), "s"))
    message(paste("\tLabelled processing:", object@labelled, "\n"))
})
