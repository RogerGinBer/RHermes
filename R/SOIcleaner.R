#' @export
setGeneric("SOIcleaner", function(struct, soiid, minint, isofidelity) {
    standardGeneric("SOIcleaner")
})
setMethod("SOIcleaner", signature = c("RHermesExp", "numeric",
    "ANY", "ANY"), function(struct, soiid, minint, isofidelity) {
    soiobject <- struct@data@SOI[[soiid]]
    fname <- soiobject@filename
    PLid <- which(struct@metadata@filenames == fname)
    PL <- struct@data@PL[[PLid]]
    ppm <- struct@metadata@ExpParam@ppm
    soilist <- soiobject@SoiList
    BiocParallelParam <- struct@metadata@cluster

    ##Filter by maximum intensity
    intense_enough <- which(soilist$MaxInt > minint)
    soilist <- soilist[intense_enough, ]

    ##Filter by isotopic fidelity
    if (isofidelity) {
        # Isotopic elution similarity
        message("Computing isotopic elution similarity:")
        soilist <- isoCos(soilist, PL, isothr = 0.85, BiocParallelParam)
        good <- which(!(soilist$MaxInt > 1e+06 & soilist$isofound == 0))
        soilist <- soilist[good, ]
        with_isos <- intense_enough[good]

        # Isotopic pattern similarity
        message("Calculating isotopic fidelity metrics:")

        # isodata <- bplapply(with_isos, IsoFidelity, struct = struct,
        #                     soilist = soiid, plot = FALSE,
        #                     BPPARAM = BiocParallelParam)
        
        isodata <- bplapply(with_isos, IsoFidelity, struct = struct,
                            soilist = soiid, plot = FALSE,
                            BPPARAM = SerialParam(progressbar = TRUE))

        cos <- vapply(isodata, function(x) {
            x[[3]]
        }, numeric(1))
        soilist <- soilist[cos > 0.5, ]

        rtmargin <- 20
        # Removing confirmed isotopic signals
        soilist <- soilist[order(-soilist$MaxInt), ]
        message("Removing confirmed isotope entries:")
        toRemove <- numeric()
        for (i in seq_len(nrow(soilist))) {
            isomasses <- soilist[i, ]$isodf[[1]][, 2]
            st <- soilist[i, ]$start
            end <- soilist[i, ]$end
            idx <- which(between(soilist$start, st - rtmargin, end) &
                        between(soilist$end, st, end + rtmargin))  #Entries in window
            entrymass <- soilist[idx, ]$mass
            overlaps <- lapply(isomasses, function(x) {
                #Multiply ppm range to cover possible isotope mishaps
                thr <- c(x - 3 * ppm * 1e-06 * x,
                         x + 3 * ppm * 1e-06 * x) 
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
    plist <- bplapply(unique(soilist$formula), recalculateDF, soilist,
                      BPPARAM = BiocParallelParam)
    plist <- do.call(rbind, plist)
    plist$isov <- rep("M0", nrow(plist))

    ##Annotate adducts by cosine similarity
    soilist <- RHermes:::adCos(soilist, FATable = struct@metadata@ExpParam@ionF[[2]],
        adthr = 0.8, BiocParallelParam = BiocParallelParam)
    struct@data@SOI[[soiid]]@SoiList <- soilist
    struct@data@SOI[[soiid]]@PlotDF <- as.data.table(plist)
    return(struct)
})

recalculateDF <- function(i, soilist){
    data <- soilist[.(i), ]
    res <- data.frame(rt = numeric(), rtiv = numeric(),
                      form = character(), stringsAsFactors = FALSE)
    for (j in seq_len(nrow(data))) {
        peaks <- data[j, ]$peaks[[1]]
        if (is.null(dim(peaks)[1])) {
            next
        }
        form <- rep(unlist(data[j, ]$formula), nrow(peaks))
        peaks <- cbind(peaks, form)
        res <- rbind(res, peaks)
    }
    return(res)
}

