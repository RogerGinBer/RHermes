#'@import data.table
parallelFilter <- function(j, ScanResults, bins, timebin){
    lapply(j, function(i) {
        data <- as.vector(ScanResults[.(i), rt])
        mint <- min(data)
        maxt <- max(data)
        res <- rep(0, length(bins))
        goodbins <- bins[between(bins, mint-timebin, maxt+timebin)]
        if (length(goodbins) == 0) {
            return(res)
        }
        l <- mapply(function(b1, b2) {
            return(length(which(data > b1 & data < b2)))
        }, goodbins[-length(goodbins)], goodbins[-1])  #Number of entries within the time bin
        if (length(l) == 0) {
            return(res)
        }
        st <- which(bins == goodbins[1])
        res[st:(st + length(goodbins) - 2)] <- l
        return(res)
    })
}


densityFilter <- function(ScanResults, h, timebin, iso = "M0", tshift = 0,
                          BiocParallelParam) {
    #Setting time bins-------------------------------------------------
    rtmin <- floor(min(h$retentionTime))
    rtmax <- ceiling(max(h$retentionTime))
    bins <- seq(from = rtmin + tshift, to = rtmax + tshift, by = timebin)

    #Getting how many scans were taken on each bin---------------------
    scans <- lapply(bins, function(curbin) {
        return(h %>% filter(., retentionTime > curbin & retentionTime <=
            curbin + timebin) %>% dim(.) %>% .[1])
    })

    #Counting how many scan entries are on each bin--------------------
    ## Parallelized approach
    nwork <- BiocParallelParam$workers
    if(!is.numeric(nwork)) nwork <- 1

    idx <- split(unique(ScanResults$formv), seq_len(nwork) * 5)
    setkeyv(ScanResults, c("formv", "rt"))
    RES <- bplapply(idx, parallelFilter, ScanResults, bins, timebin,
             BPPARAM = BiocParallelParam)

    RES <- unlist(RES, recursive = FALSE)
    #Adding bins used and scan numbers for future analysis
    #reference (eg: DensityInterpreter) -------------------------------
    names(RES) <- unlist(idx)

    RES <- c(list(unlist(scans)), list(bins), RES)
    names(RES)[c(1,2)] <- c("Bins","Scans")

    return(RES)
}

