OptScanSearch <- function(DB, raw, mzList, ppm, IsoList, labelled = FALSE,
    minhit = 0, BiocParallelParam) {
    message(paste0("This process can take quite a bit of time, depending on ",
                   "the processing power and RAM your computer has"))

    setkey(raw, mz)
    DB <- as.data.table(DB)  #Making sure its a DT
    setkey(DB, f)

    if(labelled){
        DB$numC <- as.numeric(lapply(DB$envi, function(x){
            elements <- CHNOSZ::count.elements(x)
            if("C" %in% names(elements)){return(elements[["C"]])}
            else{return(0)}
        }))
    }

    ncores <- ifelse(is.numeric(BiocParallelParam$workers[[1]]),
                     yes = BiocParallelParam$workers, no = 1)
    #Splitting the formulas into a list (with l = number of workers) to reduce
    #time loss associated with variable loading (in SOCK only)
    flist <- split(DB, f = seq_len(ncores))
    PLresults <- bplapply(seq_along(flist), PLparallelfun, flist, raw, IsoList,
                          labelled, ppm, minhit, BPPARAM = BiocParallelParam)
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
            RHermes:::regularProc(curDB, mass, formula, pmz, curiso, ppm,
                                  IsoList, minhit, i, raw)
        } else {
            numC <- curDB$numC[i]
            RHermes:::labelledProc(curDB, mass, formula, pmz, curiso, ppm,
                                   IsoList, minhit, i, raw, numC)
        }
    })
    return(do.call(rbind, localRES))
}

regularProc <- function(curDB, mass, formula, pmz, curiso, ppm, IsoList, minhit,
                        i, raw){
    ss <- lapply(mass, function(m) {
        if (m < pmz[1, 1] | m > pmz[nrow(pmz), 1]) {
            return()
        }
        RHermes:::binarySearch(pmz, m, ppm)
    })
    l <- vapply(ss, length, numeric(1))

    adidx <- which(l > minhit)
    if (length(adidx) == 0L)
    {
        return()
    }  #Return nothing if no M0 hit

    isofactors <- curiso[[i]]  #If hit, let's find the isotopologues
    # if(any(is.na(isofactors))){return()}
    isodf <- IsoList[[2]][isofactors, ]
    ch <- curDB[i, 3]  #Charge to normalize isotope deltam's


    output <- lapply(adidx, function(j) {
        scanid <- raw[ss[[j]]]
        scanid$formv <- formula
        scanid$isov <- "M0"

        isom <- mass + (isodf$deltam/abs(ch[[1]]))
        # isom <- isom[!is.na(isom)]

        isoss <- lapply(isom, function(m) {
            if (m < pmz[1, 1] | m > pmz[nrow(pmz), 1]) {
                return()
            }
            RHermes:::binarySearch(pmz, m, ppm * 1.5)
            #Added small multiplicative factor to ppm. We've seen that
            #isotope peaks may have a bit more error than M0
        })
        isol <- vapply(isoss, length, FUN.VALUE = numeric(1))

        isoidx <- which(isol != 0)
        if (length(isoidx) != 0) {
            isoentries <- do.call(rbind, lapply(isoidx,
                                                function(x) {
                                                    isoid <- raw[isoss[[x]]]
                                                    isoid$formv <- formula
                                                    isoid$isov <- as.character(isodf$ID[x])
                                                    return(isoid)
                                                }))
            isoentries <- isoentries[isoentries$rt %in%
                                         unique(scanid$rt), ]
            scanid <- rbind(scanid, isoentries)
        }
        return(scanid)
    })
    return(do.call(rbind, output))
}

labelledProc <- function(curDB, mass, formula, pmz, curiso, ppm,
                         IsoList, minhit, i, raw, numC){

    if (mass < pmz[1, 1] | mass > pmz[nrow(pmz), 1]) {
        ss <- numeric()
    } else {
        ss <- RHermes:::binarySearch(pmz, mass, ppm)
    }
    if(length(ss) == 0){return()}
    
    ##Note that we don't need any M0 hit to find the other signals
    isofactors <- curiso[[i]]
    isodf <- IsoList[[2]][isofactors, ]
    if(numC != 0){
      isodf <- rbind(isodf, data.frame(ID = paste0("M", seq_len(numC)),
                                       deltam = 1.003355*seq_len(numC)))
      isodf <- dplyr::distinct(isodf)
    }
    ch <- curDB[i, 3]  #Charge to normalize isotope deltam's


    scanid <- raw[ss, ]
    if(nrow(scanid) != 0){
        scanid$formv <- formula
        scanid$isov <- "M0"
    }

    isom <- mass + (isodf$deltam/abs(ch[[1]]))
    isoss <- lapply(isom, function(m) {
        if (m < pmz[1, 1] | m > pmz[nrow(pmz), 1]) {
            return()
        }
        RHermes:::binarySearch(pmz, m, ppm * 1.5)
        #Added small multiplicative factor to ppm. We've seen that
        #isotope peaks may have a bit more error than M0
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
                                }))
        if(nrow(scanid) == 0){return(isoentries)}
        return(rbind(scanid, isoentries))
    }
    if(nrow(scanid) == 0){return()}
    return(scanid)
}

binarySearch <- function(plist, m, ppm) {
    i <- round((m - plist[1, 1])/(plist[nrow(plist), 1] - plist[1,1]) *
                   nrow(plist))[[1]]  #Biased start based on target
    if (i == 0) {
        i <- 1
    }
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
                  if (top + j > ub) {
                    break
                  }
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
    if (!out) {
        return()
    }
    return(bot:top)
}

