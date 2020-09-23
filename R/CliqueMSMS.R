#' @title CliqueMSMS
#' @description Purify MSMS spectra detecting consistent mass signals
#' @details First it parses each IL entry annotation to store the annotated formulas
#' that correspond to each entry. Then, it reads the MSMS info and finds consistent
#' mass traces by sorting all data points by mz and using a Centwave-like algorithm.
#'  The derived 'pure' spectra are called superspectra. A single IL entry can yield different
#' superspectra. This function is NOT to be used by itself. It forms part of the MSMS processing
#' workflow.
#' @return List of formulas and a list of purified MSMS (list of dataframes,
#' one for each superspectra).
#' @import ggplot2
CliqueMSMS <- function(MS2Exp, idx, contaminant = 173.5, delta = 0.1,
                       BPPARAM = BiocParallel::SerialParam()) {
    IL <- MS2Exp@IL
    MS2list <- MS2Exp@MS2Data[idx]
    ## Parsing formulas from annotation for every IL entry
    fs <- strsplit(IL@IL$entrynames[idx], split = "#") %>% lapply(.,
        function(x) {
            lapply(x, function(y) {
                res <- strsplit(y, split = "$", fixed = TRUE)[[1]][[1]]
                return(sub(pattern = " ", replacement = "", x = res))
            }) %>% unlist() %>% unique()
        })

    #Main function    
    suppressWarnings(
        RES <- bplapply(seq_along(idx), generate_ss, BPPARAM = BPPARAM,
                        MS2list, contaminant, delta, fs, idx, FALSE)
    )
    RES <- do.call(rbind, RES)
    
    RES <- RHermes:::purify_ss(RES, minint = 30000, minpks = 2)
    
    RES$start <- as.numeric(RES$start)
    RES$end <- as.numeric(RES$end)
    RES$apex <- as.numeric(RES$apex)
    return(RES)
}

generate_ss <- function(curentry, MS2list, contaminant, delta, fs, idx,
                        to_plot){
    header <- MS2list[[curentry]][[1]]
    data <- MS2list[[curentry]][[2]]
    ## Remove weak signals (<0.5% of the most intense point)
    data <- do.call(rbind, lapply(unique(data$rt), function(t) {
        maxi <- max(data[data$rt == t, "int"])
        return(data[data$rt == t & data$int > 0.005 * maxi,
        ])
    }))
    if (nrow(data) == 0) {return()}
    ##Remove known contaminant signals
    for (x in contaminant) {
        data <- data[!between(data$mz, contaminant - delta,
                              contaminant + delta), ]
    }
    if (nrow(data) == 0) {return()}
    colnames(data)[colnames(data) == "int"] <- "rtiv"  #To match with cosineSim definition

    ##Find mass traces like CentWave
    rts <- unique(data$rt)
    soi <- list()
    avgmz <- c()
    mu <- 20
    pmin <- 5
    for (i in rts) {
        if (i == rts[1]) {
            soi <- split(data[data$rt == i, ], data$mz[data$rt == i])
            avgmz <- as.numeric(names(soi))
        } else {
            curdata <- data[data$rt == i, ]
            for (j in unique(curdata$mz)) {
                dist <- abs(avgmz - j)/j * 1e+06
                if (any(dist < mu)) {
                    id <- which.min(dist)
                    soi[[id]] <- rbind(soi[[id]], curdata[curdata$mz == j, ])
                    avgmz[[id]] <- mean(soi[[id]]$mz)
                } else {
                    soi <- c(soi, list(curdata[curdata$mz == j, ]))
                    avgmz <- c(avgmz, j)
                }
            }
        }
    }
    good <- vapply(soi, function(x) {nrow(x) > pmin}, logical(1))
    soi <- soi[good]
    avgmz <- avgmz[good]
    if (length(soi) == 0) return()
    # if (length(soi) == 1) {
    #     return(list(data.frame(mz = avgmz[[1]], int = max(soi[[1]]$rtiv))))
    # }
    ##Tidying the regions
    for (i in seq_along(avgmz)) {
        soi[[i]]$mz <- round(avgmz[i], 3)
        soi[[i]] <- soi[[i]][order(soi[[i]]$rt), ]
    }

    ##Centwave peak-picking
    soipks <- lapply(soi, function(x) {
        suppressWarnings(df <- xcms::peaksWithCentWave(int = x$rtiv,
                                                       rt = x$rt, snthresh = 0, firstBaselineCheck = FALSE,
                                                       prefilter = c(0, 0), peakwidth = c(5, 60)))
        df <- as.data.frame(df)
        if (nrow(df) == 0) {
            df <- data.frame(rt = 0, rtmin = min(x$rt), rtmax = max(x$rt),
                             into = 0, intb = 0, maxo = max(x$rtiv), sn = 0)
            df$categ <- "Putative"
        } else {
            df$categ <- "Peak"
            ##When a peak is found, consider the rest of fragment RT intervals as putative regions
            dfinterv <- unlist(df[, c("rtmin", "rtmax")])
            times <- c(min(x$rt), dfinterv[seq(1, length(dfinterv),2)],
                       dfinterv[seq(2, length(dfinterv), 2)],
                       max(x$rt))
            iter <- seq(1, length(times), 2)
            maxo <- mapply(function(t1, t2) {
                max(x$rtiv[between(x$rt, t1, t2)], na.rm = TRUE)
            }, times[iter], times[iter + 1]) %>% unlist()
            df <- rbind(df, data.frame(rt = 0, rtmin = times[iter],
                                       rtmax = times[iter + 1], into = 0, intb = 0,
                                       maxo = maxo, sn = 0, categ = "Putative"))
        }
        df$mz <- rep(x$mz[1], nrow(df))
        return(df)
    })
    pks <- do.call(rbind, soipks)
    if (nrow(pks) == 0) return()

    ##Matching data to the peaks found
    soi <- do.call(rbind, soi)
    soi$peak <- 0
    for (i in seq_len(nrow(pks))) {
        soi$peak[soi$mz == pks$mz[i] & between(soi$rt, pks$rtmin[i],
                                               pks$rtmax[i])] <- i
    }
    ##Deconvolution based on previous peak-picking -- Centwave-Propagation
    results <- RHermes:::centwavePropag(pks, soi)
    pks <- results[[1]]
    soi <- results[[2]]
    if (nrow(pks) == 0) return()

    ##Now calculate all similarities
    cos <- lapply(seq_len(nrow(pks)), function(x) {
        lapply(seq_len(nrow(pks)), function(y) {
            score <- RHermes:::pearsonSim(soi[soi$peak == x, ],
                                          soi[soi$peak == y, ])
            if (is.na(score)) {score <- 0}
            ifelse(score > 0.3, score, 0)
        })
    }) %>% unlist(.)
    cos <- matrix(cos, nrow = nrow(pks))

    net <- igraph::graph_from_adjacency_matrix(cos, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
    net <- igraph::simplify(net, remove.multiple = TRUE,
                            remove.loops = TRUE)
    
    ncomp <- igraph::components(net) %>% igraph::groups() %>% length()
    dens <- igraph::edge_density(net)
    if(is.nan(dens)){dens <- 0}
    
    if(ncomp == 1 & dens > 0.5){
        members <- rep(1, V(net) %>% length())
    } else {
        wc <- igraph::cluster_fast_greedy(net)
        members <- igraph::membership(wc)
    }


    n_mem <- length(unique(members))
    ss <- data.frame(start = numeric(n_mem), end = numeric(n_mem),
                     apex = numeric(n_mem), ILentry = numeric(n_mem),
                     precmass = numeric(n_mem))
    ss$start <- lapply(unique(members), function(id) {
        min(pks$rtmin[which(members == id)])
    })
    ss$end <- lapply(unique(members), function(id) {
        max(pks$rtmax[which(members == id)])
    })
    ss$apex <- lapply(unique(members), function(id) {
        apex <- which.max(soi$rtiv[soi$peak %in% which(members == id)])
        soi$rt[which(soi$peak %in% which(members == id))[apex]]
    })
    ss$ssdata <- lapply(unique(members), function(id) {
        return(data.frame(mz = pks$mz[which(members == id)],
                          int = pks$maxo[which(members == id)]))
    })
    ss$precmass <- rep(header$precursorMZ[1], nrow(ss))
    
    if(to_plot){
        return(list(net = net, members = members, data = data,
                    soi = soi, pks = pks, ss = ss))
    }
    
    ss$ILentry <- rep(idx[curentry], nrow(ss))
    ss$anot <- rep(fs[curentry], nrow(ss))

    return(ss)
}

reassign_and_check <- function(pks, soi) {
    if (nrow(soi) == 0) {
        return(list(data.frame(), soi))
    }
    ##Reassign scan points to peaks
    soi$peak <- 0
    for (i in seq_len(nrow(pks))) {
        soi$peak[soi$mz == pks$mz[i] & between(soi$rt, pks$rtmin[i],
            pks$rtmax[i])] <- i
    }

    ##Check that in each peak there are at least 5 scans
    tokeep <- vapply(seq_len(nrow(pks)), function(x) {
        length(which(soi$peak == x)) > 5
    }, logical(1))
    pks <- pks[which(tokeep), ]
    soi <- soi[soi$peak %in% which(tokeep), ]

    if (nrow(soi) == 0) {
        return(list(pks, soi))
    }

    #Reassign again
    soi$peak <- 0
    for (i in seq_len(nrow(pks))) {
        soi$peak[soi$mz == pks$mz[i] & between(soi$rt, pks$rtmin[i],
            pks$rtmax[i])] <- i
    }
    return(list(pks, soi))
}

centwavePropag <- function(pks, soi){
    repeat({
    rows_to_add <- data.frame(rt = numeric(0), rtmin = numeric(0),
                              rtmax = numeric(0), into = numeric(0),
                              intb = numeric(0), maxo = numeric(0),
                              sn = numeric(0), categ = character(0),
                              mz = numeric(0))
    for (i in which(pks$categ == "Peak")) {
        cur_cos <- vapply(seq_len(nrow(pks)), function(y) {
            score <- cosineSim(soi[soi$peak == i, ], soi[soi$peak == y, ])
            if (is.na(score)) {
                score <- 0
            }
            score
        }, numeric(1))
        #Putative peaks that match really well with a confirmed Centwave peak
        toconvert <- which(cur_cos > 0.8 & pks$categ == "Putative")
        if (length(toconvert) != 0) {
            for (j in toconvert) {
                dfinterv <- unlist(pks[i, c("rtmin", "rtmax")])
                x <- soi[soi$peak == j, ]
                #Define start-end pairs
                times <- c(min(x$rt), max(dfinterv[1], min(x$rt)),
                           min(dfinterv[2], max(x$rt)), max(x$rt))
                iter <- seq(1, length(times), 2)
                maxo <- mapply(function(t1, t2) {
                    max(x$rtiv[between(x$rt, t1, t2)], na.rm = TRUE)
                }, times[iter], times[iter + 1]) %>% unlist()
                pks[j, c("rtmin", "rtmax")] <- c(times[2], times[3])
                pks$maxo[j] <- max(x$rtiv[between(x$rt, times[2], times[3])])
                pks$categ[j] <- "Peak"
                rows_to_add <- rbind(rows_to_add,
                                        data.frame(rt = 0, rtmin = times[iter],
                                                   rtmax = times[iter + 1],
                                                   into = 0, intb = 0,
                                                   maxo = maxo, sn = 0,
                                                   categ = "Putative",
                                                   mz = pks$mz[j]))
            }
        }
    }
    rows_to_add <- rows_to_add[rows_to_add$maxo != -Inf, ]
    pks <- rbind(pks, rows_to_add)

    res <- RHermes:::reassign_and_check(pks, soi)
    pks <- res[[1]]
    soi <- res[[2]]

    ##Split putative peaks with >5s gaps
    for (i in which(pks$categ == "Putative")) {
        cur <- soi[soi$peak == i, c("rt", "rtiv")]
        times <- diff(cur$rt)
        jumps <- times > 3
        if (any(jumps)) {
            jumps <- which(jumps)
            splits <- data.frame(rt = numeric(0), rtmin = numeric(0),
                                 rtmax = numeric(0), into = numeric(0),
                                 intb = numeric(0), maxo = numeric(0),
                                 sn = numeric(0), categ = character(0),
                                 mz = numeric(0))
            for (j in seq_along(jumps)) {
                if (jumps[j] == jumps[1]) {
                    pks[i, c("rtmin", "rtmax")] <- c(cur$rt[1],
                                                     cur$rt[jumps[j]])
                    pks[i, "maxo"] <- max(cur$rtiv[seq_len(jumps[j])])
                } else {
                    splits <- rbind(splits,
                                    data.frame(rt = 0,
                                                rtmin = cur$rt[jumps[j-1]+1],
                                                rtmax = cur$rt[jumps[j]],
                                                into = 0, intb = 0,
                                                maxo = max(cur$rtiv[
                                                    seq_len(jumps[j])]),
                                                sn = 0, categ = "Putative",
                                               mz = pks$mz[i]))
                }
            }
            pks <- rbind(pks,
                         splits,
                         data.frame(rt = 0,
                                    rtmin = cur$rt[jumps[j] + 1],
                                    rtmax = max(cur$rt), into = 0, intb = 0,
                                    maxo = max(cur$rtiv[seq((jumps[j] + 1),
                                                            nrow(cur))]),
                                    sn = 0, categ = "Putative",
                                    mz = pks$mz[i]))
        }
    }
    res <- RHermes:::reassign_and_check(pks, soi)
    pks <- res[[1]]
    soi <- res[[2]]
    if (nrow(pks) == 0) {break}
    if (nrow(rows_to_add) == 0) {break}
    })
    return(list(pks, soi))
}

purify_ss <- function(sslist, minint = 30000, minpks = 2){
    good_ss <- vapply(sslist$ssdata, function(data){
        has_int <- any(data$int > minint)
        has_many_pks <- (nrow(data) > minpks)
        return(has_int | has_many_pks)
    }, logical(1))
    sslist <- sslist[good_ss, ]
    return(sslist)
}

