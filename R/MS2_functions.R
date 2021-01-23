#### Main MS2Proc routine ####

#' @import tidyverse
#' @import igraph
#' @import methods
#' @import plotly
#' @importFrom dplyr distinct
#'
#'@title MS2Proc
#'@description MS2Proc processes the MSMS files generated for a given IL and performs
#'compound identification using an MSMS reference spectra database. Only the struct
#'id and MS2files parameters are mandatory.
#'
#'@param struct The RHermesExp object to add the processed data.
#'@param id Numeric, the id of the inclusion list to process.
#'@param MS2files Character vector of the MS2 files address.
#'@param referenceDB Character, where to find the MSMS database to query against.
#'@param mincos Numeric, between 0 and 1, the minimum spectral cosine between query and
#'reference spectrum for a match to be reported (see the paper for more detail on this).
#'@return An RHermesExp object with the identifications set in the used MS2Exp slot.
#'@examples
#'MS2Proc(myHermes, 1, c('C:/myFolder/File1.mzML', 'C:/myFolder/File2.mzML'))
#'
#' @export
setGeneric("MS2Proc", function(struct, id, MS2files, sstype = "regular", useDB = FALSE,
                               referenceDB = "D:/sp_MassBankEU_20200316_203615.RData", mincos = 0.8) {
  standardGeneric("MS2Proc")
})
setMethod("MS2Proc", c("RHermesExp", "numeric", "character", "character",
                       "ANY", "ANY"), function(struct, id, MS2files, sstype = "regular", useDB = FALSE,
                                               referenceDB = "D:/sp_MassBankEU_20200316_203615.RData", mincos = 0.8) {
 ## MS2 Data Importation and Sorting within IL 
 validObject(struct)
 MS2Exp <- struct@data@MS2Exp[[id]]
 message("Starting MS/MS data importation, merging and sorting within the IL entries")
 
 MS2Exp@MS2Data <- MSMSimporter(MS2Exp@IL, MS2files)  #Fills in the MS2data slot
 idx <- vapply(MS2Exp@MS2Data, function(x) {
   return(length(x) != 0)
 }, logical(1))  #Which IL entries are covered by at least 1 scan
 idx <- which(idx)
 if (length(idx) == length(MS2Exp@MS2Data)) {
   message("All IL entries were covered. Nice!")
 } else {
   message("A total of ", length(MS2Exp@MS2Data) - length(idx)," (",
           round((length(MS2Exp@MS2Data) - length(idx))/
                   length(MS2Exp@MS2Data)*100, digits = 2),
           "%) entries were not covered in the experiment")
 }
 
 if (length(idx) == 0) {
   stop("No entries were covered in the IL. Please check the MS2 file paths and the IL object are correct")
 }
 
 ##Superspectra Generation 
 message("Starting superspectra generation. This may take a while...")
 purifiedSpectra <- RHermes:::CliqueMSMS(MS2Exp, idx, sstype = sstype)
 
 if(!useDB){
   message("No spectral matching was performed. Done!")
   purifiedSpectra$results <- rep("No DB matching",
                                  times = nrow(purifiedSpectra))
   MS2Exp@Ident <- list(MS2Features = purifiedSpectra,
                        DatabaseSpectra = list(),
                        MS2_correspondance = list())
   struct@data@MS2Exp[[id]] <- MS2Exp
   struct <- setTime(struct, paste("Processed MS2Exp object", id,
                                   "and generated superspectra but didn't perform",
                                   "matching"))
   return(struct)
 }
 
 ##Database Query
 # Catching possible load error
 tryCatch({
   obj <- readRDS(referenceDB)
   # Data.tables for faster indexing
   metaesp <- as.data.table(obj$df_spectra)
   metamet <- as.data.table(obj$df_metabolite)
   espmet <- as.data.table(obj$df_spectraMetabolite)
   list_fragments <- obj$list_fragments
   rm(obj)
 }, error = function(cond) {
   stop("The referenceDB input isn't valid")
 })
 
 metamet$REFformula <- unlist(metamet$REFformula)
 
 setkeyv(metaesp, c("ID_spectra"))
 setkeyv(espmet, c("ID_metabolite"))
 setkeyv(metamet, c("REFformula"))
 
 polarity <- ifelse(struct@metadata@ExpParam@ion == "+", 1, 0)
 
 # Retrieving spectra
 message("Retrieving MS2 spectra from the reference database")
 allf <- purifiedSpectra$anot %>% unlist() %>% unique()
 retrievedMSMS <- lapply(allf, function(f) {
   idmet <- metamet[.(f), ]$ID_metabolite
   if (is.na(idmet[1])) return(list())
   RES <- lapply(idmet, function(id) {
     idesp <- espmet[.(id), ]$ID_spectra
     metadata <- metaesp[.(idesp), ]
     which_polarity <- metadata$REFpolarity == polarity
     spec <- list_fragments[["spectra"]][list_fragments[["ID_spectra"]] %in%
                                           idesp[which_polarity]]
     if(length(spec)==0){return(spec)}
     energies <- apply(metadata[which_polarity, c("REFCE", "REFnature")], 1,
                       function(x){paste(x, collapse = "_")})
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
                         metamet[.(f),]$REFsmiles[!novalidspec], sep = "#")
   } else {
     names(RES) <- paste(rep(f, times = nrow(metamet[.(f),])),
                         metamet[.(f), ]$ID_metabolite,
                         lapply(metamet[.(f), ]$REFname,
                                function(x) {x[[1]]}
                         ) %>% unlist(),
                         metamet[.(f), ]$REFsmiles, sep = "#")
   }
   return(RES)
 })
 names(retrievedMSMS) <- allf
 
 message("Calculating Cosine similarities")
 
 corresponding <- lapply(purifiedSpectra$anot, function(x){
   which(allf %in% x)
 })
 cos_list <- mapply(function(entry, reference) {
   curspec <- retrievedMSMS[reference]
   n <- unlist(lapply(curspec, names))
   curspec <- unlist(curspec, recursive = FALSE, use.names = FALSE)
   names(curspec) <- n
   cos <- rapply(curspec, function(compound) {
     MSMScosineSim(entry, t(compound), minhits = 1)
   }, classes = "matrix", how = "replace")
   cos[which(lapply(cos, length) != 0)]
 }, purifiedSpectra$ssdata, corresponding)
 
 purifiedSpectra$results <- lapply(cos_list, function(x) {
   if (length(x) == 0) {return("Missing reference spectra")}
   isValidhit <- any(rapply(x, function(cos) {cos > mincos},
                            "numeric", "unlist"))
   if (!isValidhit) {return("No significant hits")}
   
   # Extracting cosines for each superspectra
   withcosines <- which(vapply(x, function(x) {
     length(x) != 0
   }, logical(1)))
   
   #Formula entries that have a good match
   goodf <- which(vapply(x, function(f){
     any(unlist(f) > mincos)
   }, logical(1)))
   
   RES <- lapply(goodf, function(f) {
     #Detect which spectra have a hit for that formula
     fnames <- names(x[f])
     resdf <- data.frame(formula = fnames, cos = numeric(1),
                         id = numeric(1),
                         stringsAsFactors = FALSE)
     
     coslist <- x[[f]]
     resdf$id <- which.max(vapply(coslist, function(cos) {
       cos[1]
     }, numeric(1)))
     resdf$cos <- coslist[[resdf$id]]
     
     return(resdf)
   }) %>% do.call(rbind, .)
   return(RES)
 })
 MS2Exp@Ident <- list(MS2Features = purifiedSpectra,
                      DatabaseSpectra = retrievedMSMS,
                      MS2_correspondance = corresponding)
 message("Done!")
 struct@data@MS2Exp[[id]] <- MS2Exp
 struct <- setTime(struct, paste("Processed MS2Exp object", id,
                                 ", generated superspectra and matched against",
                                 "database", referenceDB))
 return(struct)
})


MSMSimporter <- function(IL, filelist) {
  ##Extract all MSMS file data (header and peaks)
  filesdata <- mapply(function(f, n) {
    hand <- mzR::openMSfile(f)
    pks <- mzR::peaks(hand)
    h <- mzR::header(hand)
    pks <- pks[which(h$msLevel == 2)]
    h <- h[which(h$msLevel == 2), ]
    h <- cbind(h, runnum = rep(n, times = nrow(h)))
    return(list(h, pks))
  }, filelist, seq_along(filelist))
  
  heads <- do.call(rbind, filesdata[seq(1, length(filesdata),
                                        2)])
  heads$totIonCurrent <- heads$totIonCurrent*heads$injectionTime/50
  
  peaks <- unlist(filesdata[seq(2, length(filesdata), 2)],
                  recursive = FALSE)
  ##Organize MSMS data into the different IL entries
  mzfilter <- IL@ILParam@filtermz
  RES <- apply(IL@IL[, c(1, 2, 3)], 1, function(entry) {
    idx <- between(heads$retentionTime, entry[1], entry[2]) &
      between(heads$precursorMZ, entry[3] - mzfilter, entry[3] +
                mzfilter)
    if (!any(idx)) {
      return(list())
    }
    subh <- heads[idx, ]
    subpks <- peaks[idx]
    subpks <- do.call(rbind, lapply(seq_along(subpks), function(x) {
      s <- subpks[[x]]
      s <- cbind(s, rep(as.numeric(subh[x, "retentionTime"]),
                        times = nrow(s)), rep(as.numeric(subh[x, "runnum"]),
                                              times = nrow(s)))
      s[,2] <- s[,2] * subh[x, "injectionTime"] / 50 #IT scaling, 50ms as reference
      return(s)
    }))
    subpks <- as.data.frame(subpks)
    names(subpks) <- c("mz", "int", "rt", "ID")
    subpks <- filter(subpks, int > 1)  #Remove zeros
    if (nrow(subpks) == 0) {
      return(list())
    }
    return(list(subh, subpks))
  })
  return(RES)
}

#### MS2 spectra processing functions ####

#' @title CliqueMSMS
#' @description Purify continuous MSMS spectra detecting consistent mass signals
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
                       BPPARAM = BiocParallel::SerialParam(), sstype = "regular") {
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
  if(sstype == "regular"){
    suppressWarnings(
      RES <- bplapply(seq_along(idx), generate_ss, BPPARAM = BPPARAM,
                      MS2list, contaminant, delta, fs, idx, FALSE)
    )
    RES <- do.call(rbind, RES)
    RES <- RHermes:::purify_ss(RES, minint = 30000, minpks = 2)
  } else {
    RES <-  lapply(seq_along(idx), wrapper_failsafe_ss, idx = idx,
                   MS2list = MS2list, fs = fs)
    RES <- do.call(rbind, RES)
  }
  
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
  
  #Trying to salvage a spectra from suboptimal data
  if (length(soi) == 0){
    if(any(data$rtiv > 30000)){
      return(failsafe_ss(data, header, idx[curentry], fs[curentry]))
    } else {return()}
  } else if (length(soi) < 3) {
    #In case we have a few traces but we're missing the most intense signal
    #registered in the data
    if(max(vapply(soi, function(x){max(x$rtiv)},
                  FUN.VALUE = numeric(1))) < max(data$rtiv)){
      return(failsafe_ss(data, header, idx[curentry], fs[curentry]))
    }
  }
  
  
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
  
  ##Calculate all similarities
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
  
  
  ##Network partitioning
  comp <- components(net)
  members <- comp$membership
  for(i in unique(comp$membership)){
    vertices <- which(comp$membership == i)
    if(length(vertices) == 1){next}
    current_net <- induced_subgraph(net, vertices)
    partitioning <- cluster_fast_greedy(current_net)
    dens <- edge_density(current_net)
    mod <- modularity(current_net, membership(partitioning))
    if(is.nan(dens)){dens <- 0}
    if(dens < 0.5 & mod > 0.3){
      members[vertices] <- (1e3*i + membership(partitioning)) #1e3 as arbitrary number to avoid membership collisions
    }
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

wrapper_failsafe_ss <- function(curentry, idx, MS2list, fs){
  header <- MS2list[[curentry]][[1]]
  data <- MS2list[[curentry]][[2]]
  names(data)[2] <- "rtiv"
  
  #Select by TIC
  besttic_rt <- header$retentionTime[which.max(header$totIonCurrent)]
  data <- data[data$rt == besttic_rt, ]
  header <- header[which.max(header$totIonCurrent), ]
  return(failsafe_ss(data, header, idx[curentry], fs[curentry]))
}

failsafe_ss <- function(data, header, idx, fs){
  maxt <- data$rt[which.max(data$rtiv)]
  ss <- data.frame(start = numeric(1), end = numeric(1),
                   apex = numeric(1), ILentry = numeric(1),
                   precmass = numeric(1))
  ss$start <- min(data$rt)
  ss$end <- max(data$rt)
  data <- data[data$rt == maxt, c("mz", "rtiv")]
  data <- data[data$rtiv > 0.005*max(data$rtiv),]
  ss$ssdata <- list(data.frame(mz=data$mz, int = data$rtiv))
  ss$apex <- maxt
  ss$precmass <- header$precursorMZ[1]
  ss$ILentry <- idx
  ss$anot <- fs
  return(ss)
}
