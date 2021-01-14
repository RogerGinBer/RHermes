library(data.table)
library(magrittr)
library(visNetwork)
library(igraph)
library(BiocParallel)

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
      idesp <- espmet[.(id), ]$ID_spectra[metadata$REFCE %in% c("0","10","20","10eV")]
      metadata <- metadata[metadata$REFCE %in% c("0","10","20","10eV"), ]
      
      #Filtering by adduct
      idesp <- espmet[.(id), ]$ID_spectra[metadata$REFadduct %in% c("M+H","M-H", "M+")]
      metadata <- metadata[metadata$REFadduct %in% c("M+H","M-H", "M+"), ]
      which_polarity <- which(metadata$REFpolarity == polarity)
      spec <- list_fragments[.(idesp[which_polarity]), "spectra"] %>% unlist(., recursive = FALSE)
      if(length(spec)==0){return(spec)}
      energies <- apply(metadata[which_polarity, c("REFprecursor_mz","REFCE", "REFnature")], 1,
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
    common_deltas <- lapply(RES, function(structure){
      prec_m <- sapply(names(structure), function(x){as.numeric(strsplit(x, "_")[[1]][1])})
      structure <- structure[!is.na(prec_m)]
      prec_m <- prec_m[!is.na(prec_m)]
      observed_mzs <- c()
      for(i in seq_along(structure)){
        spectrum <- as.matrix(t(structure[[i]]))
        maxi <- max(spectrum[,2])
        spectrum <- spectrum[spectrum[,2] > 0.3*maxi, , drop = FALSE]
        deltas <- spectrum[,1] - prec_m[i]
        observed_mzs <- c(observed_mzs, spectrum[abs(deltas) > 0.01 & deltas < 0, 1])
      }
      if(length(observed_mzs) == 0){return()}
      return(sort(unique(observed_mzs)))
    })
    common_deltas <- common_deltas[sapply(common_deltas, length) != 0]
    unlist(common_deltas) %>% unique() %>% sort()
    
  })
}

anotateISF <- function(SL, DB, polarity = 1, BiocParallelParam = SerialParam()){
  DB$df_spectra <- as.data.table(DB$df_spectra)
  DB$df_metabolite <- as.data.table(DB$df_metabolite)
  DB$df_spectraMetabolite <- as.data.table(DB$df_spectraMetabolite)
  DB$list_fragments <- data.table(ID_spectra = DB$list_fragments[[1]],
                                  spectra = DB$list_fragments[[2]])
  setkeyv(DB$df_spectra, c("ID_spectra"))
  setkeyv(DB$df_spectraMetabolite, c("ID_metabolite"))
  setkeyv(DB$df_metabolite, c("REFformula"))
  setkeyv(DB$list_fragments, c("ID_spectra"))
  
  SL$ISF <- bplapply(seq_len(nrow(SL)), RHermes:::anotateParallelISF,
                     SL = SL, DB = DB, polarity = polarity, 
                     BiocParallelParam = BiocParallelParam)
  return(SL)
}

anotateParallelISF <- function(i, SL, DB, polarity){
  cur <- SL[i, ]
  masses <- generateDiffDB(DB, cur$f[[1]], polarity)
  masses <- masses[sapply(masses, length) != 0]
  if(length(masses) == 0){return(integer())}
  ISF <- lapply(masses, function(x){
    diffs <- do.call(cbind, lapply(x, function(mass){
      SL$mass - mass}))
    candidates <- which(apply(diffs, 1, function(m){any(abs(m) < 0.02)}))
    cosines <- sapply(candidates, function(cand){
      RHermes:::cosineSim(cur$peaks[[1]], SL$peaks[[cand]], nscans = 5)
    }) 
    candidates[cosines > 0.95 & candidates != i]
  })
  return(unlist(ISF))
}

plotISF <- function(SOIlist){
  net <- graph_from_adj_list(SOIlist$ISF)
  browser()
  net <- set.vertex.attribute(net, "value", value = log10(SOIlist$MaxInt),
                              index = which(seq_along(V(net)) %in% unique(SOIlist$originalID))) 
  net <- set.vertex.attribute(net, "title", value = paste0("<p>Mz: ",SOIlist$mass,"</p>",
                                                           "<p>Intensity: ", round(SOIlist$MaxInt,digits = 2), "</p>",
                                                           "<p>Anot: ", sapply(SOIlist$anot, function(an){paste(an, collapse = ",")}), "</p>"),
                              index = which(seq_along(V(net)) %in% unique(SOIlist$originalID)))
  col <- colorRamp(c(rgb(1,0,0),rgb(0,1,0), rgb(0,0,1)), bias = 0.1)
  defaultcol <- rep("#777777", length(V(net)))  
  node_color <- apply(col(SOIlist$mass/max(SOIlist$mass)), 1, function(x){
    x <- x/255
    rgb(x[1],x[2],x[3])})
  defaultcol[which(seq_along(V(net)) %in% unique(SOIlist$originalID))] <- node_color
  
  net <- set.vertex.attribute(net, "color", value = defaultcol)
  visnet <- toVisNetworkData(net)
  
  # distances <- apply(visnet[[2]], 1, function(x){
  #   SOIlist$mass[[x[2]]] - SOIlist$mass[[x[1]]]
  # })
  # edge_color <- apply(col(abs(distances)/max(abs(distances))), 1, function(x){
  #   x <- x/255
  #   rgb(x[1],x[2],x[3])})
  # 
  # visnet[[2]]$label <-  as.character(round(distances, 5))
  
  visNetwork(nodes = visnet[[1]], edges = visnet[[2]]) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5))) %>% 
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>% 
    visPhysics(stabilization = FALSE)
  
}

removeISF <- function(SL){
  net <- igraph::graph_from_adj_list(SL$ISF)
  cl <- igraph::groups(igraph::cluster_walktrap(net))
  SL$originalID <- seq(nrow(SL))
  do.call(rbind, lapply(cl, function(group){
    subnet <- subgraph(net, group)
    df <- data.frame(out = degree(subnet, mode = "out"),
                     into = degree(subnet, mode = "in"))
    isf <- which(df$int > df$out)
    if(length(isf) == 0){
      sois <- SL[group,]
      sois$group <- rep(list(group), nrow(sois))
      return(sois)
    }
    mInt <- max(SL$MaxInt[group[-isf]])
    isf <- isf[SL$MaxInt[group[isf]] < mInt]
    if(length(isf) != 0){group <- group[-isf]}
    sois <- SL[group,]
    sois$group <- rep(list(group), nrow(sois))
    return(sois)
  }))
}

#' @title ISFproc
#' @description Detect and remove ISF signals from a SOI list using low
#' collision energy MS2 data
#' @export
ISFproc <- function(struct, id, DBpath = "D:/MS2ID_B2R_20201113_083214.rds"){
  DB <- readRDS(DBpath)
  BiocParallelParam <- struct@metadata@cluster
  polarity <- ifelse(struct@metadata@ExpParam@ion == "+", 1, 0)
  #Anotate and remove ISF
  SoiObj <- struct@data@SOI[[id]]
  SL <- SoiObj@SoiList
  message("Anotating ISF:")
  SL <- anotateISF(SL, DB, polarity, BiocParallelParam)
  message("Creating ISF network and cleaning:")
  SL <- removeISF(SL)
  SL <- as.data.table(SL)
  setkeyv(SL, "formula")
  
  #Recalculate plotting DF
  plist <- bplapply(unique(SL$formula), RHermes:::preparePlottingDF,
                    SL, BPPARAM = BiocParallelParam)
  plist <- do.call(rbind, plist)
  plist$isov <- rep("M0", nrow(plist))
  
  #Update and return
  SoiObj@SoiList <- SL
  SoiObj@PlotDF <- as.data.table(plist)
  struct@data@SOI[[id]] <- SoiObj
  struct <- RHermes:::setTime(struct, paste("Removed ISF from SOI list", id,
                                  "using database file", "DBpath"))
  return(struct)
}


calculateOverlap <- function(file){
  myHermes <- readRDS(file)
  m <- sort(myHermes@metadata@ExpParam@ionF[[1]]$m)
  diffs <- lapply(seq_along(m), function(x){
    deltas <- (m[seq(max(x-200,0), min(x+200, length(m)))] - m[x])/ m[x] * 1e6
    return(abs(deltas[abs(deltas) < 10]))
  })
  collisiondata <- lapply(seq(0,5), function(colisions){
    message(colisions)
    lapply(seq(0,10,0.1), function(ppm){
      length(which(sapply(diffs, function(deltas){
        length(which(deltas < ppm & deltas != 0)) <= colisions
      })))
    })
  })
  df_to_plot <- lapply(seq_along(collisiondata), function(x){
    cur <- unlist(collisiondata[[x]])
    df <- data.frame(x=cur, ppm = seq(0,10,0.1), overlaps = rep(x-1, length(cur)))
  })
  df_to_plot <- do.call(rbind, df_to_plot)
  print(ggplot(df_to_plot) + geom_point(aes(x=ppm, y = x/length(m)*100, color = overlaps)) +
          scale_y_continuous(limits = c(0,100)) + scale_x_continuous(limits = c(0,10))+
          theme_minimal()+
          ylab("% of uniquely distinguishable ionic formulas"))
  return(unlist(diffs))  
}
