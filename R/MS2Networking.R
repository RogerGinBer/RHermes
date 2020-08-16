#
# MS2Networking <- function(struct, ID, costhr) {
#     validObject(struct)
#     MS2Exp <- struct@data@MS2Exp[[ID]]
#     identdf <- MS2Exp@Ident[[1]]
#     withScans <- MS2Exp@Ident[[4]]
#
#     #Tidying info
#     tags <- lapply(seq_along(MS2Exp@Ident[[2]][[2]]), function(ILid) {
#         l <- length(MS2Exp@Ident[[2]][[2]][[ILid]])
#         if (l == 0) {
#             return()
#         }
#         return(paste(rep(withScans[ILid], times = l), seq_len(l),
#             sep = "_"))
#     }) %>% unlist()
#
#     ident_tags <- lapply(seq_len(nrow(identdf)), function(row) {
#         cur <- identdf[row, ]
#         return(data.frame(id = paste(cur$IL_ID, cur$ss, sep = "_"),
#             comp = cur$compound, smiles = cur$smiles, stringsAsFactors = FALSE))
#     }) %>% do.call(rbind, .)
#
#     smiles <- c()
#     for (i in seq_len(nrow(ident_tags))) {
#         cur <- ident_tags[i, ]
#         tomatch <- cur$id
#         smiles <- c(smiles, cur$smiles)
#         tags[tags == cur$id] <- paste(cur$comp, cur$id, sep = "@")
#     }
#     identified <- which(grepl("@", tags, fixed = TRUE))
#
#     subset <- seq_along(tags)
#     #Calculating cosines
#     allspec <- unlist(MS2Exp@Ident[[2]][[2]], recursive = FALSE)[subset]
#     cos <- lapply(allspec, function(pattern) {
#         lapply(allspec, function(query) {
#             MSMScosineSim(pattern, query, minhits = 1, mode = "full")
#         })
#     }) %>% unlist() %>% matrix(nrow = length(allspec), byrow = TRUE)
#     colnames(cos) <- tags[subset]
#     rownames(cos) <- tags[subset]
#
#     #Graph and output
#     gr <- igraph::graph_from_adjacency_matrix(ifelse(cos > costhr,
#         1, 0))
#     wc <- igraph::cluster_walktrap(gr)
#     members <- igraph::membership(wc)
#     gr <- igraph::simplify(gr, remove.multiple = TRUE, remove.loops = TRUE)
#
#     gr <- networkD3::igraph_to_networkD3(gr, group = members)
#
#     #SMILES parsing and plotting to generate images for the network
#     txt <- lapply(smiles, function(smile) {
#         img <- view.image.2d(parse.smiles(smile)[[1]])
#         png("test.png")  #Temp file
#         grid::grid.raster(img)
#         dev.off()
#         return(RCurl::base64Encode(readBin("./test.png", "raw",
#             file.info("./test.png")[1, "size"]), "txt"))
#     })
#     file.remove("test.png")
#
#
#     gr$nodes$shape <- "circle"
#     gr$nodes$image <- NA
#     for (i in seq_along(identified)) {
#         id <- identified[i]
#         gr$nodes$shape[id] <- "image"
#         gr$nodes$image[id] <- paste("data:image/png;base64",
#             txt, sep = ",")
#     }
#     net <- visNetwork(gr$nodes %>% rename(label = name) %>% mutate(id = seq_len(nrow(gr$nodes)) -
#         1), gr$links %>% rename(from = source, to = target),
#         width = "1200", height = "1200")
#
#     net %<>% visNodes(color = list(background = "lightblue")) %>%
#         visEdges(smooth = FALSE) %>% visPhysics(solver = "barnesHut",
#         stabilization = TRUE)
#
#     net
#     # networkD3::forceNetwork(Links = gr$links, Nodes = gr$nodes,
#     #                         Source = 'source', Target = 'target',
#     #                         NodeID = 'name', Group = 'group',
#     #                         height = 3000, width = 3000, zoom = TRUE,
#     #                         opacity = 0.8)
# }
