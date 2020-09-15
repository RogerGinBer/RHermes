
#  PLPlotBlank <- function(Sample, Blank, filename = 'PL_Sample+Blank'){
#   if(class(Sample) != 'rHERMES_PL' | class(Blank) != 'rHERMES_PL'){
#     stop('Please input valid Sample and Blank objects')}
#   cairo_pdf(filename = filename, width = 15, onefile = TRUE)
#   setkey(Sample@peaklist, formv)
#   setkey(Blank@peaklist, formv)
#   pb <- txtProgressBar(min = 0, max = length(unlist(unique(Sample@peaklist[,formv]))),style = 3)
#   count <- 1
#   for(i in unlist(unique(Sample@peaklist[,formv]))){
#     setTxtProgressBar(pb, count)
#     count <- count + 1
#     curSample <- Sample@peaklist[.(i)]
#     curBlank <- Blank@peaklist[.(i)]
#     curSample$ID <- 'Sample'
#     curBlank$ID <- 'Blank'
#     data <- rbind(curSample, curBlank)
#     goodad <- vapply(unique(data$adduct), function(x){
#       cur <- data[data$adduct == x, ]
#       if(nrow(cur)>30){return(TRUE)}
#       else{return(FALSE)}
#     })
#     goodad <- unique(data$adduct)[goodad]
#     data <- data[data$adduct %in% goodad, ]
#     if(nrow(data) < 30){next}
#     print(ggplot(data) + geom_point(aes(x = rt, y = rtiv, color = isov),
#                                     size = 0.05, alpha = 0.8) +
#             labs(x = 'RT', y = 'Log Intensity', title = i) + scale_y_log10()+
#             facet_grid(rows = vars(adduct), cols = vars(ID)) + theme_minimal())
#   }
#   dev.off()
# }
#
# PLLabelled <- function(Labelled, Unlabelled, filename = 'PL_Labelled+Unlabelled.pdf'){
#   if(class(Labelled) != 'rHERMES_PL' | class(Unlabelled) != 'rHERMES_PL'){
#     stop('Please input valid Labelled and Unlabelled objects')}
#   cairo_pdf(filename = filename, onefile = TRUE)
#   setkey(Labelled@peaklist, formv)
#   setkey(Unlabelled@peaklist, formv)
#   pb <- txtProgressBar(min = 0, max = length(unlist(unique(Labelled@peaklist[,formv]))),style = 3)
#   count <- 1
#   for(i in unlist(unique(Labelled@peaklist[,formv]))){
#     setTxtProgressBar(pb, count)
#     count <- count + 1
#     curLabelled <- Labelled@peaklist[.(i)]
#     curUnlabelled <- Unlabelled@peaklist[.(i)]
#     curLabelled$ID <- 'Labelled'
#     curUnlabelled$ID <- 'Unlabelled'
#     data <- rbind(curLabelled, curUnlabelled)
#     goodad <- vapply(unique(data$adduct), function(x){
#       cur <- data[data$adduct == x, ]
#       if(nrow(cur)>30){return(TRUE)}
#       else{return(FALSE)}
#     })
#     goodad <- unique(data$adduct)[goodad]
#     data <- data[data$adduct %in% goodad, ]
#     if(nrow(data) < 30){next}
#     for(ad in goodad){
#       addata <- data[data$adduct == ad, ]
#       goodiso <-vapply(unique(addata$isov), function(x){
#         cur <- addata[addata$isov == x, ]
#         if(nrow(cur)>10){return(TRUE)}
#         else{return(FALSE)}
#       })
#       goodiso <- unique(as.character(addata$isov))[goodiso]
#       addata <- filter(addata, isov %in% goodiso)
#       print(ggplot(addata) + geom_point(aes(x = rt, y = rtiv, color = ID),
#                                         size = 0.05, alpha = 0.8) +
#               scale_color_manual(values = c('#2d3561','#c05c7e','#f3826f'))+
#               scale_y_log10()+
#               facet_grid(rows = vars(isov)) + theme_minimal()+
#               labs(x = '\nRetention time (s)', y = 'Log 10 (Intensity)\n',  title = paste(i, ad))+
#               theme_minimal()+
#               theme(plot.margin = unit(c(1,0.7,1,0.8), 'cm'), text = element_text(size = 12, family = 'Roboto Light')))
#     }
#   }
#   dev.off()
# }
#
# PLPlot <- function(Sample, filename = 'PL_Sample'){
#   if(class(Sample) != 'rHERMES_PL'){
#     stop('Please input a valid Sample object')}
#   cairo_pdf(filename = filename, width = 10, onefile = TRUE)
#   for(i in seq_len(nrow(Sample@params@DB))){
#     cat(i, ' out of: ', nrow(Sample@params@DB),'\n')
#     f <- as.character(Sample@params@DB[i, 2])
#     data <- dplyr::filter(Sample@peaklist, formv == f)
#     goodad <- vapply(unique(data$adduct), function(x){
#       cur <- data[data$adduct == x, ]
#       if(nrow(cur)>30){return(TRUE)}
#       else{return(FALSE)}
#     })
#     goodad <- unique(data$adduct)[goodad]
#     data <- data[data$adduct %in% goodad, ]
#     if(nrow(data) < 30){next}
#     print(ggplot(data) + geom_point(aes(x = rt, y = rtiv, color = isov),
#                                     size = 0.05, alpha = 0.8) +
#             labs(x = 'RT', y = 'Log Intensity', title = f) + scale_y_log10()+
#             facet_grid(rows = vars(adduct)) + theme_minimal())
#   }
#   dev.off()
# }
# SOIPlotBlank <- function(SampleSOI, SamplePL, BlankPL, filename = 'SOI_Blank+Sample'){
#   cairo_pdf(filename = filename, onefile = TRUE)
#
#   plist <- SampleSOI@PlotDF
#   datafile <- SamplePL@peaklist
#   blankfile <- BlankPL@peaklist
#   setkey(plist, form)
#   setkey(datafile, formv)
#   setkey(blankfile, formv)
#
#   ##Main loop-----------------------------------------------------------
#
#   iter <- unique(plist$form)
#   pb <- txtProgressBar(min = 0, max = length(iter), style = 3)
#   count <- 1
#   for (i in iter){
#     setTxtProgressBar(pb, count)
#     count <- count + 1
#
#     #Values from groups----------------------------------------------
#     pks <- plist[.(i)]
#     if(dim(pks)[1]==0){next}
#     pks$categ <- 'GR'
#     pks$sample <- 'Sample'
#
#     #Values from raw data--------------------------------------------
#     praw <- datafile[.(i)] %>% filter(., isov == 'M0')
#     praw$snr <- 0
#     praw$categ <-'RAW'
#     praw$sample <- 'Sample'
#     praw <- praw[ ,c('rt','rtiv','formv','adduct','snr','isov','categ','sample')]
#     names(praw) <- names(pks)
#
#
#     #Values from blank raw data--------------------------------------
#     bpraw <- blankfile[.(i)] %>% filter(., isov == 'M0')
#     if(nrow(bpraw) != 0){
#       bpraw$snr <- 0
#       bpraw$categ <-'RAW'
#       bpraw$sample <- 'Blank'
#       bpraw <- bpraw[ ,c('rt','rtiv','formv','adduct','snr','isov','categ','sample')]
#       names(bpraw) <- names(pks)
#       pfull <- rbind(pks, praw, bpraw)
#     }else{
#       pfull <- rbind(pks, praw) #Case when there aren't any scans on the blank
#     }
#
#     for(j in unique(pfull$ad)){ #Clear adducts with few scans
#       if(length(which(pfull$ad == j)) < 5){
#         pfull <- pfull[pfull$ad != j, ]
#       }
#     }
#
#
#     if(dim(pfull)[1]>5){ #Check for minimum number of scans
#       print(ggplot()+
#               geom_point(data = pfull[pfull$categ == 'GR' & pfull$sample == 'Sample',],
#                          mapping = aes(x = rt/60, y = log10(rtiv), colour = snr),
#                          alpha=0.4, size=5,shape=20)+
#               geom_point(data = pfull[pfull$categ == 'RAW' & pfull$sample == 'Sample',],
#                          mapping = aes(x = rt/60, y = log10(rtiv)), colour = '#000000',
#                          alpha=0.5,size=0.8,shape=15)+
#               geom_point(data = pfull[pfull$sample == 'Blank',],
#                          mapping = aes(x = rt/60, y = log10(rtiv)), colour = '#CC0000',
#                          alpha=0.5,size=0.8,shape=15)+
#               facet_grid(rows = vars(ad))+
#               xlab('RT (min)')+ylab('Log 10 Intensity\n')+
#               scale_color_continuous(type = 'viridis')+
#               ggtitle(label = i)+
#               theme_minimal()+
#               theme(plot.margin = unit(c(1,0.4,1,1), 'cm'),
#                     text = element_text(size = 12, family = 'Roboto Light')))
#     }
#     else {print(i)}
#   }
#   dev.off()
# }
#
# SOIPlot <- function(SampleSOI, SamplePL, filename = 'SOI_Sample'){
#   cairo_pdf(filename = filename, onefile = TRUE)
#
#   plist <- SampleSOI@PlotDF
#   datafile <- SamplePL@peaklist
#
#   ##Main loop-----------------------------------------------------------
#
#   iter <- unique(plist$form)
#   pb <- txtProgressBar(min = 0, max = length(iter), style = 3)
#   count <- 1
#   for (i in iter){
#     setTxtProgressBar(pb, count)
#     count <- count + 1
#
#     #Values from groups----------------------------------------------
#     pks <- plist[plist$form == i, ]
#     if(dim(pks)[1]==0){next}
#     pks$categ <- 'GR'
#
#     #Values from raw data--------------------------------------------
#     praw <- datafile[formv == i  & isov == 'M0',]
#     praw$snr <- 0
#     praw$categ <-'RAW'
#     praw <- praw[ ,c('rt','rtiv','formv','adduct','snr','isov','categ')]
#     names(praw) <- names(pks)
#
#     pfull <- rbind(pks, praw)
#
#     for(j in unique(pfull$ad)){ #Clear adducts with few scans
#       if(length(which(pfull$ad == j)) < 5){
#         pfull <- pfull[pfull$ad != j, ]
#       }
#     }
#     pfull$sample <- 'Sample'
#
#     if(dim(pfull)[1]>5){ #Check for minimum number of scans
#       print(ggplot()+
#               geom_point(data = pfull[pfull$categ == 'GR' & pfull$sample == 'Sample',],
#                          mapping = aes(x = rt/60, y = log10(rtiv), colour = snr),
#                          alpha=0.4, size=5,shape=20)+
#               geom_point(data = pfull[pfull$categ == 'RAW' & pfull$sample == 'Sample',],
#                          mapping = aes(x = rt/60, y = log10(rtiv)), colour = '#000000',
#                          alpha=0.5,size=0.8,shape=15)+
#               facet_grid(rows = vars(ad))+
#               xlab('RT (min)')+ylab('Log 10 Intensity\n')+
#               scale_color_continuous(type = 'viridis')+
#               ggtitle(label = i)+
#               theme_minimal()+
#               theme(plot.margin = unit(c(1,0.4,1,1), 'cm'),
#                     text = element_text(size = 12, family = 'Roboto Light')))
#     }
#     else {print(i)}
#   }
#   dev.off()
# }
#
#
# ThermoPlot <- function(markedPL, SamplePL, SampleSOI, filename){
#   cairo_pdf(filename = filename, onefile = TRUE)
#   plist <- SampleSOI@PlotDF
#   datafile <- SamplePL@peaklist
#   markedfile <- markedPL@peaklist
#
#   ##Main loop-----------------------------------------------------------
#   setkey(plist, form)
#   setkey(datafile, formv)
#   setkey(markedfile, formv)
#
#   formulas <- unique(plist$form)
#   adducts <- unique(plist$ad)
#   pb <- txtProgressBar(min = 0, max = length(formulas), style = 3)
#   count <- 0
#   for (i in formulas){
#     setTxtProgressBar(pb, count)
#     count <- count + 1
#     for (j in adducts){
#
#
#       #Values from groups----------------------------------------------
#       pks <- plist[.(i)]
#       pks <- pks[pks$ad == j & pks$isov == 'M0']
#       if(dim(pks)[1]==0){next}
#       pks$categ <- 'GR'
#
#       #Values from raw data--------------------------------------------
#       praw <- datafile[.(i)]
#       praw <- praw %>% filter(., adduct == j & isov == 'M0')
#       praw$snr <- 0
#       praw$categ <-'RAW'
#       praw <- praw[ ,c('rt','rtiv','formv','adduct','snr','isov','categ')]
#       names(praw) <- names(pks)
#
#       #Values from marked data-----------------------------------------
#       pmarked <- markedfile[.(i)]
#       pmarked <- pmarked %>% filter(., adduct == j & isov != 'M0')
#       if(dim(pmarked)[1]==0){next}
#       pmarked$snr <- 0
#       pmarked$categ <- 'MAR'
#       pmarked <- pmarked[ ,c('rt','rtiv','formv','adduct','snr','isov','categ')]
#       for(j in unique(pmarked$isov)){ #Clear isotopes with few scans
#         if(length(which(pmarked$isov == j)) < 5){
#           pmarked <- pmarked[pmarked$isov != j, ]
#         }
#       }
#
#       pfull <- rbind(pks, praw, pmarked, use.names = FALSE)
#
#       for(j in unique(pfull$ad)){ #Clear adducts with few scans
#         if(length(which(pfull$ad == j)) < 5){
#           pfull <- pfull[pfull$ad != j, ]
#         }
#       }
#       # pfull$sample <- 'Sample'
#
#       if(dim(pfull)[1]>5){ #Check for minimum number of scans
#         print(ggplot()+
#                 geom_point(data = pfull[pfull$categ == 'GR',],
#                            mapping = aes(x = rt/60, y = log10(rtiv), colour = snr),
#                            alpha=0.4, size=5,shape=20)+
#                 geom_point(data = pfull[pfull$categ == 'RAW',],
#                            mapping = aes(x = rt/60, y = log10(rtiv)), colour = '#000000',
#                            alpha=0.5,size=0.8,shape=15)+
#                 geom_point(data = pfull[pfull$categ == 'MAR',],
#                            mapping = aes(x = rt/60, y = log10(rtiv)), colour = '#f05c15',
#                            alpha=0.5,size=0.8,shape=15)+
#                 facet_grid(rows = vars(isov))+
#                 xlab('RT (min)')+ylab('Log 10 Intensity\n')+
#                 scale_color_continuous(type = 'viridis')+
#                 ggtitle(label = paste(i, j, sep = '-'))+
#                 theme_minimal()+
#                 theme(plot.margin = unit(c(1,0.4,1,1), 'cm'),
#                       text = element_text(size = 12, family = 'Roboto Light')))
#       }
#       else {print(i)}
#     }}
#   dev.off()
# }

#'@export
setGeneric("MirrorPlot", function(struct, ms2id, ssnumber, patform) {
    standardGeneric("MirrorPlot")
})
setMethod("MirrorPlot", c("RHermesExp", "numeric", "numeric"),
    function(struct, ms2id, ssnumber, patform) {
        # function(query, pattern, title = 'Mirror plot', subtitle = '', baseline = 1000, maxint = 0, molecmass = 200){
        entry <- struct@data@MS2Exp[[ms2id]]@Ident$MS2Features[ssnumber,]
        query <- entry$ssdata[[1]]
        ref <- struct@data@MS2Exp[[ms2id]]@Ident$MS2_correspondance[[ssnumber]]
        pattern <- struct@data@MS2Exp[[ms2id]]@Ident$DatabaseSpectra[ref]
        pattern <- unname(pattern) #Avoid "name.subname" when unlisting
        pattern <- unlist(pattern, recursive = FALSE, use.names = TRUE)
        maxint <- max(query$int)
        query$int <- query$int/max(query$int) * 100
        bestdf <- query[query$int > 10, ]
        bestdf$mz <- round(bestdf$mz, 4)

        molecmass <- entry$precmass
        baseline <- 1000
        subtitle <- ""
        title <- ""
        baseline <- baseline/maxint * 100

        mirrplot <- lapply(patform, function(x){
            refspec <- pattern[[x]][[entry$results[[1]]$id[which(entry$results[[1]]$formula == x)]]]
            if(is.null(refspec)){return(ggplotly(ggplot()))}
            refspec <- as.data.frame(t(refspec))

            colnames(refspec) <- c("mz", "int")
            refspec$int <- refspec$int/max(refspec$int) * 100

             bldf <- data.frame(xmin = min(c(refspec$mz, query$mz,
                                molecmass)) - 5, xmax = max(c(refspec$mz,
                                                              query$mz,
                                molecmass) + 5), y = baseline)

             moldf <- data.frame(mz = molecmass)


             pl <- ggplot() + geom_segment(data = query, aes(x = mz,
                xend = mz, y = 0, yend = int), color = "black") +
                geom_segment(data = refspec, aes(x = mz, xend = mz,
                             y = 0, yend = -int), color = "red") +
                geom_segment(data = bldf, aes(x = xmin, xend = xmax, y = y,
                                              yend = y), linetype = "dashed",
                             color = "black", alpha = 0.3) +
                geom_segment(data = bldf, aes(x = xmin, xend = xmax, y = -y,
                                              yend = -y), linetype = "dashed",
                             color = "red", alpha = 0.3) +
                geom_point(data = moldf, aes(x = mz, y = 0), shape = 17,
                           size = 2) + theme_minimal() + ylab("% Intensity") +
                scale_x_continuous(limits = c(min(c(refspec$mz,
                                                    query$mz, molecmass)) - 20,
                                              max(c(refspec$mz, query$mz,
                                                    molecmass)) + 20)) +
                theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
                      text = element_text(size = 11, family = "Segoe UI Light"),
                      plot.title = element_text(hjust = 0.5)) +
                geom_text(data = bestdf, aes(x = mz, y = int + 5, label = mz),
                          family = "Segoe UI Light", check_overlap = TRUE) +
                ggtitle(title, subtitle)

              ggplotly(pl, height = ifelse(length(patform)<5, 300, 200)*
                           length(patform))
            })
        return(subplot(mirrplot, nrows = length(mirrplot), shareX = TRUE,
                       which_layout = 1))
    })

#'@export
setGeneric("PLPlot", function(struct, id, formula, rtrange,
    dynamicaxis, ads) {
    standardGeneric("PLPlot")
})
setMethod("PLPlot", c("RHermesExp", "numeric", "character",
    "numeric", "logical", "character"), function(struct, id,
    formula, rtrange, dynamicaxis, ads) {
    # plist <- struct@data@SOI[[id]]@PlotDF
    # filen <- struct@data@SOI[[id]]@filename
    # plid <- which(vapply(struct@data@PL, function(x){return(x@filename == filen)}, logical(1)))
    datafile <- struct@data@PL[[id]]@peaklist
    FA_to_ion <- struct@metadata@ExpParam@ionF[[2]]
    setkey(FA_to_ion, "f")

    fs <- FA_to_ion[f == formula, ]
    ions <- fs$ion
    datafile <- filter(datafile, formv %in% ions)
    if (nrow(datafile) == 0) {return()}
    datafile <- filter(datafile, between(rt, rtrange[1], rtrange[2]))
    if (nrow(datafile) == 0) {return()}

    # toKeep <- vapply(unique(datafile$formv), function(f){
    #   ifelse(length(which(datafile$formv == f)) > 10, T, F)
    # }, logical(1))
    # datafile <- datafile[datafile$formv %in% unique(datafile$formv)[toKeep], ]
    # if(nrow(datafile) == 0){return()}

    datafile$ad <- ""
    for (f in unique(datafile$formv)) {
        ad <- fs$an[fs$ion == f]
        datafile$ad[datafile$formv == f] <- ad
    }
    datafile <- filter(datafile, ad %in% ads)
    if (nrow(datafile) == 0) {return()}

    return(ggplotly(ggplot() + geom_point(data = datafile, mapping = aes(x = rt,
        y = rtiv, color = isov), size = 0.5, alpha = 0.6) + facet_grid(rows = vars(ad)) +
        theme_minimal() + ggtitle(formula) + scale_color_brewer(palette = "Dark2"),
        dynamicTicks = dynamicaxis))

})


#'@export
setGeneric("SoiPlot", function(struct, id, formula, rtrange,
    dynamicaxis, ads, blankid = NA) {
    standardGeneric("SoiPlot")
})
setMethod("SoiPlot", c("RHermesExp", "numeric", "character",
    "numeric", "logical", "character", "ANY"), function(struct, id,
    formula, rtrange, dynamicaxis, ads, blankid = NA) {
    plist <- struct@data@SOI[[id]]@PlotDF
    filen <- struct@data@SOI[[id]]@filename
    plid <- which(vapply(struct@data@PL, function(x) {
        return(x@filename == filen)
    }, logical(1)))[1]



    datafile <- struct@data@PL[[plid]]@peaklist
    datafile <- datafile[datafile$isov == "M0", ]
    datafile$Class <- "Sample"

    if(!is.na(blankid)){
        blankfile <- struct@data@PL[[blankid]]@peaklist
        blankfile <- blankfile[blankfile$isov == "M0", ]
        blankfile$Class <- "Blank"
        datafile <- rbind(datafile, blankfile)
    }

    FA_to_ion <- struct@metadata@ExpParam@ionF[[2]]
    setkey(FA_to_ion, "f")
    fs <- FA_to_ion[f == formula, ]
    ions <- fs$ion
    datafile <- filter(datafile, formv %in% ions)
    datafile <- filter(datafile, between(rt, rtrange[1], rtrange[2]))
    soiinfo <- filter(plist, form %in% ions)
    soiinfo <- filter(soiinfo, between(rt, rtrange[1], rtrange[2]))
    if (nrow(soiinfo) == 0) {return()}
    names(soiinfo)[names(soiinfo) == "form"] <- "formv"
    names(soiinfo)[names(soiinfo) == "rtiv"] <- "Intensity"
    names(datafile)[names(datafile) == "rtiv"] <- "Intensity"


    # toKeep <- vapply(unique(datafile$formv), function(f){
    #   ifelse(length(which(datafile$formv == f)) > 30, T, F)
    # }, logical(1))
    # datafile <- datafile[datafile$formv %in% unique(datafile$formv)[toKeep], ]

    if (nrow(soiinfo) == 0) {return()}

    datafile$ad <- ""
    for (f in unique(datafile$formv)) {
        ad <- fs$an[fs$ion == f]
        datafile$ad[datafile$formv == f] <- ad
    }

    soiinfo$ad <- ""
    for (f in unique(soiinfo$formv)) {
        ad <- fs$an[fs$ion == f]
        soiinfo$ad[soiinfo$formv == f] <- ad
    }

    datafile <- filter(datafile, ad %in% ads)
    soiinfo <- filter(soiinfo, ad %in% ads)

    if (nrow(soiinfo) == 0) {return()}
    soiinfo$Class <- "Sample-SOI"
    if(!is.na(blankid)){
        return(ggplotly(ggplot() +
                            geom_point(data = datafile[datafile$Class == "Sample", ],
                                       mapping = aes(x = rt, y = Intensity, color = Class),
                                       alpha = 0.4) +
                            geom_point(data = datafile[datafile$Class == "Blank", ],
                                       mapping = aes(x = rt, y = Intensity, color = Class),
                                       alpha = 0.5)+
                            geom_point(data = soiinfo, mapping = aes(x = rt, y = Intensity, color = Class),
                                       alpha = 0.8) +
                            scale_color_manual(breaks = c("Blank", "Sample", "Sample-SOI"),
                                               values = c("#6D9503", "#8E032B", "#370B6B"))+
                            facet_grid(rows = vars(ad)) + theme_minimal() +
                            ggtitle(formula), dynamicTicks = dynamicaxis))
    } else {
        return(ggplotly(ggplot() +
                            geom_point(data = datafile[datafile$Class == "Sample", ],
                                       mapping = aes(x = rt, y = Intensity, color = Class),
                                        alpha = 0.3)+
                            geom_point(data = soiinfo, mapping = aes(x = rt, y = Intensity, color = Class),
                                        alpha = 0.8) +
                            scale_color_manual(breaks = c("Sample", "Sample-SOI"),
                                               values = c("#8E032B", "#370B6B"))+
                            facet_grid(rows = vars(ad)) + theme_minimal() +
                            ggtitle(formula), dynamicTicks = dynamicaxis))
    }


})

#'@export
setGeneric("RawMS2Plot", function(struct, ms2id, entryid, bymz) {
    standardGeneric("RawMS2Plot")
})
setMethod("RawMS2Plot", c("RHermesExp", "ANY", "ANY", "ANY"),
          function(struct, ms2id, entryid, bymz = TRUE) {
    if (is.na(entryid)) { 
      return(list(p_bymz = ggplotly(ggplot()), p_bygroup = ggplotly(ggplot()),
                  net = NA, pks = data.frame()))
    }
    ms2data <- struct@data@MS2Exp[[ms2id]]@MS2Data
    withdata <- which(vapply(ms2data, function(x){length(x) != 0}, logical(1)))
    if(!entryid %in% withdata){
      return(list(p_bymz = ggplotly(ggplot()), p_bygroup = ggplotly(ggplot()),
                  net = NA, pks = data.frame()))
    }
    ss <- RHermes:::generate_ss(entryid, MS2list = ms2data, contaminant = 173.5,
                      delta = 0.1, fs = chracter(), idx = numeric(),
                      to_plot = TRUE)
    soi <- ss$soi
    members <- ss$members
    net <- ss$net
    data <- ss$data
    pks <- ss$pks
    pks$members <- members
    ss <- ss$ss
    mem_to_keep <- vapply(members, function(x){
      curpks <- pks[pks$members == x, ]
      is_intense <- any(curpks$maxo > 3e4)
      has_many_pks <- nrow(curpks) > 2
      return(is_intense | has_many_pks)
    }, logical(1))
    
    members <- members[mem_to_keep]
    
    pks <- pks[pks$members %in% members, ]
    
    oldpoints <- soi
    oldpoints$member <- "Original points"
    xlim <- c(min(soi$rt)-8, max(soi$rt)+8)
    
    res <- reassign_and_check(pks, soi)
    pks <- res[[1]]
    soi <- res[[2]]
    if(nrow(soi) == 0){
      return(list(p_bymz = ggplotly(ggplot()), p_bygroup = ggplotly(ggplot()),
                  net = NA, pks = data.frame()))
    }
    
    if(any(!mem_to_keep)) net <- igraph::delete.vertices(net, which(!mem_to_keep))

    soi$member <- "Not considered"
    for (i in unique(members)) {
        soi$member[soi$peak %in% which(members == i)] <- paste("Superspec.", as.character(i))
    }

    p_bymz <- ggplotly(ggplot(rbind(soi, oldpoints)) + geom_point(aes(x = rt, y = rtiv,
                color = as.factor(mz))) + xlim(xlim))
    p_bygroup <- ggplotly(ggplot(rbind(soi, oldpoints)) + geom_point(aes(x = rt, y = rtiv,
                color = as.factor(member))) + xlim(xlim))

    net <- networkD3::igraph_to_networkD3(net, group = members)

    net <- visNetwork(net$nodes %>% rename(label = name) %>%
        mutate(id = seq_len(nrow(net$nodes)) - 1), net$links %>%
        rename(from = source, to = target))

    net %<>% visNodes(color = list(background = "lightblue")) %>%
        visEdges(smooth = FALSE) %>% visPhysics(solver = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -100),
        stabilization = FALSE)
    return(list(p_bymz = p_bymz, p_bygroup = p_bygroup, net = net, pks = pks))
})


#'@export
setGeneric("IsoFidelity", function(struct, soilist, entry,
    plot = TRUE) {
    standardGeneric("IsoFidelity")
})
setMethod("IsoFidelity", c("RHermesExp", "numeric", "numeric",
    "ANY"), function(struct, soilist, entry, plot = TRUE) {
    #Extract SOI and PL information from the selected SOI list
    SOI <- struct@data@SOI[[soilist]]
    fname <- SOI@filename
    correspondingPL <- which(struct@metadata@filenames == fname)
    PL <- struct@data@PL[[correspondingPL]]@peaklist

    #Filter the PL to the selected SOI region
    curSOI <- SOI@SoiList[entry, ]
    PL <- filter(PL, formv == curSOI$formula)
    PL <- filter(PL, data.table::between(rt, curSOI$start, curSOI$end))

    if (plot) {
        p <- ggplotly(ggplot(PL) +
                        geom_point(aes(x = rt, y = rtiv, color = isov)) +
                        theme_minimal() +
                        scale_color_brewer(palette = "Set2") +
                        xlab("RT(s)"))
    }

    #Carbon number check with M1 peak intensity pattern
    maxpoint <- which.max(PL$rtiv[PL$isov == "M0"])
    rt_at_max <- PL$rt[PL$isov == "M0"][maxpoint]
    if (any(with(PL, isov == "M1" & rt == rt_at_max))) {
        numC <- 100 * PL$rtiv[with(PL, isov == "M1" & rt == rt_at_max)][1]/(1.1 *
            PL$rtiv[PL$isov == "M0"][maxpoint])
    } else {
        numC <- "Could not be defined. M1 not found below M0 max value"
    }

    f <- gsub(x = curSOI$formula, pattern = "[", replacement = "",
        fixed = TRUE) %>% gsub(x = ., pattern = "]", replacement = "",
        fixed = TRUE)

    makeup <- CHNOSZ::makeup(f)
    if (is.numeric(numC)) {
        carbonCheck <- ifelse(round(numC) == makeup[["C"]], "M1 peak intensity pattern matches with the formula",
            "M1 pattern doesn't quite match")
    } else {
        carbonCheck <- numC
    }

    #Generate the isopattern plots
    isotopecode <- data.frame(name = c("13C", "17O", "18O", "2H",
        "15N", "33S", "34S", "36S", "37Cl", "81Br", "41K", "6Li",
        "10B", "21Ne", "22Ne", "25Mg", "26Mg", "29Si", "30Si",
        "42Ca", "43Ca", "44Ca", "48Ca", "46Ti", "47Ti", "49Ti",
        "50Ti", "50Cr", "53Cr", "54Cr", "54Fe", "57Fe", "58Fe",
        "60Ni", "61Ni", "62Ni", "64Ni", "65Cu", "66Zn", "67Zn",
        "68Zn", "70Zn", "76Se", "77Se", "78Se", "82Se", "84Sr",
        "86Sr", "87Sr", "91Zr", "92Zr", "94Zr", "96Zr"), code = c("M",
        "[17O]", "[18O]", "D", "[15N]", "[33S]", "[34S]", "[36S]",
        "[37Cl]", "[81Br]", "[41K]", "[6Li]", "[10B]", "[21Ne]",
        "[22Ne]", "[25Mg]", "[26Mg]", "[29Si]", "[30Si]", "[42Ca]",
        "[43Ca]", "[44Ca]", "[48Ca]", "[46Ti]", "[47Ti]", "[49Ti]",
        "[50Ti]", "[50Cr]", "[53Cr]", "[54Cr]", "[54Fe]", "[57Fe]",
        "[58Fe]", "[60Ni]", "[61Ni]", "[62Ni]", "[64Ni]", "[65Cu]",
        "[66Zn]", "[67Zn]", "[68Zn]", "[70Zn]", "[76Se]", "[77Se]",
        "[78Se]", "[82Se]", "[84Sr]", "[86Sr]", "[87Sr]", "[91Zr]",
        "[92Zr]", "[94Zr]", "[96Zr]"), stringsAsFactors = FALSE)
    data("isotopes")
    # pat <- enviPat::isopattern(isotopes = isotopes, chemforms = f, threshold = 10000/PL$rtiv[PL$isov == 'M0'][maxpoint])[[1]]

    pat <- enviPat::isopattern(isotopes = isotopes, chemforms = f,
        threshold = 0.1, verbose = F)[[1]]


    pat <- rbind(t(apply(pat, 1, function(x) {
        x[3:length(x)] <- x[3:length(x)] - pat[1, 3:ncol(pat)]
        x
    })))
    pat[pat < 0] <- 0
    cols <- which(colnames(pat) %in% isotopecode$name)
    #Sort the cols in the order the iso appear on the list
    idx <- vapply(cols, function(x) {
        which(colnames(pat)[x] == isotopecode$name)[1]
    }, numeric(1))
    cols <- cols[order(idx)]

    pat <- as.data.frame(pat)
    pat$code <- ""
    for (i in cols) {
        different <- which(pat[, i] != 0)
        pat$code[different] <- paste0(pat$code[different], rep(isotopecode$code[colnames(pat)[i] ==
            isotopecode$name], length(different)), pat[different,
            i])
    }
    pat$code[1] <- "M0"
    pat$abundance <- pat$abundance * PL$rtiv[PL$isov == "M0"][maxpoint]/100
    pat$class <- "Theoretical"

    #Calculate isointensities at max M0 intensity peak
    exp_int <- vapply(pat$code, function(iso) {
        x <- PL$rtiv[PL$rt == rt_at_max & PL$isov == iso]
        if (length(x) == 0) {
            return(0)
        } else {
            return(x[1])
        }
    }, numeric(1)) %>% as.data.frame(x = .)
    colnames(exp_int) <- "abundance"
    pat$code <- factor(pat$code, levels = unique(pat$code[order(-pat$abundance)]))  #Sort by abundance

    exp_int$code <- pat$code
    exp_int$class <- "Experimental"

    df <- rbind(pat[, c("abundance", "code", "class")], exp_int)

    if (plot) {
        p2 <- ggplotly(ggplot(df) + geom_col(aes(x = code, y = abundance,
            fill = class), position = "dodge") + scale_fill_brewer(palette = "Set2") +
            ylab("Expected intensities") + theme_minimal() +
            theme(axis.text.x = element_text(angle = 90), plot.margin = unit(c(1,
                1, 1, 1), "cm")))

        # p2 <- p2 %>% style(p2, showlegend = FALSE, traces = 1:nrow(pat))
        p3 <- subplot(p, p2, nrows = 2, which_layout = 2) %>%
            layout(title = curSOI$formula)
    }

    #Calculate isotopic cosine score -We drop M0 to avoid bias-
    tint <- pat$abundance[-1]
    eint <- exp_int$abundance[-1]
    tooWeak <- which(tint < 3000)
    if (length(tooWeak) != 0) {
        tint <- tint[-tooWeak]
        eint <- eint[-tooWeak]
    }
    if (length(tint) == 0) {
        cos <- 1  #Nothing that we could observe. Doesn't get penalized
    } else if (sum(eint) == 0) {
        cos <- 0  #No isotopes detected when there should be
    } else {
        cos <- sum((tint * eint))/(sqrt(sum(tint^2)) * sqrt(sum(eint^2)))
    }

    ##Other atom checks
    toCheck <- c("37Cl", "81Br", "34S")
    otherChecks <- lapply(toCheck, function(atom) {
        if (any(grepl(pattern = atom, x = df$code))) {
            isoname <- paste0("[", atom, "]", 1)
            cur <- df[df$code == isoname, ]
            th <- cur$abundance[cur$class == "Theoretical"]
            if (th < 30000) {
                return(NA)
            }
            exp <- cur$abundance[cur$class == "Experimental"]
            return(data.table::between(exp, th * 0.5, th * 1.5))
        } else {
            return(NA)
        }
    })
    names(otherChecks) <- toCheck

    if (plot) {
        return(list(p3, numC, carbonCheck, cos, otherChecks))
    } else {
        return(list(numC, carbonCheck, cos, otherChecks))
    }
})



setGeneric("coveragePlot", function(struct, entry){standardGeneric("coveragePlot")})
setMethod("coveragePlot", signature = c("RHermesExp", "numeric"), function(struct, entry){
    pl <- struct@data@PL[[entry]]@peaklist[ ,c("rt", "rtiv", "mz")]
    distinct_pl <- nrow(distinct(pl))
    raw <- nrow(struct@data@PL[[entry]]@raw)
    pieplot <- data.frame(Class = c("Covered", "Non-covered"),
                          Value = c(distinct_pl, raw-distinct_pl))
    colors <- c('rgb(211,94,96)', 'rgb(128,133,133)')
    p1 <- plot_ly(data = pieplot, labels = ~Class, values = ~Value,
                  type = "pie",
                  insidetextfont = list(color = '#FFFFFF', size = 18),
                  marker = list(colors = colors,
                                line = list(color = '#FFFFFF', width = 1))
    )
    p1 <- p1 %>% layout(title = "Raw data points covered as PL entries")

    barplot <- data.frame(Class = c("Total entries", "Distinct entries"),
                          Value = c(nrow(pl), distinct_pl))
    p2 <- plot_ly(data = barplot, x = ~Class, y = ~Value,
                  type = "bar",
                  insidetextfont = list(color = '#FFFFFF', size = 18),
                  marker = list(color = colors,
                                line = list(color = '#FFFFFF', width = 1))
    )
    p2 <- p2 %>% layout(title = "Distinct PL entries vs Total Number")

    return(list(p1,p2))
})












