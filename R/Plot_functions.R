#'@export
setGeneric("PLPlot", function(struct, id, formula, rtrange, dynamicaxis, ads) {
    standardGeneric("PLPlot")
})
setMethod("PLPlot", c("RHermesExp", "numeric", "character",
    "numeric", "logical", "character"),
function(struct, id, formula, rtrange, dynamicaxis, ads) {
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
        ad <- as.character(fs$an[fs$ion == f])
        datafile$ad[datafile$formv == f] <- ad
    }
    datafile <- filter(datafile, ad %in% ads)
    if (nrow(datafile) == 0) {return()}

    pl <- ggplot() +
        geom_point(data = datafile,
                    mapping = aes(x = rt, y = rtiv, color = isov),
                    size = 0.5, alpha = 0.6) +
        facet_grid(rows = vars(ad)) +
        theme_minimal() + ggtitle(formula)

    return(ggplotly(pl, dynamicTicks = dynamicaxis))

})


#'@export
setGeneric("coveragePlot", function(struct, entry){
    standardGeneric("coveragePlot")
})
setMethod("coveragePlot", signature = c("RHermesExp", "numeric"),
function(struct, entry) {
    pl <- struct@data@PL[[entry]]@peaklist[, c("rt", "rtiv", "mz")]
    distinct_pl <- nrow(distinct(pl))
    noise <- struct@metadata@ExpParam@nthr
    raw <- nrow(filter(struct@data@PL[[entry]]@raw, .data$rtiv > noise))
    pieplot <- data.frame(Class = c("Covered", "Non-covered"),
                            Value = c(distinct_pl, raw - distinct_pl))
    colors <- c("rgb(211,94,96)", "rgb(128,133,133)")

    p1 <- plot_ly(data = pieplot, labels = ~Class, values = ~Value,
                    type = "pie",
                    insidetextfont = list(color = "#FFFFFF", size = 18),
                    marker = list(colors = colors,
                                line = list(color = "#FFFFFF", width = 1)))
    p1 <- p1 %>% layout(title = "Raw data points covered as PL entries")

    barplot <- data.frame(Class = c("Total entries", "Distinct entries"),
                            Value = c(nrow(pl), distinct_pl))
    p2 <- plot_ly(data = barplot, x = ~Class, y = ~Value, type = "bar",
                    insidetextfont = list(color = "#FFFFFF", size = 18),
                    marker = list(color = colors,
                                line = list(color = "#FFFFFF", width = 1)))
    p2 <- p2 %>% layout(title = "Distinct PL entries vs Total Number")

    return(list(p1, p2))
})


#'@export
setGeneric("SoiPlot", function(struct, id, formula, rtrange,
                                dynamicaxis = TRUE, ads = NA, blankid = NA,
                                interactive = TRUE) {
    standardGeneric("SoiPlot")
})
setMethod("SoiPlot", c("RHermesExp", "numeric", "character",
    "numeric", "ANY", "ANY", "ANY", "ANY"),
function(struct, id, formula, rtrange, dynamicaxis = TRUE, ads = NA,
            blankid = NA, interactive = TRUE) {

    if(is.na(ads[1])){
        ads <- struct@metadata@ExpParam@adlist$adduct
    }

    plist <- struct@data@SOI[[id]]@PlotDF
    filen <- struct@data@SOI[[id]]@filename
    plid <- which(vapply(struct@data@PL, function(x) {
        return(x@filename == filen)
    }, logical(1)))[1]

    datafile <- struct@data@PL[[plid]]@peaklist
    datafile <- datafile[isov == "M0", ]
    datafile[, Class := "Sample"]

    if(!is.na(blankid)){
        blankfile <- struct@data@PL[[blankid]]@peaklist
        blankfile <- blankfile[isov == "M0", ]
        blankfile[, Class := "Blank"]
        datafile <- rbindlist(list(datafile, blankfile))
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
        ad <- as.character(fs$an[fs$ion == f])
        datafile$ad[datafile$formv == f] <- ad
    }

    soiinfo$ad <- ""
    for (f in unique(soiinfo$formv)) {
        ad <- as.character(fs$an[fs$ion == f])
        soiinfo$ad[soiinfo$formv == f] <- ad
    }

    datafile <- filter(datafile, ad %in% ads)
    soiinfo <- filter(soiinfo, ad %in% ads)

    if (nrow(soiinfo) == 0) {return()}
    soiinfo$Class <- "Sample-SOI"
    if(!is.na(blankid)){
        plot <- ggplot() +
                geom_point(data = datafile[datafile$Class == "Sample", ],
                        mapping = aes(x = rt, y = Intensity, color = Class),
                        alpha = 0.4) +
                geom_point(data = datafile[datafile$Class == "Blank", ],
                        mapping = aes(x = rt, y = Intensity, color = Class),
                        alpha = 0.5)+
                geom_point(data = soiinfo,
                        mapping = aes(x = rt, y = Intensity, color = Class),
                        alpha = 0.8) +
                scale_color_manual(breaks = c("Blank", "Sample", "Sample-SOI"),
                                values = c("#6D9503", "#8E032B", "#370B6B"))+
                facet_grid(rows = vars(ad)) + theme_minimal() +
                ggtitle(formula)
    } else {
        plot <- ggplot() +
                geom_point(data = datafile[datafile$Class == "Sample", ],
                        mapping = aes(x = rt, y = Intensity, color = Class),
                        alpha = 0.3)+
                geom_point(data = soiinfo,
                        mapping = aes(x = rt, y = Intensity, color = Class),
                        alpha = 0.8) +
                scale_color_manual(breaks = c("Sample", "Sample-SOI"),
                                values = c("#8E032B", "#370B6B"))+
                facet_grid(rows = vars(ad)) + theme_minimal() +
                ggtitle(formula)
    }
    if(interactive){
        return(ggplotly(plot, dynamicTicks = dynamicaxis))
    } else {
        return(plot)
    }
})

#'@export
setGeneric("IsoFidelity", function(struct, soilist, entry,
                                    plot = TRUE) {
    standardGeneric("IsoFidelity")
})
setMethod("IsoFidelity", c("RHermesExp", "numeric", "numeric", "ANY"),
function(struct, soilist, entry, plot = TRUE) {
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
                        xlab("RT(s)"))
    }

    #Carbon number check with M1 peak intensity pattern
    maxpoint <- which.max(PL$rtiv[PL$isov == "M0"])
    rt_at_max <- PL$rt[PL$isov == "M0"][maxpoint]
    if (any(with(PL, isov == "M1" & rt == rt_at_max))) {
        numC <- 100 * PL$rtiv[with(PL, isov == "M1" & rt == rt_at_max)][1] /
            (1.1 * PL$rtiv[PL$isov == "M0"][maxpoint])
    } else {
        numC <- "Could not be defined. M1 not found below M0 max value"
    }

    f <- gsub(x = curSOI$formula, pattern = "[", replacement = "",
                fixed = TRUE) %>%
                gsub(x = ., pattern = "]", replacement = "", fixed = TRUE)

    makeup <- CHNOSZ::makeup(f)
    if (is.numeric(numC)) {
        carbonCheck <- ifelse(round(numC) == makeup[["C"]],
                            "M1 intensity pattern matches with the formula",
                            "M1 pattern doesn't quite match")
    } else {
        carbonCheck <- numC
    }

    #Generate the isopattern plots
    isodf <- data.frame(
        name = c("13C", "17O", "18O", "2H",
            "15N", "33S", "34S", "36S", "37Cl", "81Br", "41K", "6Li",
            "10B", "21Ne", "22Ne", "25Mg", "26Mg", "29Si", "30Si",
            "42Ca", "43Ca", "44Ca", "48Ca", "46Ti", "47Ti", "49Ti",
            "50Ti", "50Cr", "53Cr", "54Cr", "54Fe", "57Fe", "58Fe",
            "60Ni", "61Ni", "62Ni", "64Ni", "65Cu", "66Zn", "67Zn",
            "68Zn", "70Zn", "76Se", "77Se", "78Se", "82Se", "84Sr",
            "86Sr", "87Sr", "91Zr", "92Zr", "94Zr", "96Zr"),
        code = c("M",
            "[17O]", "[18O]", "D", "[15N]", "[33S]", "[34S]", "[36S]",
            "[37Cl]", "[81Br]", "[41K]", "[6Li]", "[10B]", "[21Ne]",
            "[22Ne]", "[25Mg]", "[26Mg]", "[29Si]", "[30Si]", "[42Ca]",
            "[43Ca]", "[44Ca]", "[48Ca]", "[46Ti]", "[47Ti]", "[49Ti]",
            "[50Ti]", "[50Cr]", "[53Cr]", "[54Cr]", "[54Fe]", "[57Fe]",
            "[58Fe]", "[60Ni]", "[61Ni]", "[62Ni]", "[64Ni]", "[65Cu]",
            "[66Zn]", "[67Zn]", "[68Zn]", "[70Zn]", "[76Se]", "[77Se]",
            "[78Se]", "[82Se]", "[84Sr]", "[86Sr]", "[87Sr]", "[91Zr]",
            "[92Zr]", "[94Zr]", "[96Zr]"),
        stringsAsFactors = FALSE)

    data("isotopes", package = "enviPat")

    pat <- enviPat::isopattern(isotopes = isotopes, chemforms = f,
                            threshold = 0.1, verbose = F)[[1]]


    pat <- apply(pat, 1, function(x){
        x[3:length(x)] <- x[3:length(x)] - pat[1, 3:ncol(pat)]
        x
    }) %>% t() %>% rbind()
    pat[pat < 0] <- 0
    cols <- which(colnames(pat) %in% isodf$name)

    #Sort the cols in the order the iso appear on the list
    idx <- vapply(cols, function(x) {
        which(colnames(pat)[x] == isodf$name)[1]
    }, numeric(1))
    cols <- cols[order(idx)]

    pat <- as.data.frame(pat)
    pat$code <- ""
    for (i in cols) {
        different <- which(pat[, i] != 0)
        pat$code[different] <- paste0(pat$code[different],
                                        rep(isodf$code[colnames(pat)[i] ==
                                                            isodf$name],
                                        length(different)),
                                        pat[different, i])
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

    #Sort by abundance
    colnames(exp_int) <- "abundance"
    pat$code <- factor(pat$code,
                        levels = unique(pat$code[order(-pat$abundance)]))

    exp_int$code <- pat$code
    exp_int$class <- "Experimental"

    df <- rbind(pat[, c("abundance", "code", "class")], exp_int)

    if (plot) {
    p2 <- ggplot(df)+
            geom_col(aes(x = code, y = abundance, fill = class),
                        position = "dodge") +
            ylab("Expected intensities") + theme_minimal() +
            theme(axis.text.x = element_text(angle = 90),
                    plot.margin = unit(c(1, 1, 1, 1), "cm"))
    p2 <- ggplotly(p2)
    p3 <- subplot(p, p2, nrows = 2, which_layout = 2) %>%
            layout(title = curSOI$formula)
    }

    #Calculate isotopic cosine score -We drop M0 to avoid bias-
    tint <- pat$abundance[-1]
    eint <- exp_int$abundance[-1]
    tooWeak <- which(tint < 20000)
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

#'@export
setGeneric("MirrorPlot", function(struct, ms2id, ssnumber, patform,
                                    mode = "hits") {
    standardGeneric("MirrorPlot")
})
setMethod("MirrorPlot", c("RHermesExp", "numeric", "numeric"),
function(struct, ms2id, ssnumber, patform, mode = "hits") {
    entry <- struct@data@MS2Exp[[ms2id]]@Ident$MS2Features[ssnumber,]
    query <- entry$ssdata[[1]]
    maxint <- max(query$int)
    query$int <- query$int/max(query$int) * 100
    bestdf <- query[query$int > 10, ]
    bestdf$mz <- round(bestdf$mz, 4)
    molecmass <- entry$precmass
    baseline <- 1000
    subtitle <- ""
    title <- ""
    baseline <- baseline/maxint * 100

    f <- list(
        family = "Open Sans",
        size = 16,
        color = "black"
    )


    if(mode == "hits"){
        ref <- struct@data@MS2Exp[[ms2id]]@Ident$MS2_correspondance[[ssnumber]]
        pattern <- struct@data@MS2Exp[[ms2id]]@Ident$DatabaseSpectra[ref]
        pattern <- unname(pattern) #Avoids "name.subname" when unlisting
        pattern <- unlist(pattern, recursive = FALSE, use.names = TRUE)
    }

    mirrplot <- lapply(patform, function(x){
        pl <- ggplot()
        if(mode == "hits"){
            index <- entry$results[[1]]$id[entry$results[[1]]$formula == x]
            if(length(index) == 0){
                refspec <- NULL
            } else {
                refspec <- pattern[[x]][[index]]
                spec_energy <- names(pattern[[x]])[index]
                name <- strsplit(patform, split = "#")[[patform == x]][[3]]
                a <- list(
                    text = paste(name, spec_energy),
                    font = f,
                    xref = "paper",
                    yref = "paper",
                    yanchor = "bottom",
                    xanchor = "center",
                    align = "center",
                    x = 0.5,
                    y = 1,
                    showarrow = FALSE
                )
            }
            if(is.null(refspec)){return(ggplotly(pl))}
            refspec <- as.data.frame(t(refspec))
            pl <- pl + scale_x_continuous(limits = c(min(c(refspec$mz,
                                                    query$mz, molecmass)) - 20,
                                                    max(c(refspec$mz, query$mz,
                                                        molecmass)) + 20))
        }
        if(mode == "versus"){
            refmass<-struct@data@MS2Exp[[ms2id]]@Ident$MS2Features$precmass[[x]]
            refspec <- struct@data@MS2Exp[[ms2id]]@Ident$MS2Features$ssdata[[x]]
            if(is.null(refspec)){return(ggplotly(pl))}

            a <- list(
                text = paste("SS comparison between", ssnumber, "and", x),
                font = f,
                xref = "paper",
                yref = "paper",
                yanchor = "bottom",
                xanchor = "center",
                align = "center",
                x = 0.5,
                y = 1,
                showarrow = FALSE
            )

            refdf <- data.frame(mz = refmass)
            pl <- pl +
                geom_point(data = refdf, aes(x = mz, y = 0), shape = 25,
                            size = 4, color = "red", fill = "red") +
                scale_x_continuous(limits = c(min(c(refspec$mz, query$mz,
                                                    molecmass, refmass)) - 20,
                                    max(c(refspec$mz, query$mz, molecmass,
                                            refmass)) + 20))
        }
            colnames(refspec) <- c("mz", "int")
            refspec$int <- refspec$int/max(refspec$int) * 100

            moldf <- data.frame(mz = molecmass)
            bldf <- data.frame(xmin = min(c(refspec$mz, query$mz,
                                        molecmass)) - 5,
                            xmax = max(c(refspec$mz,query$mz,
                                        molecmass) + 5), y = baseline)

            pl <- pl + geom_segment(data = query, aes(x = mz,
                        xend = mz, y = 0, yend = int), color = "black") +
                    geom_segment(data = refspec, aes(x = mz, xend = mz,
                        y = 0, yend = -int), color = "red") +
                    geom_segment(data = bldf, aes(x = xmin, xend = xmax,
                        y = y, yend = y), linetype = "dashed", color = "black",
                        alpha = 0.3) +
                    geom_segment(data = bldf, aes(x = xmin, xend = xmax,
                        y = -y, yend = -y), linetype = "dashed", color = "red",
                        alpha = 0.3) +
                    geom_point(data = moldf, aes(x = mz, y = 0), shape = 17,
                        size = 2) +
                    theme_minimal() + ylab("% Intensity") +
                    theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
                            text = element_text(size = 11,
                                                family = "Segoe UI Light"),
                            plot.title = element_text(hjust = 0.5)) +
                    geom_text(data = bestdf, aes(x = mz, y = int + 5,
                                                    label = mz),
                            family = "Segoe UI Light", check_overlap = TRUE)

            base_height <- ifelse(length(patform)<5, 850/length(patform), 200)
            ggplotly(pl, height = base_height * length(patform)) %>%
                layout(annotations = a)
        })

    return(subplot(mirrplot, nrows = length(mirrplot), shareX = TRUE,
                    which_layout = 1))
})

#'@export
ILplot <- function(struct, ILnumber){
    ggplotly(ggplot(struct@data@MS2Exp[[ILnumber]]@IL@IL) +
            geom_segment(aes(x = .data[["start"]], xend = .data[["end"]],
                            y = .data[["mass"]], yend = .data[["mass"]],
                            color = log10(.data[["MaxInt"]]))))
}

#'@export
setGeneric("SSPlot", function(struct, ms2id, ssnumber) {
    standardGeneric("SSPlot")
})
setMethod("SSPlot", c("RHermesExp", "numeric", "numeric"),
function(struct, ms2id, ssnumber) {
    entry <- struct@data@MS2Exp[[ms2id]]@Ident$MS2Features[ssnumber,]
    if(!is.data.frame(entry)){return(ggplotly(ggplot()))}
    query <- entry$ssdata[[1]]
    maxint <- max(query$int)
    #query$int <- query$int/max(query$int) * 100 #No normalization for this plot
    bestdf <- query[query$int > 10, ]
    bestdf$mz <- round(bestdf$mz, 4)
    molecmass <- entry$precmass
    moldf <- data.frame(mz = molecmass)

    subtitle <- ""
    title <- ""

    pl <- ggplot() +
        geom_segment(data = query, aes(x = mz, xend = mz, y = 0, yend = int),
                        color = "black") +
        geom_point(data = moldf, aes(x = mz, y = 0), shape = 17, size = 2) +
        theme_minimal() + ylab("Intensity (Counts)") +
        theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
            text = element_text(size = 11, family = "Segoe UI Light"),
            plot.title = element_text(hjust = 0.5)) +
        geom_text(data = bestdf, aes(x = mz, y = int + 5, label = mz),
                family = "Segoe UI Light", check_overlap = TRUE) +
        scale_x_continuous(limits = c(min(query$mz, molecmass) - 20,
                                    max(query$mz, molecmass) + 20))+
        ggtitle(title, subtitle)

    # ggplotly(pl, height = 400)
    pl
})


#'@export
#'@import networkD3
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
    ss <- generate_ss(entryid, MS2list = ms2data, contaminant = 173.5,
                        delta = 0.1, fs = character(), idx = numeric(),
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
    if (nrow(soi) == 0) {
        return(list(p_bymz = ggplotly(ggplot()), p_bygroup = ggplotly(ggplot()),
                    net = NA, pks = data.frame()))
    }

    if (any(!mem_to_keep)) {
        net <- igraph::delete.vertices(net, which(!mem_to_keep))
    }
    soi$member <- "Not considered"
    for (i in unique(members)) {
        soi$member[soi$peak %in% which(members == i)] <- paste("Superspec.",
                                                            as.character(i))
    }

    p_bymz <- ggplot(rbind(soi, oldpoints)) +
                geom_point(aes(x = rt, y = rtiv, color = as.factor(mz)))+
                xlim(xlim)
    p_bygroup <- ggplot(rbind(soi, oldpoints)) +
                geom_point(aes(x = rt, y = rtiv, color = as.factor(member)))+
                xlim(xlim)

    net <- networkD3::igraph_to_networkD3(net, group = members)

    net <- visNetwork(net$nodes %>% rename(label = name) %>%
        mutate(id = seq_len(nrow(net$nodes)) - 1), net$links %>%
        rename(from = source, to = target))

    net %<>% visNodes(color = list(background = "lightblue")) %>%
        visEdges(smooth = FALSE) %>%
        visPhysics(solver = "forceAtlas2Based",
            forceAtlas2Based = list(gravitationalConstant = -100),
            stabilization = FALSE)
    return(list(p_bymz = p_bymz, p_bygroup = p_bygroup, net = net, pks = pks))
})

