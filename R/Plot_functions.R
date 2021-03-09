#' @title plotPL
#' @author Roger Gine
#' @family Plots
#' @description Plots the raw data points anotated with a given formula and
#'   adducts.
#' @param struct An RHermesExp object
#' @param id Number of the file to plot
#' @param formula Formula annotation to search for (eg. "C6H12O6")
#' @param ads Adducts to plot. Defaults to NA, which plots all of them by
#'   default.
#' @param rtrange The retention time interval to plot, in seconds (eg.
#'   c(0,1000)). Defaults to a 0-10000s interval, which will cover all points.
#' @param dynamicaxis Whether to use a fixed y scale for all adducts or to adapt
#'   the scale according to each adduct intensity
#' @return An interactive plot_ly object
#' @examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'p <- plotPL(struct, 1, "C5H11NO2", c("M+H", "M+Na", "M+K"), c(80, 120))
#'@export
setGeneric("plotPL", function(struct, id, formula, ads = NA, rtrange = c(0,1e4),
                                dynamicaxis = TRUE) {
    standardGeneric("plotPL")
})

#' @rdname plotPL
setMethod("plotPL", c("RHermesExp", "numeric", "character", "ANY",
    "ANY", "ANY" ),
function(struct, id, formula, ads= NA, rtrange = c(0,1e4), dynamicaxis = TRUE) {
    if(is.na(ads[1])){
        ads <- struct@metadata@ExpParam@adlist$adduct
    } else {
        if(!all(ads %in% struct@metadata@ExpParam@adlist$adduct)){
            ads <- ads %in% struct@metadata@ExpParam@adlist$adduct
            if(length(ads) == 0) {
                stop(paste0("No valid adducts selected, please check your",
                            "selection or leave ads=NA to select them all"))
            }
            warning("Some of the adducts are invalid and will not be included")
        }
    }

    datafile <- struct@data@PL[[id]]@peaklist
    FA_to_ion <- struct@metadata@ExpParam@ionF[[2]]
    setkey(FA_to_ion, "f")

    fs <- FA_to_ion[f == formula, ]
    ions <- fs$ion
    datafile <- filter(datafile, .data$formv %in% ions)
    if (nrow(datafile) == 0) {return()}
    datafile <- filter(datafile, between(.data$rt, rtrange[1], rtrange[2]))
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
                    mapping = aes(x = .data$rt, y = .data$rtiv,
                                    color = .data$isov),
                    size = 0.5, alpha = 0.6) +
        facet_grid(rows = vars(ad)) +
        theme_minimal() + ggtitle(formula)

    return(ggplotly(pl, dynamicTicks = dynamicaxis))
})

#' @title plotCoverage
#' @author Roger Gine
#' @family Plots
#' @description Plots a representation of the raw data points covered by the
#' annotations and the redundancy of those annotations (that is, how many times
#' does the same data point get annotated as two different things).
#' @param struct An RHermesExp object
#' @param id Number of the PL to use in the plot
#' @return A list of two interactive plot_ly objects
#' @examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#' if(FALSE){plotCoverage(struct, 1)}
#'@export
setGeneric("plotCoverage", function(struct, id){
    standardGeneric("plotCoverage")
})

#' @rdname plotCoverage
setMethod("plotCoverage", signature = c("RHermesExp", "numeric"),
function(struct, id) {
    pl <- struct@data@PL[[id]]@peaklist[, c("rt", "rtiv", "mz")]
    distinct_pl <- nrow(distinct(pl))
    noise <- struct@metadata@ExpParam@nthr
    if(nrow(struct@data@PL[[id]]@raw) == 0){
        return(list(plot_ly(), plot_ly()))
    }
    raw <- nrow(filter(struct@data@PL[[id]]@raw, .data$rtiv > noise))
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

#'@title plotSOI
#'@author Roger Gine
#'@family plots
#'@description Plots the SOI data points, the non-SOI points with the same
#'  annotation and, if blank subtraction was performed, it also plots the blank
#'  data points.
#'@param struct An RHermesExp object
#'@param id Number of the SOI list to plot
#'@param formula Formula annotation to search for (eg. "C6H12O6")
#'@param ads Adducts to plot. Defaults to NA, which plots all of them by
#'  default.
#'@param rtrange The retention time interval to plot, in seconds (eg.
#'  c(0,1000)). Defaults to a 0-10000s interval, which will cover all points.
#'@param dynamicaxis Whether to use a fixed y scale for all adducts or to adapt
#'  the scale according to each adduct intensity
#'@param interactive Whether to return a plotly object or a ggplot. Defaults to
#'  TRUE (plotly).
#'@return An interactive plot_ly object or a static ggplot, depending on the
#'  value of interactive
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'p <- plotSOI(struct, 1, "C5H11NO2", c("M+H", "M+Na", "M+K"), c(80, 120))
#'@export
setGeneric("plotSOI", function(struct, id, formula, ads = NA,
                                rtrange = c(0,1e4), dynamicaxis = TRUE,
                                interactive = TRUE) {
    standardGeneric("plotSOI")
})

#' @rdname plotSOI
setMethod("plotSOI", c("RHermesExp", "numeric", "character",
    "ANY", "ANY", "ANY", "ANY"),
function(struct, id, formula, ads = NA, rtrange= c(0, 1e4), dynamicaxis = TRUE,
         interactive = TRUE) {

    #Adduct selection
    if(is.na(ads[1])){
        ads <- struct@metadata@ExpParam@adlist$adduct
    } else {
        if(!all(ads %in% struct@metadata@ExpParam@adlist$adduct)){
            ads <- ads %in% struct@metadata@ExpParam@adlist$adduct
            if(length(ads) == 0) {
                stop(paste0("No valid adducts selected, please check your",
                            "selection or leave ads=NA to select them all"))
            }
            warning("Some of the adducts are invalid and will not be included")
        }
    }

    #Importing SOI list data
    plist <- struct@data@SOI[[id]]@PlotDF
    filen <- struct@data@SOI[[id]]@filename
    plid <- which(vapply(struct@data@PL, function(x) {
        return(x@filename == filen)
    }, logical(1)))[1]

    datafile <- struct@data@PL[[plid]]@peaklist
    datafile <- datafile[.data$isov == "M0", ]
    Class <- NULL #To appease R CMD Check "no visible binding"
    datafile[, Class := "Sample"]

    #Importing blank data if blank subtraction was performed on the SOI list
    if(struct@data@SOI[[id]]@SOIParam@blanksub){
        blankid <- which(struct@metadata@filenames ==
                            struct@data@SOI[[id]]@SOIParam@blankname)
        if(length(blankid) > 1) blankid <- blankid[1]
        blankfile <- struct@data@PL[[blankid]]@peaklist
        blankfile <- blankfile[.data$isov == "M0", ]
        blankfile[, Class := "Blank"]
        datafile <- rbind(datafile, blankfile)
    } else {blankid <- NA}

    #Filter by selected adducts and RT interval
    FA_to_ion <- struct@metadata@ExpParam@ionF[[2]]
    setkeyv(FA_to_ion, c("f"))
    fs <- FA_to_ion[f == formula, ]
    ions <- fs$ion
    datafile <- filter(datafile, .data$formv %in% ions)
    datafile <- filter(datafile, between(.data$rt, rtrange[1], rtrange[2]))
    soiinfo <- filter(plist, .data$form %in% ions)
    soiinfo <- filter(soiinfo, between(.data$rt, rtrange[1], rtrange[2]))

    #No SOI data in the RT interval
    if (nrow(soiinfo) == 0) {return()}

    names(soiinfo)[names(soiinfo) == "form"] <- "formv"
    names(soiinfo)[names(soiinfo) == "rtiv"] <- "Intensity"
    names(datafile)[names(datafile) == "rtiv"] <- "Intensity"

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

    datafile <- filter(datafile, .data$ad %in% ads)
    soiinfo <- filter(soiinfo, .data$ad %in% ads)

    if (nrow(soiinfo) == 0) {return()}
    soiinfo$Class <- "Sample-SOI"
    if(!is.na(blankid)){
        plot <- ggplot() +
                geom_point(data = datafile[datafile$Class == "Sample", ],
                        mapping = aes(x = .data$rt, y = .data$Intensity,
                                        color = .data$Class),
                        alpha = 0.4) +
                geom_point(data = datafile[datafile$Class == "Blank", ],
                        mapping = aes(x = .data$rt, y = .data$Intensity,
                                        color = .data$Class),
                        alpha = 0.5)+
                geom_point(data = soiinfo,
                        mapping = aes(x = .data$rt, y = .data$Intensity,
                                        color = .data$Class),
                        alpha = 0.8) +
                scale_color_manual(breaks = c("Blank", "Sample", "Sample-SOI"),
                                values = c("#6D9503", "#8E032B", "#370B6B"))+
                facet_grid(rows = vars(ad)) + theme_minimal() +
                ggtitle(formula)
    } else {
        plot <- ggplot() +
                geom_point(data = datafile[datafile$Class == "Sample", ],
                        mapping = aes(x = .data$rt, y = .data$Intensity,
                                        color = .data$Class),
                        alpha = 0.3)+
                geom_point(data = soiinfo,
                        mapping = aes(x = .data$rt, y = .data$Intensity,
                                        color = .data$Class),
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


#'@title plotFidelity
#'@author Roger Gine
#'@family plots
#'@description Plots the selected SOI isotopic profile and a
#'  comparison with its theoretical abundances calculated from
#'  the molecular formula.
#'@param struct An RHermesExp object.
#'@param id Number of the SOI to plot.
#'@param entry The SOI entry to check.
#'@param plot Default to TRUE. The parameter is used for
#'  consistency with the function internal use in filterSOI(). If
#'  set to FALSE, returns some statistics about the isotopic
#'  fidelity.
#'@return An interactive plot_ly object. If plot set to FALSE,
#'  returns a list of isotopic fidelity metrics.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'p <- plotFidelity(struct, 1, 9)
#'@export
setGeneric("plotFidelity", function(struct, id, entry, plot = TRUE) {
    standardGeneric("plotFidelity")
})

#' @rdname plotFidelity
setMethod("plotFidelity", c("RHermesExp", "numeric", "numeric", "ANY"),
function(struct, id, entry, plot = TRUE) {
    #Extract SOI and PL information from the selected SOI list
    SOI <- struct@data@SOI[[id]]
    fname <- SOI@filename
    correspondingPL <- which(struct@metadata@filenames == fname)
    PL <- struct@data@PL[[correspondingPL]]@peaklist

    #Filter the PL to the selected SOI region
    curSOI <- SOI@SOIList[entry, ]
    PL <- filter(PL, .data$formv == curSOI$formula)
    PL <- filter(PL, data.table::between(.data$rt, curSOI$start, curSOI$end))

    if (plot) {
        p <- ggplotly(ggplot(PL) +
                        geom_point(aes(x = .data$rt, y = .data$rtiv,
                                        color = .data$isov)) +
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

    isotopes <- NULL #To appease R CMD Check "no visible binding"
    data("isotopes", package = "enviPat", envir = environment())

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
            geom_col(aes(x = .data$code, y = .data$abundance,
                            fill = .data$class),
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


#'@title plotIL
#'@author Roger Gine
#'@description Draws a mz-rt representation of the inclusion list entries,
#'coloured by their intensity value.
#'@param struct An RHermesExp object
#'@param id Number of the inclusion list to plot
#'@return An interactive plot_ly object
#'@family Plots
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'plotIL(struct, 1)
#'@export
setGeneric("plotIL", function(struct, id) {
    standardGeneric("plotIL")
})
#' @rdname plotIL
setMethod("plotIL", c("RHermesExp", "numeric"),
    function(struct, id){
        ggplotly(ggplot(struct@data@MS2Exp[[id]]@IL@IL) +
                geom_segment(aes(x = .data[["start"]], xend = .data[["end"]],
                                y = .data[["mass"]], yend = .data[["mass"]],
                                color = log10(.data[["MaxInt"]]))))
})



#'@title plotSS
#'@author Roger Gine
#'@description Plots an Hermes-cleaned MS2 spectum (mz/int).
#'@param struct An RHermesExp object
#'@param ms2id Number of the MS2Exp object where the spectrum is
#'@param ssnumber Number of the spectrum to plot
#'@return A ggplot object
#'@family plots
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'if(FALSE){plotSS(struct, 1, 1)}
#'@export
setGeneric("plotSS", function(struct, ms2id, ssnumber) {
    standardGeneric("plotSS")
})

#' @rdname plotSS
setMethod("plotSS", c("RHermesExp", "numeric", "numeric"),
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
        geom_segment(data = query, aes(x = .data$mz, xend = .data$mz,
                                        y = 0, yend = .data$int),
                        color = "black") +
        geom_point(data = moldf, aes(x = .data$mz, y = 0), shape = 17,
                    size = 2) +
        theme_minimal() + ylab("Intensity (Counts)") +
        theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
            text = element_text(size = 11, family = "Segoe UI Light"),
            plot.title = element_text(hjust = 0.5)) +
        geom_text(data = bestdf, aes(x = .data$mz, y = .data$int + 5,
                                        label = .data$mz),
                family = "Segoe UI Light", check_overlap = TRUE) +
        scale_x_continuous(limits = c(min(query$mz, molecmass) - 20,
                                    max(query$mz, molecmass) + 20))+
        ggtitle(title, subtitle)

    # ggplotly(pl, height = 400)
    pl
})

#'@title plotRawMS2
#'@author Roger Gine
#'@description Plots an all MS2 fragments distribution along the rt of a given
#'  inclusion list entry.
#'@param struct An RHermesExp object
#'@param ms2id Number of the MS2Exp object where the spectrum is
#'@param entryid Number of the inclusion list entry to plot
#'@return A list of four objects: a plot of the fragments grouped by their mz, a
#'  plot of the fragments grouped by their shape similarity, the similarity
#'  network and a summary table of peaks found and the group they belong to.
#'@family plots
#'@export
#'@import networkD3
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'plotRawMS2(struct, 1, 2)
setGeneric("plotRawMS2", function(struct, ms2id, entryid) {
    standardGeneric("plotRawMS2")
})

#' @rdname plotRawMS2
setMethod("plotRawMS2", c("RHermesExp", "numeric", "ANY"),
function(struct, ms2id, entryid) {
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

    # oldpoints <- soi
    # oldpoints$member <- "Original points"
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

    p_bymz <- ggplot(soi) +
                geom_point(aes(x = .data$rt, y = .data$rtiv,
                                color = as.factor(.data$mz)))+
                xlim(xlim) + theme_minimal()
    p_bygroup <- ggplot(soi) +
                geom_point(aes(x = .data$rt, y = .data$rtiv,
                                color = as.factor(.data$member)))+
                xlim(xlim) + theme_minimal()

    net <- networkD3::igraph_to_networkD3(net, group = members)

    net <- visNetwork(net$nodes %>% rename(label = .data$name) %>%
        mutate(id = seq_len(nrow(net$nodes)) - 1), net$links %>%
        rename(from = .data$source, to = .data$target))

    net %<>% visNodes(color = list(background = "lightblue")) %>%
        visEdges(smooth = FALSE) %>%
        visPhysics(solver = "forceAtlas2Based",
            forceAtlas2Based = list(gravitationalConstant = -100),
            stabilization = FALSE)
    return(list(p_bymz = p_bymz, p_bygroup = p_bygroup, net = net, pks = pks))
})

