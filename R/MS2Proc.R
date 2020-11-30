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
    #### MS2 Data Importation and Sorting within IL ####------------------------
    validObject(struct)
    MS2Exp <- struct@data@MS2Exp[[id]]
    message("Starting MS/MS data importation, merging and sorting within the IL entries")

    MS2Exp@MS2Data <- MSMSimporter(MS2Exp@IL, MS2files)  #Fills in the MS2data slot
    idx <- vapply(MS2Exp@MS2Data, function(x) {
        return(length(x) != 0)
    }, logical(1))  #Which IL entries are covered by at least 1 scan
    idx <- which(idx)
    if (length(idx) == length(MS2Exp@MS2Data)) {
        message("It seems that all IL entries were covered. Nice!")
    } else {
        message("A total of ", length(MS2Exp@MS2Data) - length(idx)," (",
            round((length(MS2Exp@MS2Data) - length(idx))/
                length(MS2Exp@MS2Data)*100, digits = 2),
            "%) entries were not covered in the experiment")
    }

    if (length(idx) == 0) {
        stop("No entries were covered in the IL. Please check the MS2 file paths and the IL object are correct")
    }

    #### Superspectra Generation ####-------------------------------------------
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

    #### Database Query ####----------------------------------------------------
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

