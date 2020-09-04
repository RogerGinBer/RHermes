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
#'@param plotting Logical, whether you want mirror plots to be generated.
#'@param outpdf Character, where you want to store the plots (see plotting).
#'@param outcsv Character, where to store a summary csv of the identifications.
#'@return An RHermesExp object with the identifications set in the used MS2Exp slot.
#'@examples
#'MS2Proc(myHermes, 1, c('C:/myFolder/File1.mzML', 'C:/myFolder/File2.mzML'))
#'
#' @export
setGeneric("MS2Proc", function(struct, id, MS2files,
    referenceDB = "D:/sp_MassBankEU_20200316_203615.RData",
    mincos = 0.6, plotting = FALSE, outpdf = "mirror.pdf",
    outcsv = "./results.csv") {
    standardGeneric("MS2Proc")
})
setMethod("MS2Proc", c("RHermesExp", "numeric", "character",
    "ANY", "ANY", "ANY", "ANY", "ANY"), function(struct, id,
    MS2files, referenceDB = "D:/sp_MassBankEU_20200316_203615.RData",
    mincos = 0.6) {
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
    purifiedSpectra <- CliqueMSMS(MS2Exp, idx)

    #### Database Query ####----------------------------------------------------

    # Catching possible user errors
    tryCatch({
        load(referenceDB)
    }, error = function(cond) {
        stop("The referenceDB input isn't valid")
    })

    # Data.tables for faster indexing
    metaesp <- as.data.table(df_metaespectres)
    espmet <- as.data.table(df_EspectreMetabolit)
    metamet <- as.data.table(df_metametabolits)

    metamet$formula <- unlist(metamet$formula)

    setkeyv(metaesp, c("idespectre"))
    setkeyv(espmet, c("idmetabolit"))
    setkeyv(metamet, c("formula"))

    metaesp$polaritat <- toupper(metaesp$polaritat)
    metamet$nommetabolit[is.na(metamet$nommetabolit)] <- metamet$casname[is.na(metamet$nommetabolit)]
    polarity <- ifelse(struct@metadata@ExpParam@ion == "+", "POSITIVE",
        "NEGATIVE")


    allf <- purifiedSpectra$anot %>% unlist() %>% unique()

    # Retrieving spectra
    message("Retrieving MS2 spectra from the reference database")
    retrievedMSMS <- lapply(allf, function(f) {
            idmet <- metamet[.(f), ]$idmetabolit
            if (is.na(idmet[1])) {
                return(list())
            }
            RES <- lapply(idmet, function(id) {
                idesp <- espmet[.(id), ]$idespectre
                idesp <- idesp[metaesp[.(idesp), ]$polaritat == polarity]
                return(list_fragments[["espectre"]][list_fragments[["idespectre"]] %in%
                  idesp])
            })
            novalidspec <- vapply(RES, function(x) {
                length(x) != 0
            }, logical(1))
            if (!any(novalidspec)) {
                RES <- RES[-novalidspec]
                if (length(RES) == 0) {
                  return(list())
                }
                names(RES) <- paste(rep(f, times = nrow(metamet[.(f),
                  ]) - length(novalidspec)),
                  metamet[.(f), ]$idmetabolit[-novalidspec],
                  lapply(metamet[.(f), ]$nommetabolit, function(x) {
                    x[[1]]
                  }) %>% unlist()[-novalidspec], metamet[.(f),
                    ]$smiles[-novalidspec], sep = "#")

            } else {
                names(RES) <- paste(rep(f, times = nrow(metamet[.(f),])),
                                    metamet[.(f), ]$idmetabolit,
                                    lapply(metamet[.(f), ]$nommetabolit,
                                        function(x) {x[[1]]}
                                    ) %>% unlist(),
                                    metamet[.(f), ]$smiles,
                                    sep = "#")
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
        cos <- rapply(curspec, function(compound) {
            MSMScosineSim(entry, t(compound), minhits = 1, mode = "full")
        }, classes = "matrix", how = "replace")
        cos[which(lapply(cos, length) != 0)]
    }, purifiedSpectra$ssdata, corresponding)

    #### Output formatting ####------------------------------------------------
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
            goodf <- which(vapply(x[[f]], function(coslist) {
                any(vapply(coslist, function(cos) {
                  cos > mincos
                }, logical(1)))
            }, logical(1)))
            fnames <- names(x[[f]][goodf])
            resdf <- data.frame(formula = fnames, cos = numeric(length(fnames)),
                id = numeric(length(fnames)),
                stringsAsFactors = FALSE)
            for (i in seq_along(goodf)) {
                coslist <- x[[f]][[goodf[i]]]
                resdf$id[i] <- which.max(vapply(coslist, function(cos) {
                  cos[1]
                }, numeric(1)))
                resdf$cos[i] <- coslist[[resdf$id[i]]]
            }
            return(resdf)

        }) %>% do.call(rbind, .)
        return(RES)
    })
    MS2Exp@Ident <- list(MS2Features = purifiedSpectra,
                         DatabaseSpectra = retrievedMSMS,
                         MS2_correspondance = corresponding)
    message("Done!")
    struct@data@MS2Exp[[id]] <- MS2Exp
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
