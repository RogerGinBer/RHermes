#'@title generateIL
#'@description Groups the SOIs generated in the previous step and sets them up
#'  in an inclusion list ready for the MS/MS analysis. It takes into account the
#'  instrument precursor filter mz range and groups the SOIs accordingly. The
#'  user can also specify whether they want to use all SOIs or prioritize SOIs
#'  with a certain annotation (M+H adducts, prioritize SOIs with adequate
#'  isotopic patterns, etc.).
#'@param struct The RHermesExp object in which to save the resulting IL
#'@param id The SOI list ID from which to generate the IL
#'@param par An inclusion list parameter object, generated with ILParam()
#'@return An RHermes object with the IL inside a new MS2Exp entry, inside the
#'  data slot
#'@examples
#'if(FALSE){
#' par <- ILParam(filtermz = 0.5, priorization = "full")
#' myHermes <- generateIL(myHermes, 1, par)
#'}
#'@importFrom dplyr between
#'@export
setGeneric("generateIL", function(struct, id, par) {
    standardGeneric("generateIL")
})
#' @rdname generateIL
setMethod("generateIL", c("RHermesExp", "numeric", "ANY"),
    function(struct, id, par = ILParam()) {
        validObject(struct)
        validObject(par)
        if (!is(par, "ILParam")) {
            stop("Please enter a valid IL parameter set with ILParam()")
        }
        if (id > length(struct@data@SOI)) {
            stop("The entered ID is greater than the number of SOI lists")
        }
        message("Now processing the IL:")
        struct@data@MS2Exp <- c(struct@data@MS2Exp,
                                RHermesMS2Exp(
                                    IL = inclusionList(struct, par, id),
                                    MS2Data = list(),
                                    Ident = list()
                                ))
        message("Done!")
        return(struct)
})

inclusionList <- function(struct, params, id) {
    SoiList <- struct@data@SOI[[id]]
    ppm <- struct@metadata@ExpParam@ppm
    oSoi <- SoiList@SOIList
    rtmargin <- params@rtmargin
    priorization <- params@priorization
    adduct <- params@ad
    filtermz <- params@filtermz
    filterrt <- params@filterrt
    if (priorization == "yes") {
        GL <- GLprior(oSoi, adduct, rtmargin, ppm)
        low <- which(GL$MaxInt < 50000)
        rare <- which(vapply(GL$ad[low], function(x) {
            !any(unlist(x) %in% c("M+H", "M+NH4", "M+Na", "M+K", "M+", "M-H",
                                  "M+Cl", "M+Br", "M+H2O-H"))
        }, logical(1)))
        if (length(rare) != 0) {GL <- GL[-low[rare], ]}
        GL <- GL[, c("start", "end", "formula", "mass", "MaxInt", "anot")]
        GL$originalSOI <- unlist(sapply(seq(nrow(GL)),function(j){
            y <- GL[j,]
            which(oSoi$start==y$start &
                      oSoi$end==y$end &
                      oSoi$mass==y$mass &
                      oSoi$formula==y$formula)
        }))
    }else{
        GL <- oSoi[, c("start", "end", "formula", "mass", "MaxInt", "anot")]
        GL$originalSOI <- seq_len(nrow(GL))
    }
    # GL$originalSOI <- seq_len(nrow(GL))

    GL <- lapply(seq_len(nrow(GL)), function(x) {
        cur <- GL[x, ]
        df <- as.data.frame(matrix(unlist(strsplit(cur$anot[[1]], "_")),
                                ncol = 2, byrow = TRUE))
        dplyr::bind_cols(cur[rep(1, nrow(df)), ], df)
    })
    GL <- do.call(rbind, GL)
    colnames(GL)[c(8, 9)] <- c("f", "ad")
    GL$f <- as.character(GL$f)
    GL$ad <- as.character(GL$ad)

    ##Adduct Priorization
    if (priorization == "only") {
        message(paste("Selecting only", adduct))
        GL <- GL[which(GL$ad %in% unlist(adduct)), ]
    }

    #Group entries by similar characterization attributes
    GL <- GLgroup(GL, rtmargin, ppm)

    #Tidy the entries, sort them by relevance and store them into a internal DF
    #Keep the important info for MS2 acquisition (rti, rtf and mass for
    #planning injections)
    GL <- GLtidy(GL, filterrt, filtermz)

    GL$entrynames <- vapply(GL$jointentries, function(entry) {
        res <- lapply(entry$metadata, function(subentry) {
            return(paste(unique(paste(subentry$f, subentry$ad, sep = "$"),
                            collapse = "#")))
        })
        return(paste(unique(unlist(res)), sep = "$", collapse = "#"))
    }, character(1))
    GL$entrynames <- unlist(GL$entrynames)

    PL_id <- which(struct@metadata@filenames == SoiList@filename)[1]
    raw <- PL(struct, PL_id)@raw
    if(nrow(raw) != 0){
        GL$XIC <- lapply(seq_len(nrow(GL)), function(x){
            cur <- GL[x,]
            # mzs <- c(min(cur$mass) * (1 - ppm * 1e-6),
            # max(cur$mass) * (1 + ppm * 1e-6))
            # mzs <- c(min(cur$mass) - filtermz,
                     # max(cur$mass) + filtermz)
            mzs <- c(min(cur$jointentries[[1]]$mass) * (1 - ppm * 1e-6),
                     max(cur$jointentries[[1]]$mass) * (1 + ppm * 1e-6))
            calculate_XIC_estimation(raw,
                                     mzs,
                                     c(cur$start, cur$end))
        })
    }
    GL <- RHermesIL(IL = as.data.table(GL[, -5]), annotation = GL$jointentries,
                    SOInum = id, ILParam = params)
    return(GL)
}

GLprior <- function(GL, ad, rtmargin, ppm) {
    if (!"adrows" %in% colnames(GL)) {
        warning("No adduct similarity annotation has been found on this SOI list
                Adduct priorization failed.
                Please use filterSOI() before priorization.")
        return(GL)
    }
    for (i in ad) {
        message(paste("Prioritizing", i))
        priorityentries <- GL[unlist(vapply(GL$ad, function(x) {
            any(x %in% i)
        }, FUN.VALUE = logical(1))), ]
        if (nrow(priorityentries) == 0) {
            next
        }
        connections <- unlist(lapply(priorityentries$adrows, function(ads) {
            ads[[i]]
        }))
        toremove <- which(GL$originalID %in% connections)
        for (j in toremove) {
            st <- GL$start[j]
            end <- GL$end[j]
            m <- GL$mass[j]
            idx <- which(between(GL$start, st - rtmargin, end + rtmargin) &
                            between(GL$end, st - rtmargin, end + rtmargin) &
                            between(GL$mass, m - m * 2*ppm/1e+06, m + m * 2 *
                                    ppm/1e+06))
            toremove <- c(toremove, idx)
        }
        if (length(toremove) == 0) {next}
        GL <- GL[-unique(toremove), ]
    }
    return(GL)
}

GLgroup <- function(GL, rtmargin, ppm) {
    message("Grouping the entries")
    out <- list()
    repeat {
        cur <- GL[1, ]
        st <- cur$start
        end <- cur$end
        m <- cur$mass
        idx <- which(between(GL$start, st - rtmargin, end + rtmargin) &
                    between(GL$end, st - rtmargin, end + rtmargin) &
                    between(GL$mass, m - m * ppm/1e+06, m + m * ppm/1e+06))
        out <- c(out, list(GL[idx, ]))
        GL <- GL[-idx, ]
        if (nrow(GL) == 0) {
            break
        }
    }
    return(out)
}

GLtidy <- function(GL, filterrt, filtermz) {
    message("Tidying the entries")
    GL <- lapply(GL, function(entry) {
        # idx <- order(-entry$adfound, -entry$isofound)
        idx <- order(-entry$MaxInt)
        entry <- entry[idx, ]
        chosen <- entry[1, c("start", "end", "mass", "MaxInt")]
        chosen$metadata <- list(list(entry))
        return(chosen)
    })
    GL <- do.call(rbind, GL)
    GL <- GL[order(GL$mass), ]
    out <- list()
    repeat {
        cur <- GL[1, ]
        st <- cur$start
        end <- cur$end
        m <- cur$mass
        idx <- which((between(GL$start, st - filterrt, end + filterrt) |
                        between(GL$end, st - filterrt, end + filterrt)) &
                        between(GL$mass, m - filtermz, m + filtermz))
        uniondf <- data.frame(start = min(GL[idx, "start"]),
                                end = max(GL[idx, "end"]),
                                mass = mean(unlist(GL[idx, "mass"])),
                                MaxInt = max(GL[idx, "MaxInt"]))
        uniondf$jointentries <- list(GL[idx, ])
        out <- c(out, list(uniondf))
        GL <- GL[-idx, ]
        if (nrow(GL) == 0) {
            break
        }
    }
    GL <- do.call(rbind, out)
    return(GL)
}
