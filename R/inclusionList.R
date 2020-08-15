#'@title genIL
#'
#'@description Groups the SOIs generated in the previous step and sets them up in an inclusion list ready
#'for the MS/MS analysis. It takes into account the instrument precursor filter mz range and groups the
#'SOIs accordingly. The user can also specify whether they want to use all SOIs or prioritize SOIs with a
#'certain annotation (M+H adducts, prioritize SOIs with adequate isotopic patterns, etc.).
#'
#'
#'@param struct The RHermesExp object in which to save the resulting IL
#'@param id The SOI list ID from which to generate the IL
#'@param par An inclusion list parameter object, generated with ILParam()
#'
#'@return An RHermes object with the IL inside a new MS2Exp entry, inside the data slot
#'
#'@seealso \link[RHermes]{ILParam} \link[RHermes]{RHermesMS2Exp}
#'
#'@examples
#'if(FALSE){
#' par <- ILParam(filtermz = 0.5, priorization = "full")
#' myHermes <- genIL(myHermes, 1, par)
#'}
#' @importFrom dplyr between
#'@export
setGeneric("genIL", function(struct, id, par) {
    standardGeneric("genIL")
})
setMethod("genIL", c("RHermesExp", "numeric", "ANY"), function(struct,
    id, par = ILParam()) {
    validObject(struct)
    validObject(par)
    if (!is(par, "ILParam")) {
        stop("Please enter a valid IL parameter set with ILParam()")
    }
    if (id > length(struct@data@SOI)) {
        stop("The entered ID is greater than the number of SOI lists")
    }
    message("Now processing the IL:")
    struct@data@MS2Exp <- c(struct@data@MS2Exp, RHermesMS2Exp(IL = inclusionList(struct,
        par, id), MS2Data = list(), Ident = list()))
    message("Done!")
    return(struct)
})

#'@export
inclusionList <- function(struct, params, id) {
    SoiList <- struct@data@SOI[[id]]
    ppm <- struct@metadata@ExpParam@ppm
    GL <- SoiList@SoiList
    adlist <- struct@metadata@ExpParam@adlist

    rtmargin <- params@rtmargin
    priorization <- params@priorization
    adduct <- params@ad
    filtermz <- params@filtermz
    filterrt <- params@filterrt

    if (priorization == "yes") {
        message(paste("Prioritizing", adduct))
        GL <- RHermes:::GLprior(GL, adduct)
        low <- which(GL$MaxInt < 50000)
        rare <- which(vapply(GL$ad[low], function(x) {
            !any(unlist(x) %in% c("M+H", "M+NH4", "M+Na", "M+K",
                "M+CH3OH", "M-H", "M+Cl", "M+Br", "M+H2O-H"))
        }, logical(length(low))))
        if (length(rare) != 0) {
            GL <- GL[-low[rare], ]
        }
    }
    GL <- GL[, c("start", "end", "formula", "mass", "MaxInt",
        "anot")]
    GL$originalSOI <- seq_len(nrow(GL))

    GL <- lapply(seq_len(nrow(GL)), function(x) {
        cur <- GL[x, ]
        df <- as.data.frame(matrix(unlist(strsplit(cur$anot[[1]],
            "_")), ncol = 2, byrow = TRUE))
        dplyr::bind_cols(cur[rep(1, nrow(df)), ], df)
    })
    GL <- do.call(rbind, GL)
    colnames(GL)[c(8, 9)] <- c("f", "ad")
    GL$f <- as.character(GL$f)
    GL$ad <- as.character(GL$ad)

    ##1st Step: Remove confirmed isotope entries from the SoiList--------------------------
    # message('Removing isotopes:')
    # GL <- GLisorem(GL, rtmargin, ppm)

    ##Adduct Priorization ---------------------------------------------------------------
    if (priorization == "only") {
        message(paste("Selecting only", adduct))
        GL <- GL[which(GL$ad %in% unlist(adduct)), ]
    }


    #2nd Step: Group entries by similar characterization attributes------------------------
    GL <- RHermes:::GLgroup(GL, rtmargin, ppm)

    #3rd Step: Tidy the entries, sort them by relevance and store them into a internal DF
    #Keep the important info outside for MS/MS acquisition (rti, rtf and mass for
    #planning injections)------------------------------------------------------------------

    GL <- RHermes:::GLtidy(GL, filterrt, filtermz)

    # print(ggplot(GL) +
    #         geom_tile(aes(x = (start+end)/2, y = mass, width = end-(start+end)/2, height = 1)) +
    #         xlab('RT') + ylab('mz') +
    #         ggtitle('Representation of IL entries in the mz-rt space'))

    GL$entrynames <- vapply(GL$jointentries, function(entry) {
        res <- lapply(entry$metadata, function(subentry) {
            return(paste(subentry$f, subentry$ad, sep = "$",
                collapse = "#"))
        })
        return(paste(unlist(res), sep = "$", collapse = "#"))
    }, character(1))
    GL$entrynames <- unlist(GL$entrynames)
    GL <- RHermesIL(IL = as.data.table(GL[, -5]), annotation = GL$jointentries,
        SOInum = id, ILParam = params)
    return(GL)
}

GLisorem <- function(GL, rtmargin, ppm) {
    rowstoremove <- numeric()
    for (i in seq_len(nrow(GL))) {
        isomasses <- GL[i, ]$isodf[[1]][, 2]
        st <- GL[i, ]$start
        end <- GL[i, ]$end
        idx <- which(between(GL$start, st - rtmargin, st + rtmargin) &
            between(GL$end, end - rtmargin, end + rtmargin))  #Entries in window
        entrymass <- GL[idx, ]$mass
        ol <- lapply(isomasses, function(x) {
            thr <- c(x - ppm * 1e-06 * x, x + ppm * 1e-06 * x)
            if (any(entrymass > thr[1] & entrymass < thr[2])) {
                return(which(entrymass > thr[1] & entrymass <
                  thr[2]))
            }
            return()
        })
        ol <- do.call(rbind, ol)
        if (length(ol) != 0) {
            idx <- idx[ol]
            rowstoremove <- c(rowstoremove, idx)
        }
    }
    if (length(rowstoremove) != 0) {
        GL <- GL[-unique(rowstoremove)]
    }
    return(GL)
}


GLprior <- function(GL, ad) {

    originalrows <- seq_len(nrow(GL))

    for (i in ad) {
        priorityentries <- GL[GL$ad == i, ]
        if (nrow(priorityentries) == 0) {
            next
        }
        toremove <- which(originalrows %in% unique(unlist(priorityentries$adrows[[i]])))
        if (length(toremove) == 0) {
            next
        }
        originalrows <- originalrows[-toremove]
        GL <- GL[-toremove, ]
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
        idx <- which(between(GL$start, st - rtmargin, st + rtmargin) &
            between(GL$end, end - rtmargin, end + rtmargin) &
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
    GL <- GL[order(GL$start, GL$mass), ]
    out <- list()
    repeat {
        cur <- GL[1, ]
        st <- cur$start
        end <- cur$end
        m <- cur$mass
        idx <- which((between(GL$start, st - filterrt, st + filterrt) |
            between(GL$end, end - filterrt, end + filterrt)) &
            between(GL$mass, m - filtermz, m + filtermz))
        uniondf <- data.frame(start = min(GL[idx, "start"]),
            end = max(GL[idx, "end"]), mass = mean(unlist(GL[idx,
                "mass"])), MaxInt = max(GL[idx, "MaxInt"]))
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
