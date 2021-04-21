#'@title removeSOI
#'@author Roger Gine
#'@description Removes SOI lists from an RHermesExp object
#'@param struct The RHermesExp object
#'@param id The SOI lists to remove (indexes)
#'@return An RHermesExp object without the SOI lists indicated by idx.
#'@examples
#'if(FALSE){
#' myHermes <- removeSOI(myHermes, 2) #Remove the 2nd SOI list in the object
#'}
#'@export
setGeneric("removeSOI", function(struct, id) {
    standardGeneric("removeSOI")
})
#' @rdname removeSOI
setMethod("removeSOI", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    nSOI <- length(struct@data@SOI)
    if (nSOI == 0) {
        warning("No SOI to substract")
        return(struct)
    }
    tosub <- c()
    for (i in id) {
        if (i > nSOI) {
            warning(paste("The index", i,
                            "is larger than the number of SOI lists"))
        } else {
            tosub <- c(tosub, i)
        }
    }
    if (length(tosub) != 0) {
        struct@data@SOI <- struct@data@SOI[-tosub]
        struct <- setTime(struct, paste("Removed the SOI lists:", tosub))
    } else {
        warning("Couldn't remove any SOIs")
    }
    return(struct)
})

duplicateSOI <- function(struct, id){
    validObject(struct)
    if(!id %in% seq_along(struct@data@SOI))
    struct@data@SOI <- c(struct@data@SOI, SOI(struct, id))
    return(struct)
}


#' @title getSOIpar
#' @description Returns a SOIParam object for posterior SOI detection.
#' This function does NOT require the RHermesExp object as multple SOIParam
#' can be used at once. See \link[RHermes]{findSOI} for more info.
#' @param tag A character string that tells which premade SOI parameter
#' object to use. Currently the following tags are available: 'single',
#' 'double', 'triple', and their extended counterparts for longer
#' chromatography experiments, 'single-x', 'double-x' and 'triple-x'.
#' These are all stored in /app/www/SOIFilterParams.csv, feel free to
#' locally change them or add new ones for your use (if you know what
#' you're doing).
#' @return A SoiParam object
#' @examples
#' if(FALSE){
#' par <- getSOIpar('double')
#' par2 <- getSOIpar('triple-x') #Etc. etc.
#' }
#'@export
setGeneric("getSOIpar", function(tag = "double") {
    standardGeneric("getSOIpar")
})
#' @rdname getSOIpar
setMethod("getSOIpar", c("ANY"), function(tag = "double") {
    temp <- read.csv2(system.file("extdata", "SOITemplates.csv",
                                    package = "RHermes"))
    specdf <- filter(temp, .data$name == tag)[, seq(2, 6)]
    if (nrow(specdf) == 0) {
        stop("No templated was found with that ID", call. = FALSE, )
    }
    return(SOIParam(specs = specdf[, seq(1, 3)], maxlen = specdf[1, 4],
        minint = specdf[1, 5]))
})


#'@title SOIsim
#'@author Roger Gine
#'@description Calculates the elution profile similarity between SOIs.
#'@param struct The RHermesExp object.
#'@param id The SOI list to be used.
#'@param subset A subset of SOI list entries that you want to compare the
#'  profile similarities with. Defaults to NA, which means all the entries.
#'@param mode Whether to compare the similarity with "all" other SOIs or only
#'  with those determined by subset. Defaults to "all", any other value
#'  restricts the comparison to the subset.
#'@return A similarity matrix.
#'@examples
#'if(FALSE){
#' myHermes <- SOIsim(myHermes, 2) #Remove the 2nd SOI list in the object
#'}
#'@export
setGeneric("SOIsim", function(struct, id, subset = NA, mode = "all") {
    standardGeneric("SOIsim")
})

#' @rdname SOIsim
setMethod("SOIsim", c("RHermesExp", "numeric", "ANY", "ANY"),
        function(struct, id, subset = NA, mode = "all") {
    validObject(struct)
    if (length(struct@data@SOI) == 0) {
        stop("This object doesn't have any SOI lists")
    }
    if (!between(length(struct@data@SOI), 1, id)) {
        stop("Please enter a valid SOI list number")
    }
    if(is.na(subset)) subset <- seq_len(nrow(SOI))

    SOI <- struct@data@SOI[[id]]@SOIList
    if(mode != "all"){
        against <- subset
        SOI <- SOI[subset,]
    }
    SOI <- SOI[order(SOI$start), ]
    m <- matrix(0, nrow = nrow(SOI), ncol = nrow(SOI))
    for (i in subset) {
        st <- SOI$start[i]
        end <- SOI$end
        if(mode == "all") against <- seq_len(i)
        for (j in against) {
            if (i == j) {
                m[i, j] <- 1
                next
            }
            if (end < SOI$start[j]) {
                break
            }
            if (st > SOI$end[j]) {
                next
            }
            # m[i, j] <- cosineSim(SOI$peaks[[i]], SOI$peaks[[j]])
            m[i, j] <- pearsonSim(SOI$peaks[[i]], SOI$peaks[[j]])
            m[j, i] <- m[i, j]
        }
    }
    return(m)
})


#' @export
#' @rdname SOIParam-class
#' @param object A SOIParam object
setMethod("show", "SOIParam", function(object) {
    message("SOI parameters info:")
    message("\tBins used:")
    print(object@specs)
    message(paste("\tMax SOI length: ", object@maxlen))
    message(paste("\tMin data point intensity: ", object@minint))
    message(paste("\tBlank substraction performed: ", object@blanksub))
    message(paste("\tBlank filename: ", object@blankname))
})

#' @export
#' @rdname RHermesSOI-class
#' @param object A RHermesSOI object
setMethod("show", "RHermesSOI", function(object) {
    message("Info about this SOI list:")
    message(paste("\tOriginal file name:", object@filename))
    message(paste("\tNumber of SOIs:", nrow(object@SOIList)))
    show(object@SOIParam)
})
