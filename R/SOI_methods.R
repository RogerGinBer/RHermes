#'@title removeSOI
#'@description Removes SOI lists from an RHermesExp object
#'@slot struct The RHermesExp object
#'@slot idx The SOI lists to remove (indexes)
#'@return An RHermesExp object without the SOI lists indicated by idx.
#'@examples
#'if(FALSE){
#' myHermes <- removeSOI(myHermes, 2) #Remove the 2nd SOI list in the object
#'}

#'@export
setGeneric("removeSOI", function(struct, idx) {
    standardGeneric("removeSOI")
})
setMethod("removeSOI", c("RHermesExp", "numeric"), function(struct, idx) {
    validObject(struct)
    nSOI <- length(struct@data@SOI)
    if (nSOI == 0) {
        warning("No SOI to substract")
        return(struct)
    }
    tosub <- c()
    for (i in idx) {
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
setMethod("getSOIpar", c("ANY"), function(tag = "double") {
    temp <- read.csv2(system.file("extdata", "SOITemplates.csv",
                                    package = "RHermes"))
    specdf <- filter(temp, name == tag)[, seq(2, 6)]
    if (nrow(specdf) == 0) {
        stop("No templated was found with that ID", call. = FALSE, )
    }
    return(SOIParam(specs = specdf[, seq(1, 3)], maxlen = specdf[1, 4],
        minint = specdf[1, 5]))
})



#'@export
setGeneric("SOIcos", function(struct, id) {
    standardGeneric("SOIcos")
})
setMethod("SOIcos", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    if (length(struct@data@SOI) == 0) {
        stop("This object doesn't have any SOI lists")
    }
    if (!between(length(struct@data@SOI), 1, id)) {
        stop("Please enter a valid SOI list number")
    }
    SOI <- struct@data@SOI[[id]]@SoiList
    SOI <- SOI[order(SOI$start), ]
    m <- matrix(0, nrow = nrow(SOI), ncol = nrow(SOI))
    for (i in seq_len(nrow(SOI))) {
        st <- SOI$start[i]
        end <- SOI$end[i]
        for (j in seq_len(i)) {
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
            m[i, j] <- RHermes:::pearsonSim(SOI$peaks[[i]], SOI$peaks[[j]])
            m[j, i] <- m[i, j]
        }
    }
    return(m)
})


#'@export
setMethod("show", "SOIParam", function(object) {
    message("SOI parameters info:")
    message("\tBins used:")
    print(object@specs)
    message(paste("\tMax SOI length: ", object@maxlen))
    message(paste("\tMin data point intensity: ", object@minint))
    message(paste("\tBlank substraction performed: ", object@blanksub))
    message(paste("\tBlank filename: ", object@blankname))
})

#'@export
setMethod("show", "RHermesSOI", function(object) {
    message("Info about this SOI list:")
    message(paste("\tOriginal file name:", object@filename))
    message(paste("\tNumber of SOIs:", nrow(object@SoiList)))
    show(object@SoiParam)
})
