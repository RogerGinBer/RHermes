#'@export
setMethod("show", "RHermesExp", function(object) {
    show(object@metadata@ExpParam)
    nfiles <- length(object@metadata@filenames)
    if (nfiles != 0) {
        message(paste0("Number of processed MS1 files: ", nfiles ))
        lapply(object@metadata@filenames, function(x){message(paste0("\t", x))})
        show(object@data@PL)
        nsoi <- length(object@data@SOI)
        if (nsoi != 0) {
            show(object@data@SOI)
            nMS2 <- length(object@data@MS2Exp)
            if (nMS2 != 0) {
                show(object@data@MS2Exp)
            }
        } else {
            message("No SOI lists generated yet.")
        }
    } else {
        message("No files processed yet.")
    }
    message("Processing Timestamps ---------------------------------------")
    readTime(object)
})


setGeneric("setTime", function(struct, message) {
    standardGeneric("setTime")
})
setMethod("setTime", c("RHermesExp", "character"), function(struct,
    message) {
    validObject(struct)
    struct@metadata@timestamps <- c(struct@metadata@timestamps,
        paste(date(), message, sep = " - "))
    return(struct)
})



#' @title readTime
#' @description Prints all timestamps of a given RHermesExp object. Useful to
#' keep track of all changes made to the object (added files, generated SOI
#'  lists, changed parameters, etc.)
#'
#' @param struct The RHermesExp object you want to read the timestamps from.
#'
#' @examples
#' if(FALSE){
#' readTime(myHermes)
#' }
#'
#' @return None, just prints the timestamps
#'@export
setGeneric("readTime", function(struct) {
    standardGeneric("readTime")
})
setMethod("readTime", "RHermesExp", function(struct) {
    validObject(struct)
    lapply(struct@metadata@timestamps, function(m) {
        message(m)
    })
    message("")
})

## Object getters ##
#'@export
setGeneric("DB", function(struct) {
    standardGeneric("DB")
})
setMethod("DB", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@ExpParam@DB
})

#'@export
setGeneric("adlist", function(struct) {
    standardGeneric("adlist")
})
setMethod("adlist", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@ExpParam@adlist
})

#'@export
setGeneric("PL", function(struct, id) {
    standardGeneric("PL")
})
setMethod("PL", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@PL[[id]]
})


#'@export
setGeneric("SOI", function(struct, id) {
    standardGeneric("SOI")
})
setMethod("SOI", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@SOI[[id]]
})


#'@export
setGeneric("IL", function(struct, id) {
    standardGeneric("IL")
})
setMethod("IL", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@MS2Exp[[id]]@IL
})


#'@export
setGeneric("MS2Data", function(struct, id) {
    standardGeneric("MS2Data")
})
setMethod("MS2Data", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@MS2Exp[[id]]@MS2Data
})

#'@export
setGeneric("Ident", function(struct, id) {
    standardGeneric("Ident")
})
setMethod("Ident", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@MS2Exp[[id]]@Ident[[1]]
})

#'@export
setGeneric("Cluster", function(struct) {
    standardGeneric("Cluster")
})
setMethod("Cluster", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@cluster
})


