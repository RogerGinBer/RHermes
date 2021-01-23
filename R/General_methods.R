
#'@export
setMethod("show", "RHermesExp", function(object) {
    show(object@metadata@ExpParam)
    nfiles <- length(object@metadata@filenames)
    if(nfiles != 0){
        message(paste0("Number of processed MS1 files: ", nfiles ))
        lapply(object@metadata@filenames, function(x){message(paste0("\t", x))})
        show(object@data@PL)
        nsoi <- length(object@data@SOI)
        if(nsoi != 0){
            show(object@data@SOI)
            nMS2 <- length(object@data@MS2Exp)
            if(nMS2 != 0){
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
#' keep track of all changes made to the object (added files, generated SOI lists,
#' changed parameters, etc.)
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
