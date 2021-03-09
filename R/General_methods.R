#'@export
#' @rdname RHermesExp-class
#' @param object An RHermesExp object
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

#' @title setTime
#' @author Roger Gine
#' @description Adds a timestamp to the RHermesExp object with a specified
#' message. Useful for adding info to the object for reproducibility purposes.
#' @param struct The RHermesExp object you want to add the timestamp to.
#' @param message A character string to append.
#' @examples setTime(RHermesExp(), "Important info")
#' @return A copy of struct with the updated timestamp list.
#' @export
setGeneric("setTime", function(struct, message) {
    standardGeneric("setTime")
})
#' @rdname setTime
setMethod("setTime", c("RHermesExp", "character"), function(struct,
    message) {
    validObject(struct)
    struct@metadata@timestamps <- c(struct@metadata@timestamps,
        paste(date(), message, sep = " - "))
    return(struct)
})



#' @title readTime
#' @author Roger Gine
#' @description Prints all timestamps of a given RHermesExp object. Useful to
#' keep track of all changes made to the object (added files, generated SOI
#'  lists, changed parameters, etc.)
#' @param struct The RHermesExp object you want to read the timestamps from.
#' @examples readTime(RHermesExp())
#' @return None, just prints the timestamps
#'@export
setGeneric("readTime", function(struct) {
    standardGeneric("readTime")
})
#' @rdname readTime
setMethod("readTime", "RHermesExp", function(struct) {
    validObject(struct)
    lapply(struct@metadata@timestamps, function(m) {
        message(m)
    })
    message("")
})

## Object getters ##
#'@title DB
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@return A data.frame with all molecular formulas in struct.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'DB(struct)
#'@export
setGeneric("DB", function(struct) {
    standardGeneric("DB")
})
#' @rdname DB
setMethod("DB", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@ExpParam@DB
})

#'@title adlist
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@return The current data.frame of adducts in struct.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'adlist(struct)
#'@export
setGeneric("adlist", function(struct) {
    standardGeneric("adlist")
})
#' @rdname adlist
setMethod("adlist", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@ExpParam@adlist
})

#'@title PL
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@param id The number corresponding to the PL you want to access
#'@return The RHermes_PL object at position id.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'PL(struct, 1)
#'@export
setGeneric("PL", function(struct, id) {
    standardGeneric("PL")
})
#' @rdname PL
setMethod("PL", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    validIds <- seq_along(struct@data@PL)
    if(id %in% validIds) return(struct@data@PL[[id]])
    warning(paste("You have input a wrong PL id, valid ids are:",
                min(validIds), "-", max(validIds)))
})

#'@title SOI
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@param id The number corresponding to the SOI list you want to access
#'@return The RHermes_SOI object at position id.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'SOI(struct, 1)
#'@export
setGeneric("SOI", function(struct, id) {
    standardGeneric("SOI")
})
#' @rdname SOI
setMethod("SOI", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    validIds <- seq_along(struct@data@SOI)
    if(id %in% validIds) return(struct@data@SOI[[id]])
    warning(paste("You have input a wrong SOI id, valid ids are:",
                min(validIds), "-", max(validIds)))
})

#'@title IL
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@param id The number corresponding to the inclusion list you want to access
#'@return The RHermes_IL object at position id.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'IL(struct, 1)
#'@export
setGeneric("IL", function(struct, id) {
    standardGeneric("IL")
})
#' @rdname IL
setMethod("IL", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    validIds <- seq_along(struct@data@MS2Exp)
    if(id %in% validIds) return(struct@data@MS2Exp[[id]]@IL)
    warning(paste("You have input a wrong IL id, valid ids are:",
                min(validIds), "-", max(validIds)))
})

#'@title MS2Data
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@param id The number corresponding to the inclusion list you want to access
#'@return A list of all IL entries of the RHermes_IL object at position id with
#'their corresponding MS2 scan data
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'MS2Data(struct, 1)
#'@export
setGeneric("MS2Data", function(struct, id) {
    standardGeneric("MS2Data")
})
#' @rdname MS2Data
setMethod("MS2Data", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@MS2Exp[[id]]@MS2Data
    if(length(struct@data@MS2Exp) > 0){
        validIds <- seq_along(struct@data@MS2Exp)
        if(id %in% validIds) return(struct@data@MS2Exp[[id]]@MS2Data)
        warning(paste("You have input a wrong MS2Data id, valid ids are:",
                min(validIds), "-", max(validIds)))
    } else {
        warning("No MS2Exp in the object")
    }
})

#'@title Ident
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@param id The number corresponding to the inclusion list you want to access
#'@return The identification data.frame with the RHermes_MS2Exp at position id.
#'@examples
#'\dontshow{struct <- readRDS(system.file("extdata", "exampleObject.rds",
#'                              package = "RHermes"))}
#'Ident(struct, 1)
#'@export
setGeneric("Ident", function(struct, id) {
    standardGeneric("Ident")
})
#' @rdname Ident
setMethod("Ident", c("RHermesExp", "numeric"), function(struct, id) {
    validObject(struct)
    struct@data@MS2Exp[[id]]@Ident[[1]]
})


#'@title Cluster
#'@author Roger Gine
#'@family Getters
#'@param struct An RHermesExp object
#'@return The BiocParallel backend object associated with struct
#'@examples Cluster(RHermesExp())
#'@export
setGeneric("Cluster", function(struct) {
    standardGeneric("Cluster")
})
#' @rdname Cluster
setMethod("Cluster", c("RHermesExp"), function(struct) {
    validObject(struct)
    struct@metadata@cluster
})


