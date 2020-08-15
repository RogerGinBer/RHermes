#Binds hmdb formulas by ppm proximity
#' @export
Binded_hmdb <- function(hmdb, ppm = 1) {
    umasses <- unique(hmdb$EnviPatMass)
    hmdbM <- hmdb$EnviPatMass
    
    hmdb3 <- lapply(umasses, function(m) {
        idx <- abs((hmdbM - m)) * 1e+06/m
        idx <- which(idx <= ppm)
        fms <- unique(hmdb$MolecularFormula[idx])
        fms <- paste(fms, collapse = ", ")
        data.frame(m, fms)
    })
    hmdb3 <- do.call("rbind", hmdb3)
    return(hmdb3)
    # hmdb3[grep(',',hmdb3$fms),]
}

#' @export
Binded_norman <- function(norman, ppm = 1) {
    umasses <- unique(norman$EnviPatMass)
    normanM <- as.numeric(norman$EnviPatMass)
    
    norman3 <- lapply(umasses, function(m) {
        idx <- abs((normanM - m)) * 1e+06/m
        idx <- which(idx <= ppm)
        fms <- unique(norman$MOLECULAR_FORMULA[idx])
        fms <- paste(fms, collapse = ", ")
        data.frame(m, fms)
    })
    norman3 <- do.call("rbind", norman3)
    return(norman3)
    # hmdb3[grep(',',hmdb3$fms),]
}
