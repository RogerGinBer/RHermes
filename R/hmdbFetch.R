hmdbFetch <- function(form) {
    return(hmdb$Name[hmdb$MolecularFormula == toupper(form)])
}
