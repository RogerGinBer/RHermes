
#'@export
#'
#'
#'
#Counts carbon atoms in a given formula or list of formulas. Useful for calculating possible 13C
#isotopologues in ScanSearch

CarbonCount <- function(FormulaList) {
    lapply(FormulaList, function(x) {
        x <- strsplit(x = x[[1]][1], split = ", ")[1]
        x <- x[[1]][1]
        # print(x)
        f <- CHNOSZ::makeup(x)
        if (any(names(f) == "C")) {
            return(f[which(names(f) == "C")])
        } else return(0)
    })
}
