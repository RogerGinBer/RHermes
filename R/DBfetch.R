#Allows to extract formulas or formula-adducts that present overlap with others from the Overlap_xx scripts
DBfetch <- function(DB_results, form, adduct = "") {
    if (adduct != "") {
        query <- paste(form, adduct, sep = "_")
        return(DB_results[as.character(DB_results[, 1]) == query, 
            ])
    } else {
        query <- paste(form, "_", sep = "")
        return(DB_results[grepl(x = DB_results[, 1], pattern = query), 
            ])
    }
}

