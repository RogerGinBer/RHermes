adductTables <- function(ch_max = 1, mult_max = 1) {
    ads <- adducts_MetaboCore_to_Hermes()
    ads <- filter(ads, abs(Charge) <= ch_max)
    ads <- filter(ads, abs(Mult) <= mult_max)
    return(list(ads[ads$Ion_mode == "negative", ],
                ads[ads$Ion_mode == "positive", ]))
}

#' @importFrom MetaboCoreUtils adducts
adducts_MetaboCore_to_Hermes <- function(){
    ads <- rbind(MetaboCoreUtils::adducts("positive"),
                 MetaboCoreUtils::adducts("negative"))
    
    #Reorder columns
    ads <- ads[,c("name", "charge", "mass_multi", "mass_add", "positive",
                  "formula_add", "formula_sub")]
    
    #Adapt names and format
    names(ads) <- c("adduct", "Charge", "Mult", "massdiff", "Ion_mode",
                    "Formula_add", "Formula_ded")
    ads$Ion_mode <- ifelse(ads$Ion_mode, "positive", "negative")
    ads$Mult <- ads$Mult * abs(ads$Charge)
    ads$massdiff <- ads$massdiff * abs(ads$Charge)
    
    #Adapt atoms to add and subtract
    incomplete <- grepl("[A-Za-z]$", ads$Formula_add)
    ads$Formula_add[incomplete] <- paste0(ads$Formula_add[incomplete], "1")
    ads$Formula_add[ads$Formula_add == ""] <- "C0"
    
    
    incomplete <- grepl("[A-Za-z]$", ads$Formula_ded)
    ads$Formula_ded[incomplete] <- paste0(ads$Formula_ded[incomplete], "1")
    ads$Formula_ded[ads$Formula_ded == ""] <- "C0"
    
    ads
}

#' @importFrom readxl read_excel
#' @importFrom KEGGREST keggGet
#' @importFrom utils data read.csv read.csv2 write.csv
database_importer <- function(template = "custom",
                                filename,
                                minmass = 50, maxmass = 1200, keggpath = "") {
    if (template == "hmdb") {
        db <- read.csv(system.file("extdata", "hmdb.csv", package = "RHermes"))
        db <- db[!grepl("\\.", db$MolecularFormula), ]
    } else if (template == "norman") {
        db <- read_excel(system.file("extdata", "norman.xlsx",
                                                package = "RHermes"))
        names(db)[names(db) == "MOLECULAR_FORMULA"] <- "MolecularFormula"
        names(db)[names(db) == "NAME"] <- "Name"
    } else if (template == "custom") {
        if (grepl(pattern = "csv", x = filename)) {
            db <- read.csv(filename)
            if (ncol(db) == 1){
                db <- read.csv2(filename)
            }
        } else {
            db <- read_excel(filename)
        }
        if (!all(c("MolecularFormula", "Name") %in% names(db))) {
            stop("The colnames are not correct. The file must contain the
            'MolecularFormula' and 'Name' columns at least")
        }
    } else if (template == "kegg_p") {
        suppressWarnings({
            splitpath <- split(seq_along(keggpath),
                                seq_len(ceiling(length(keggpath)/10)))
        })
        comp <- lapply(splitpath, function(x) {
            pathdata <- keggGet(keggpath[x])
            lapply(pathdata, function(y) {
                names(y$COMPOUND)
            })
        })
        comp <- unique(unlist(comp))
        suppressWarnings({splitcomp <- split(seq_along(comp),
                                            seq_len(ceiling(length(comp)/10)))})
        db <- lapply(splitcomp, function(x) {
            data <- keggGet(comp[x])
            data <- lapply(data, function(y) {
                c(y[[1]],
                  gsub(y[2][[1]][[1]], pattern = ";", replacement = "")[[1]],
                  y[[3]])
            })
            as.data.frame(do.call(rbind, data))
        })
        db <- as.data.frame(do.call(rbind, db))
        names(db)[c(2, 3)] <- c("Name", "MolecularFormula")
    } else {
        stop("You haven't entered a valid template. Options are: 'hmdb',
        'norman' and 'custom'")
    }

    #Clean molecule names
    db$Name <- iconv(db$Name, from = "UTF-8", to = "ASCII", sub = "?")
    db$Name <- gsub('[^[:graph:][:space:]]', '', db$Name)
    
    #Clean molecular formulas that contain unknown elements
    db <- db[grepl("^C.?.?", db$MolecularFormula), ]
    db <- db[!grepl("[", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("R", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("X", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("T", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(".", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(")", db$MolecularFormula, fixed = TRUE), ]
    db$MolecularFormula <- as.character(db$MolecularFormula)
    
    # Clean deuterium in the formulae and substitute by [2H]
    db$MolecularFormula <- gsub(pattern = "D([0-9]*)",
                                replacement = "\\[2H\\1\\]",
                                x = db$MolecularFormula)
    
    #Calculate exact mass from formula
    db$ExactMass <- MetaboCoreUtils::calculateMass(db$MolecularFormula)
    if (any(is.na(db$ExactMass))) {
        invalid <- which(is.na(db$ExactMass))
        warning("Careful, some of the molecular formulas were not valid, first 10 are:",
                db$MolecularFormula[invalid[seq(1, min(10, length(invalid)))]])
        db <- db[!invalid, ]
    }
    db$ExactMass <- as.numeric(db$ExactMass)
    db <- dplyr::filter(db, dplyr::between(as.numeric(db$ExactMass),
                                    minmass, maxmass))
    return(db)
}
