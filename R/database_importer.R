database_importer <- function(template = "hmdb", filename = "./app/www/norman.xls",
    minmass = 70, maxmass = 750, keggpath = "") {
    if (template == "hmdb") {
        db <- read.csv(system.file("extdata", "hmdb.csv", package = "RHermes"))
        db <- db[!grepl("\\.", db$MolecularFormula), ]
    } else if (template == "norman") {
        db <- readxl::read_excel(system.file("extdata", "norman.xls",
            package = "RHermes"))
        names(db)[names(db) == "MOLECULAR_FORMULA"] <- "MolecularFormula"
        names(db)[names(db) == "NAME"] <- "Name"
    } else if (template == "custom") {
        if (grepl(pattern = "csv", x = filename)) {
            db <- read.csv(filename)
        } else {
            db <- readxl::read_excel(filename)
        }
        if (!all(c("MolecularFormula", "Name") %in% names(db))) {
            stop("The colnames aren't adequate. The file must contain the
           'MolecularFormula' and 'Name' columns at least")
        }
    } else if (template == "kegg_p"){
        suppressWarnings(
            splitpath <- split(seq_along(keggpath),
                               seq_len(ceiling(length(keggpath)/10)))
        )
        comp <- lapply(splitpath, function(x){
                                 pathdata <- keggGet(keggpath[x])
                                 lapply(pathdata, function(y){
                                     names(y$COMPOUND)
                                 })
                             })
        comp <- unique(unlist(comp))
        db <- lapply(split(seq_along(comp),
                            seq_len(ceiling(length(comp)/10))), function(x){
                                data <- keggGet(comp[x])
                                data <- lapply(data, function(y){
                                    c(y[1],
                                      gsub(y[2][[1]][[1]], pattern = ";",
                                           replacement = ""),
                                      y[3])
                                })
                                as.data.frame(do.call(rbind,data))
                            })
        db <- as.data.frame(do.call(rbind, db))
        names(db)[c(2,3)] <- c("Name", "MolecularFormula")
    } else {
        stop("You haven't entered a valid template. Options are: 'hmdb',
        'norman' and 'custom'")
    }

    db <- db[grepl("^C.?.?", db$MolecularFormula), ]
    db <- db[!grepl("[", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("R", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("X", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("T", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(".", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(")", db$MolecularFormula, fixed = TRUE), ]
    db$MolecularFormula <- as.character(db$MolecularFormula)

    data(isotopes, package = "enviPat", envir = environment())  # for formula search without mass
    db$EnviPatMass <- lapply(db$MolecularFormula, function(x) {
        res <- try(enviPat::isopattern(isotopes, x, threshold = 99,
            verbose = FALSE))
        # Isotopes should be loaded first
        if (class(res[[1]]) != "try-error") {
            res <- res[[1]][1, 1]
            res <- as.numeric(unname(res))
            return(res)
        } else {
            return(NA)
        }
    })
    if (any(is.na(db$EnviPatMass))) {
        idx <- is.na(db$EnviPatMass)
        db <- db[-idx, ]
    }
    db$EnviPatMass <- as.numeric(db$EnviPatMass)

    suppressWarnings({
        db <- db[which(dplyr::between(as.numeric(db$EnviPatMass),
                                      minmass, maxmass)), ]
    })
    return(db)
}



