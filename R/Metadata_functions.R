adductTables <- function(ch_max = 1, mult_max = 1) {
    data(adducts, package = "enviPat", envir = environment())
    adducts$Mass[49] <- adducts$Mass[49] * (-1)  #Fixed wrong one

    negative.envi <- adducts[which(adducts$Ion_mode == "negative"), ]
    positive.envi <- adducts[which(adducts$Ion_mode == "positive"), ]

    ad_positive <- which(
        positive.envi$Charge %in% as.character(seq(ch_max)) &
        positive.envi$Mult %in% seq(mult_max)
    )
    ad_negative <- which(
        negative.envi$Charge %in% as.character(seq(from = -1, to = -ch_max)) &
        negative.envi$Mult %in% seq(mult_max)
    )

    positive.envi <- positive.envi[ad_positive, ]
    negative.envi <- negative.envi[ad_negative, ]

    negative.ad <- negative.envi[, -c(2, 9)]
    colnames(negative.ad)[c(1, 4)] <- c("adduct", "massdiff")
    negative.ad$massdiff <- negative.ad$massdiff * abs(negative.ad$Charge)

    positive.ad <- positive.envi[, -c(2, 9)]
    colnames(positive.ad)[c(1, 4)] <- c("adduct", "massdiff")
    positive.ad$massdiff <- positive.ad$massdiff * abs(positive.ad$Charge)

    return(list(negative.ad, positive.ad))
}

#' @importFrom utils data read.csv read.csv2 write.csv
database_importer <- function(template = "hmdb",
                                filename = "./app/www/norman.xls",
                                minmass = 70, maxmass = 750, keggpath = "") {
    if (template == "hmdb") {
        db <- read.csv(system.file("extdata", "hmdb.csv", package = "RHermes"))
        db <- db[!grepl("\\.", db$MolecularFormula), ]
    } else if (template == "norman") {
        db <- readxl::read_excel(system.file("extdata", "norman.xlsx",
                                                package = "RHermes"))
        names(db)[names(db) == "MOLECULAR_FORMULA"] <- "MolecularFormula"
        names(db)[names(db) == "NAME"] <- "Name"
    } else if (template == "custom") {
        if (grepl(pattern = "csv", x = filename)) {
            db <- read.csv(filename)
            if(ncol(db) == 1){
                db <- read.csv2(filename)
            }
        } else {
            db <- readxl::read_excel(filename)
        }
        if (!all(c("MolecularFormula", "Name") %in% names(db))) {
            stop("The colnames aren't adequate. The file must contain the
            'MolecularFormula' and 'Name' columns at least")
        }
    } else if (template == "kegg_p") {
        suppressWarnings(
            splitpath <- split(seq_along(keggpath),
                                seq_len(ceiling(length(keggpath)/10)))
        )
        comp <- lapply(splitpath, function(x) {
            pathdata <- KEGGREST::keggGet(keggpath[x])
            lapply(pathdata, function(y) {
                names(y$COMPOUND)
            })
        })
        comp <- unique(unlist(comp))
        suppressWarnings(splitcomp <- split(seq_along(comp),
                                            seq_len(ceiling(length(comp)/10))))
        db <- lapply(splitcomp, function(x) {
            data <- KEGGREST::keggGet(comp[x])
            data <- lapply(data, function(y) {
                c(y[1], gsub(y[2][[1]][[1]], pattern = ";", replacement = ""),
                    y[3])
            })
            as.data.frame(do.call(rbind, data))
        })
        db <- as.data.frame(do.call(rbind, db))
        names(db)[c(2, 3)] <- c("Name", "MolecularFormula")
    } else {
        stop("You haven't entered a valid template. Options are: 'hmdb',
        'norman' and 'custom'")
    }

    #Clean molecular formulas that contain unknown elements
    db <- db[grepl("^C.?.?", db$MolecularFormula), ]
    db <- db[!grepl("[", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("R", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("X", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl("T", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(".", db$MolecularFormula, fixed = TRUE), ]
    db <- db[!grepl(")", db$MolecularFormula, fixed = TRUE), ]
    db$MolecularFormula <- as.character(db$MolecularFormula)

    #Calculate M0 mass from formula
    isotopes <- NULL #To appease R CMD Check "no visible binding"
    data(isotopes, package = "enviPat", envir = environment())
    db$EnviPatMass <- lapply(db$MolecularFormula, function(x) {
        res <- tryCatch(enviPat::isopattern(isotopes, x, threshold = 99,
                                        verbose = FALSE),
                        error = function(cond){NA})
        # Isotopes should be loaded first
        if (is.matrix(res[[1]])) {
            res <- res[[1]][1, 1]
            res <- as.numeric(unname(res))
            return(res)
        } else {
            return(NA)
        }
    })
    if (any(is.na(db$EnviPatMass))) {
        db <- db[!is.na(db$EnviPatMass), ]
    }
    db$EnviPatMass <- as.numeric(db$EnviPatMass)
    db <- filter(db, dplyr::between(as.numeric(db$EnviPatMass),
                                        minmass, maxmass))
    return(db)
}
