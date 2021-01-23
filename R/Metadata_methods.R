#'@export
setMethod("show", "ExpParam", function(object){
    message("Experimental parameters info:")
    message(paste0("\tppm error: ", object@ppm))
    message(paste0("\tMass resolution: ", object@res))
    message(paste0("\tMass range: ", object@minmz, "-",
                    object@maxmz))
    message(paste0("\tNoise cutoff value: ", object@nthr))
    message(paste0("\tIonization mode: ", object@ion))
    message(paste0("\tNumber of compounds: ", nrow(object@DB)))
    message(paste0("\tNumber of unique formulas: ",
                    length(unique(object@DB$MolecularFormula))))
    message(paste0("\tNumber of adducts: ", nrow(object@adlist)))
    if (nrow(object@adlist) != 0) {
        message(paste0("\tAdduct list: ", paste(object@adlist$adduct,
                                                collapse = ", ")))
    }
    message()
})

#' @title setExpParam
#' @description Main interface to set the experimental parameters into a
#'  RHermesExp object. The user can either provide an ExpParam object 
#' (partially or fully customized) or provide a template tag to use a pre-made
#'  parameter set.
#' @slot struct The RHermesExp object to modify.
#' @slot params Optional if using template. An ExpParam object detailing a
#'  custom set of parameters.
#' @slot template Ignored if using params. A character describing the ID of
#'  the pre-made parameter set to use. Currently you can use: 'orbi-pos',
#' 'orbi-neg', 'qtof-pos' and 'qtof-neg'. The parameters can be found,
#'  modified or even expanded in the csv at app/www/InstrumentTemplates.csv
#' @return An RHermesExp object with the desired experimental parameters set.
#' @examples
#' if(FALSE){
#' myHermes <- RHermesExp()
#' 
#' #Using template
#' myHermes <- setExpParam(myHermes, template = 'orbi-pos') 
#' 
#' #Custom parameters
#' myHermes <- setExpParam(myHermes,
#'                         ExpParam(ppm = 10, res = 135000, ion = '-')) 
#' }
#'
#'@export
setGeneric("setExpParam", function(struct, params = ExpParam(),
                                    template = character(0)) {
    standardGeneric("setExpParam")
})
setMethod("setExpParam", signature = c("RHermesExp", "ANY", "ANY"),
function(struct, params = ExpParam(), template = character(0)) {
    validObject(struct)
    if (length(template) == 0) {
        struct <- setTime(struct,
                        "Changed Experimental Parameters to custom template")
    } else {
        if (length(template) > 1) {
            warning("More than one template inserted, using only the first one")
        }
        tryCatch(expr = {
            csvfile <- read.csv2(system.file("extdata",
                                                "InstrumentTemplates.csv",
                                                package = "RHermes"),
                                    stringsAsFactors = F)
        })
        row <- which(csvfile$name == template[1])
        if (length(row) == 0) {
            stop("No template matches the inserted name ID")
        }
        params <- ExpParam(ppm = csvfile$ppm[row], res = csvfile$res[row],
                        nthr = csvfile$nthr[row], minmz = csvfile$minmz[row],
                        maxmz = csvfile$maxmz[row], ion = csvfile$ion[row],
                        instr = csvfile$instr[row])
        struct <- setTime(struct, paste0("Changed Experimental Parameters to ",
                                        template, " template"))
        }
    struct@metadata@ExpParam <- params
    return(struct)
})

#' @title setDB
#' @description It updates the RHermesExp object and sets the formula database
#'  and adduct list to use in the following steps of the workflow. The adduct
#'   list used is based on the
#'  \href{http://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00941}{EnviPat}
#' package
#' @slot struct The RHermesExp object to update
#' @slot db A character defining which database to use. Currently we only 
#' support 'hmdb' \href{https://hmdb.ca/}{HMDB} and 'norman'
#'  \href{https://www.norman-network.com/nds/susdat/}{NORMAN},
#'  but  new database support and customizability will come in future releases.
#' @slot adcharge The maximum charge (in absolute value, so polarity does not
#'  matter) the adducts in the list can have. Defaults to 1.
#' @slot admult The maximum multiplicity (M, 2M, 3M and so on) the adducts
#'  can have. Defaults to 1.
#' @slot folder Folder where the database is located.
#' @return An RHermesExp object with the formula and adduct database set.
#' @examples
#' if(FALSE){
#'     myHermes <- setDB(myHermes, 'hmdb') #Adcharge and admult default to 1
#'     myHermes <- setDB(myHermes, 'norman', 2, 1) #Charge 2, multiplicity 1
#' }
#'
#'
#'@export
setGeneric("setDB", function(struct, db = "hmdb", adcharge = 1,
                            admult = 1, filename = "./app/www", keggpath = "") {
    standardGeneric("setDB")
})
setMethod("setDB", signature = c("RHermesExp", "ANY", "ANY", "ANY", "ANY",
                                    "ANY"),
function(struct, db = "hmdb", adcharge = 1, admult = 1, filename = "./app/www",
            keggpath = "") {
    validObject(struct)
    message(paste("Parsing the", db, "formula database"))
    minmz <- struct@metadata@ExpParam@minmz
    maxmz <- struct@metadata@ExpParam@maxmz
    dbdata <- switch(db,
        hmdb = RHermes:::database_importer("hmdb",
                                            minmass = minmz,
                                            maxmass = maxmz),
        norman = RHermes:::database_importer("norman",
                                            minmass = minmz,
                                            maxmass = maxmz),
        custom = RHermes:::database_importer("custom",
                                            filename = filename,
                                            minmass = minmz,
                                            maxmass = maxmz),
        kegg_p = RHermes:::database_importer("kegg_p",
                                            keggpath = keggpath,
                                            minmass = minmz,
                                            maxmass = maxmz),
        stop("Please input a valid DB name")
    )
    struct@metadata@ExpParam@DB <- dbdata
    #Set corresponding adduct table
    ad <- RHermes:::adductTables(adcharge, admult)
    ion <- struct@metadata@ExpParam@ion
    if (ion == "+") {
        struct@metadata@ExpParam@adlist <- ad[[2]]
    } else {
        struct@metadata@ExpParam@adlist <- ad[[1]]
    }
    struct <- setTime(struct, paste("Added the", db,
                                    "formula database and an adduct list",
                                    "with charge", adcharge, "and multiplicity",
                                    admult))
    return(struct)
})

#' @title addAd
#' @description Adds new custom adducts to the RHermesExp adduct list
#' @details The function adds an entry to the Envipat-style dataframe containing
#' the adducts used for the experiment. It is important to specify what atoms to
#' add or substract to the molecular formula. Polarity is inferred directly from
#' the specified charge.
#' @return An updated RHermesExp object with the new adduct list
#'
#' @param struct The RHermesExp object
#' @param name Name of the adduct to add
#' @param ch Adduct charge. Remember to specify the symbol if negative.
#' @param mult Multiplicity (eg, M, 2M, 3M, etc.)
#' @param toadd Atoms to add to the formula when it ionizes as the adduct.
#' (ex: M+H -> H1, M+Cl -> Cl1). Note the explicit use of a 1 after the symbol.
#' @param tosub Atoms to substract to the formula. Same idea as in toadd.
#' (ex: M-H2O+H -> H2O1)
#'
#' @examples
#' if(FALSE){
#'  addAd(myHermes, 'M+H', 1.0072, ch = 1, mult = 1, toadd = 'H1')
#' }
#' @export
setGeneric("addAd", function(struct, name, deltam, ch = 1, mult = 1,
                                toadd = "FALSE", tosub = "FALSE") {
    standardGeneric("addAd")
})
setMethod("addAd", c("RHermesExp", "character", "numeric", "ANY", "ANY", "ANY",
                        "ANY"),
function(struct, name, deltam, ch = 1, mult = 1, toadd = "FALSE",
            tosub = "FALSE") {
    validObject(struct)
    if (ch == 0 | mult == 0) {
        warning("Invalid charge or multiplicity, they can't be zero")
        return(struct)
    }
    newentry <- data.frame(adduct = name, Charge = ch, Mult = mult,
                            massdiff = deltam,
                            Ion_mode = ifelse(ch > 0, "positive", "negative"),
                            Formula_add = toadd, Formula_ded = tosub)
    struct@metadata@ExpParam@adlist <- rbind(struct@metadata@ExpParam@adlist,
                                                newentry)
    message("This is the new adduct list: ")
    show(struct@metadata@ExpParam@adlist)
    struct <- setTime(struct, paste("Added the adduct", name,
                                    "with mass", deltam))
    return(struct)
})

#' @title remAd
#' @description Remove adducts from an RHermesExp adduct list
#' @details Removes a row (or many) from the adduct list slot inside the
#' RHermesExp object
#' @return An updated RHermesExp object
#'
#' @param struct The RHermesExp object
#' @param name The names of the adducts to remove. If they aren't in the
#' RHermesExp will raise a warning but no rows will be substracted.
#'
#' @examples
#' if(FALSE){
#'  remAd(myHermes, 'M+H') #Just one
#'  remAd(myHermes, c('M+DMSO+H', 'M+2K-H', 'M+CH3OH+H')) #Or many at a time
#' }
#' @export

setGeneric("remAd", function(struct, name) {
    standardGeneric("remAd")
})
setMethod("remAd", c("RHermesExp", "character"), function(struct, name) {
    validObject(struct)
    nad <- nrow(struct@metadata@ExpParam@adlist)
    if (nad == 0) {
        warning("No adducts to substract")
        return(struct)
    }
    tosub <- c()
    for (i in seq_along(name)) {
        cur <- which(struct@metadata@ExpParam@adlist$adduct ==
                        name[i])
        if (length(cur) == 0) {
            warning(paste("Adduct", name[i], "not found\n"))
        }
        tosub <- c(tosub, cur)
    }
    if (length(tosub) != 0) {
        struct@metadata@ExpParam@adlist <- 
            struct@metadata@ExpParam@adlist[-tosub,]
    }
    struct <- setTime(struct, paste("Removed the adduct/s", name))
    return(struct)
})
