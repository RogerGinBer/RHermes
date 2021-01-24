#' @title ExpParam
#' @description ExpParam holds all experimental data as well as database-related
#'  info generated during the processing steps. Use setExpParams to save it into
#'  your RHermesExp object
#' @slot ppm The mass spectrometer (MS) mass error in
#'  parts per million (ppm)
#' @slot res The MS resolution at m/z = 200
#' @slot nthr The noise threshold to consider. Any signal
#'  weaker than it will NOT be considered in the PL. Useful to filter
#'  'grass-like' peaks in a Q-TOF instrument. Defaults to 1000.
#' @slot minmz Minimum mz registered in the MS1 experiment.
#' Defaults to 80
#' @slot maxmz Maximum mz registered in the MS1 experiment.
#' Defaults to 1000.
#' @slot ion Polarity used in the experiment. RHermes currently
#'  only supports  one polarity at a time per RHermesExp, so if you're
#'  performing both + and - polarity you will have to use 2 different
#'  RHermesExp objects, one for each. Important: must be described ONLY as
#'  '+' or '-'.
#' @slot instr MS instrument type used to acquire the data. It can be
#' either QTOF' or 'Orbitrap'.
#' @slot DB Formula database
#' @slot adlist Adduct list
#' @slot ionF Ionic formula list. Formed by all formula-adduct
#'  combinations
#' @slot isoList Isotopic pattern exploration results.
#' @export ExpParam
#' @examples
#' if(FALSE){
#'  #Default object
#'  exp <- ExpParam()
#'  #Setting some values - check all possible slots
#'  exp2 <- ExpParam(ppm = 3, res = 150000, instr = "Orbitrap")
#' }
ExpParam <- setClass("ExpParam", slots = list(ppm = "numeric",
    res = "numeric", nthr = "numeric", minmz = "numeric", maxmz = "numeric",
    ion = "character", instr = "character", DB = "list", adlist = "data.frame",
    ionF = "list", isoList = "list"), prototype = list(ppm = 5,
    res = 60000, nthr = 100, minmz = 80, maxmz = 1000, ion = "+",
    instr = "QTOF", DB = list(), adlist = data.frame(), ionF = list(),
    isoList = list()))


#'@export SOIParam
SOIParam <- setClass("SOIParam", slots = list(specs = "data.frame",
    maxlen = "numeric", minint = "numeric", blanksub = "logical",
    blankname = "character"), prototype = list(specs = data.frame(),
    maxlen = 30, minint = 100, blanksub = FALSE, blankname = "None"))


#'@export RHermesPL
RHermesPL <- setClass("RHermesPL", slots = list(peaklist = "data.table",
    header = "data.frame", raw = "data.table", labelled = "logical",
    filename = "character"))


#'@export RHermesSOI
RHermesSOI <- setClass("RHermesSOI", slots = list(SoiList = "data.table",
    PlotDF = "data.table", SoiParam = "SOIParam", filename = "character"))


#' @export ILParam
ILParam <- setClass("ILParam", slots = list(filtermz = "numeric",
    filterrt = "numeric", rtmargin = "numeric", priorization = "character",
    ad = "character"), prototype = list(filtermz = 0.5, filterrt = 10,
    rtmargin = 5, priorization = "only", ad = c("M+H")))


#'@export RHermesIL
RHermesIL <- setClass("RHermesIL", slots = list(IL = "data.table",
    annotation = "list", SOInum = "numeric", ILParam = "ILParam"))


#' @export RHermesMS2Exp
RHermesMS2Exp <- setClass("RHermesMS2Exp", slots = list(IL = "RHermesIL",
    MS2Data = "list", Ident = "list"))


RHermesIdent <- setClass("RHermesIdent",
                            slots = list(IdentifiedSOI = "data.table",
    CompoundList = "list", MSMSMatchings = "list", Metadata = "list"))


#' @import BiocParallel
setRHermesCluster <- function(){
    if (Sys.info()[1] == "Windows") {
        ram <- system2("wmic", args =  "OS get FreePhysicalMemory /Value",
                        stdout = TRUE)
        ram <- ram[grepl("FreePhysicalMemory", ram)]
        ram <- gsub("FreePhysicalMemory=", "", ram, fixed = TRUE)
        ram <- gsub("\r", "", ram, fixed = TRUE)
        ram <- as.integer(ram)

        #Suppose max 2GB per worker
        nwork <- min(floor(ram/2e6), BiocParallel::snowWorkers())
        if (nwork == 1) {
            warning(paste("Maybe you have too little RAM.",
                            "Proceeding with SerialParam()"))
            return(BiocParallel::SerialParam(progressbar = TRUE))
        }
        return(BiocParallel::SnowParam(nwork, progressbar = TRUE))
    } else {
        return(BiocParallel::MulticoreParam(
            BiocParallel::multicoreWorkers() - 2,
            progressbar = TRUE)
        )
    }
}

#'@title RHermesMeta
#'@description The metadata storage class. It holds the experimental parameters,
#' all the names of the processed files and the timestamps for all operations
#' performed on the parental RHermesExp object
#'@slot ExpParam  Contains all experimental information
#'(ppm error, resolution,
#' polarity, etc.), as well as the formula and adduct databases used.
#'  See [RHermes]{ExpParam} for more info.
#'@slot filenames All filenames of the processed files. This info
#' is also
#' available for individual PL and SOI objects in their respective slot.
#'@slot timestamps Timestamps generated when any operation was
#'performed on the
#' parental object. You can easily check them with readTime().
#'@slot cluster Selected automatically based on your operating
#'system and your
#' number of cores. Can be set to any valid BiocParallelParam.
RHermesMeta <- setClass("RHermesMeta",
    slots = list(
        ExpParam = "ExpParam",
        filenames = "character",
        timestamps = "character",
        cluster = "BiocParallelParam"
    ), prototype = list(
        ExpParam = ExpParam(),
        filenames = character(0),
        timestamps = c(paste("System info:",
                    paste(Sys.info(), collapse = "/")),
                    paste("RHermes version:",
                    packageVersion("RHermes"))),
        cluster = RHermes:::setRHermesCluster()
    )
)

#'@title RHermesData
#'@description The experimental data storage class. Holds all information of the
#' peaklists, SOI lists, inclusion lists, MS2 data and identifications.
#'@slot PL A list object that holds all RHermesPL objects
#'containing a PL each.
#' See [RHermes]{RHermesPL} for more info
#'@slot SOI List that holds RHermesSOI objects with a SOI list
#' inside each.
#' See [RHermes]{RHermesSOI} for more info.
#'@slot MS2Exp List to hold inclusion lists, MSMS experimental
#'data and compound identifications.
RHermesData <- setClass("RHermesData",
                        slots = list(PL = "list", SOI = "list",
                                        MS2Exp = "list")
)

#' @title RHermesExp
#' @description  The main RHermes class: The RHermesExp. It is a container for
#' all generated information. All main RHermes functions use it and return an
#' updated version of the object.
#'
#' @slot metadata Where all the complementary info is stored
#' (experimental parameters, timestamps, databases, etc.).
#' See [RHermes]{RHermesMeta} for more info.
#' @slot data All experimental info is stored in data. It is
#' divided into PL (peaklist) SOI (scans of interest) and MS2Exp
#' (for the IL, MS2 data and identifications).
#' See [RHermes]{RHermesData} for more info.
#' @export RHermesExp
#' @examples
#' if(FALSE){
#'  myHermes <- RHermesExp() #Initializing empty object
#' }
RHermesExp <- setClass("RHermesExp",
                        slots = list(metadata = "RHermesMeta",
                                        data = "RHermesData")
)




