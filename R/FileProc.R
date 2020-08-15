#'@title FileProc
#'@description Main function to process the input MS1 mzML files. It accesses the files
#'with mzR and processes them with RHermes internal functions. For each file a PL object
#'is generated inside the object_data_PL slot.
#'
#'@param struct RHermesExp S4 object to update with the processed data. IMPORTANT: The
#'objects needs to have the experimental parameters already set before calling this
#'function.
#'@param files Character vector of all the paths to the files to be processed.
#'@param labelled Logical, whether to check for all 13C isotopic signals. Still in development. Defaults to FALSE
#'@return An RHermesExp object with the processed PLs.
#'@examples
#'FileProc(myHermes, 'C:/myFolder/myFile.mzML', FALSE) #For regular use
#'FileProc(myHermes, 'C:/myFolder/myFile.mzML', TRUE) #For labelled data
#'@export

setGeneric("FileProc", function(struct, files, labelled = FALSE) {
    standardGeneric("FileProc")
})
setMethod("FileProc", c("RHermesExp", "character", "ANY"), function(struct,
    files, labelled = FALSE) {

    ppm <- struct@metadata@ExpParam@ppm
    noise <- struct@metadata@ExpParam@nthr

    message(paste0("Preprocessing steps, calculating possible ionic formulas ",
            "and isotopic distributions"))

    F_DB <- struct@metadata@ExpParam@DB[,c("MolecularFormula", "EnviPatMass")]
    F_DB <- dplyr::distinct_at(F_DB, "MolecularFormula",
            .keep_all = T)  #Can break if colname isn't exactly that
    colnames(F_DB) <- c("fms", "m")

    IF_DB <- IonicForm(F_DB, struct@metadata@ExpParam@adlist,
                        BiocParallelParam = struct@metadata@cluster)

    IC <- IsoCalc(
                IF_DB[[1]], FWHM = struct@metadata@ExpParam@res,
                intTHR = 0.02, kTHR = 1, instr = struct@metadata@ExpParam@instr,
                refm = 200, BiocParallelParam = struct@metadata@cluster)

    message("Starting file processing...")
    struct <- setTime(struct, "Started file processing into PL")
    toAdd <- lapply(seq_along(files), function(i) {
        lf = files[i]
        message(paste0("Now processing: ", lf))
        imported <- import_and_filter(lf, 20, noise)
        ss <- RHermes:::OptScanSearch(DB = IF_DB[[1]],
            raw = imported[[3]], mzList = imported[[2]], ppm = ppm,
            labelled = labelled, IsoList = IC,
            BiocParallelParam = struct@metadata@cluster)

        #Construction of S4 Object output
        RHermesPL(peaklist = ss, header = imported[[2]], raw = imported[[1]],
                    labelled = labelled, filename = lf)
    })
    struct@data@PL <- c(struct@data@PL, toAdd)
    struct@metadata@ExpParam@ionF <- IF_DB
    struct@metadata@ExpParam@isoList <- IC
    struct@metadata@filenames <- c(struct@metadata@filenames, files)

    struct <- setTime(struct,
                      paste("Ended file processing. The following files were",
                    "successfully processed:", paste(files, collapse = "  ")))
    return(struct)
})



import_and_filter <- function(lf, minpks = 20, noise = 1000) {
    fileml <- mzR::openMSfile(lf)  # opening the connection to a single mzML file
    plist <- mzR::peaks(fileml)
    h <- mzR::header(fileml)
    h <- h[, -which(vapply(h, function(x) all(x == 0), FUN.VALUE = logical(1)))]  #Filtering header
    if (any(h$peaksCount < minpks)) {
        plist <- plist[-which(h$peaksCount < minpks)]  #Removing scans with very very few peaks detected
        h <- h[-which(h$peaksCount < minpks), ]
    }
    raw <- lapply(seq_along(plist), function(x) {
        #Extracting raw data into a DT
        rpeaks <- plist[[x]]
        rt <- h$retentionTime[x]
        return(data.table(mz = rpeaks[, 1], rtiv = rpeaks[, 2],
            rt = rt))
    })
    raw <- do.call(rbind, raw)
    filtered <- raw[raw$rtiv > noise, ]
    return(list(raw, h, filtered))
}

