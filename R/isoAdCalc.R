isoCos <- function(soilist, PL, isothr = 0.99, BiocParallelParam) {
    PL <- PL@peaklist
    setkey(PL, formv)
    message("Calculating isotope similarity:")
    clist <- bplapply(seq_len(nrow(soilist)), parallelIsoCos, soilist, PL,
                      BPPARAM = BiocParallelParam)
    hits <- lapply(clist, function(x) {return(x[[1]])})
    hitdf <- lapply(clist, function(x) {return(x[[2]])})
    soilist$isofound <- as.numeric(hits)
    soilist$isodf <- hitdf
    return(soilist)
}

adCos <- function(soilist, FATable, adthr = 0.8, BiocParallelParam) {
    message("Calculating adduct similarity:")
    parseanot <- lapply(soilist$anot, function(x) {
        strsplit(x = x, split = "_")
    })
    soilist$f <- lapply(parseanot, function(x) {
        vapply(x, function(y) {y[[1]]}, character(1))
    })
    soilist$ad <- lapply(parseanot, function(x) {
        vapply(x, function(y) {y[[2]]}, character(1))
    })
    setkey(FATable, "f")
    clist <- bplapply(seq_len(nrow(soilist)), parallelAdCos, soilist, FATable,
                      BPPARAM = BiocParallelParam)
    soilist$adrows <- clist
    return(soilist)
}

parallelIsoCos <- function(i, soilist, PL){
    SOI <- soilist[i, ]
    rti <- SOI$start[[1]]
    rtend <- SOI$end[[1]]
    cur <- PL[.(SOI$formula[[1]])]
    cur <- cur[cur$rt >= rti & cur$rt <= rtend, ]
    setkey(cur, isov, rt)
    count <- 0
    hitdf <- data.frame(iso = character(), mass = numeric(),
                        stringsAsFactors = FALSE)
    if (length(which(cur$isov == "M0")) < 5) {
        return(list(count, hitdf))
    }
    for (j in unique(cur$isov)) {
        if (j == "M0") {next}
        cos <- RHermes:::cosineSim(pattern = cur[which(cur$isov == "M0"), ],
                                   query = cur[which(cur$isov == j),])
        if (!is.na(cos) & cos > isothr) {
            count <- count + 1
            hitdf <- rbind(hitdf,
                           data.frame(iso = j,
                                      mass = mean(cur$mz[cur$isov == j, ]),
                                      stringsAsFactors = FALSE))
        }
    }
    return(list(count, hitdf))
}

parallelAdCos <- function(i, soilist, FATable){
    SOI <- soilist[i, ]
    rti <- SOI$start[[1]]
    rtend <- SOI$end[[1]]
    f <- SOI$f[[1]]
    equiv <- lapply(f, function(x) {FATable[.(x)]})
    equiv <- lapply(seq_along(f), function(x) {
        equiv[[x]][equiv[[x]]$ion != SOI$formula, ]
    })
    ids <- lapply(equiv, function(entries) {
        lapply(entries$ion, function(x) {
            candidates <- which(soilist$formula == x &
                                soilist$start >= (rti - 10) &
                                soilist$end <= (rtend + 10))
            if (length(candidates) == 0) {
                return()
            } else {
                return(
                    lapply(candidates, function(row){
                        score <- cosineSim(pattern = SOI$peaks[[1]],
                                           query = soilist$peaks[row][[1]],
                                           nscans = 5)

                        if (score > adthr) {return(row)}
                        else {return()}
                    })
                )
            }
        }) %>% unlist()
    })
    names(ids) <- unlist(SOI$ad)
    return(ids)
}
