#'@export
cosineSim <- function(pattern, query, nscans = 5) {
    #Only scans found on both pattern and query are selected
    sametime <- which(query$rt %in% pattern$rt)
    if (length(sametime) < nscans) {return(0)}
    query <- query[sametime, ]

    sametime <- which(pattern$rt %in% query$rt)
    if (length(sametime) < nscans) {return(0)}
    pattern <- pattern[sametime, ]

    pattern <- dplyr::distinct(pattern[order(pattern$rt), ])
    query <- dplyr::distinct(query[order(query$rt), ])
    dotprod <- sum(pattern$rtiv * query$rtiv)
    scaling <- sqrt(sum(pattern$rtiv^2)) * sqrt(sum(query$rtiv^2))
    return(dotprod/scaling)
}


#'@export
pearsonSim <- function(pattern, query, nscans = 5) {
    #Only scans found on both pattern and query are selected
    sametime <- which(query$rt %in% pattern$rt)
    if (length(sametime) < nscans) {return(0)}
    query <- query[sametime, ]

    sametime <- which(pattern$rt %in% query$rt)
    if (length(sametime) < nscans) {return(0)}
    pattern <- pattern[sametime, ]

    pattern <- dplyr::distinct(pattern[order(pattern$rt), ])
    query <- dplyr::distinct(query[order(query$rt), ])

    num <- sum((pattern$rtiv - mean(pattern$rtiv)) * (query$rtiv -
        mean(query$rtiv)))
    denom <- sqrt(sum((pattern$rtiv - mean(pattern$rtiv))^2) *
        sum((query$rtiv - mean(query$rtiv))^2))
    return(num/denom)
}

#' @export
MSMScosineSim <-
    function (pattern, query, minhits = 0, mzdiff = 0.05, minint = 0.1,
              do.sqrt=F) {
        if (nrow(pattern) == 0 | nrow(query) == 0) {
            stop("Invalid pattern or query")
        }
        if (is.null(dim(query))) {
            query <- matrix(query, ncol = 2, byrow = TRUE)
        }
        query <- query[query[, 2] > minint, , drop = FALSE]
        patint <- c()
        qint <- c()
        allmz <- unique(pattern$mz)
        allmz <- sort(c(pattern$mz,query[, 1]))
        if(do.sqrt){
            pattern$int <- sqrt(pattern$int)
            query[, 2] <- sqrt(query[, 2])
        }
        toignore <- c()
        for (i in seq_along(allmz)) {
            if(!i%in%toignore){
                m <- allmz[i]
                dist1 <- abs(query[, 1] - m)
                dist2 <- abs(pattern[, 1] - m)
                if( any(dist1 < mzdiff) | any(dist2 < mzdiff)   ){
                    seqmz <- which(abs(allmz-m)<mzdiff)
                    patint <- c(patint , sum(pattern[which(dist2 < mzdiff), 2]))
                    qint <- c(qint, sum(query[which(dist1 < mzdiff) , 2]))
                    toignore <- c(toignore,seqmz)
                }
            }
        }

        if (length(patint) < minhits | length(qint) < minhits) {
            return(0)
        }
        denom <- sqrt(sum(patint^2)) * sqrt(sum(qint^2))
        cos <- sum(patint * qint)/denom
        return(cos)
    }





