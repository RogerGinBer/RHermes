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
MSMScosineSim <- function(pattern, query, minhits = 1, mode = "full") {
    if (nrow(pattern) == 0 | nrow(query) == 0) {
        stop("Invalid pattern or query")
    }
    if (ncol(query) == 1) {query <- t(query)}
    patint <- c()
    qint <- c()
    allmz <- unique(pattern$mz)
    pattern$int <- sqrt(pattern$int)
    query[, 2] <- sqrt(query[, 2])
    for (i in seq_along(allmz)) {
        m <- allmz[i]
        dist <- abs(query[, 1] - m)
        if (dist[which.min(dist)] < 0.001) {
            patint <- c(patint, pattern[i, 2])
            qint <- c(qint, query[which.min(dist), 2])
        }
    }
    if (length(patint) < minhits | length(qint) < minhits) {return(0)}
    if (mode == "full") {
        denom <- sqrt(sum(pattern$int^2)) * sqrt(sum(query[, 2]^2))
    } else {
        denom <- sqrt(sum(patint^2)) * sqrt(sum(qint^2))
    }
    return(sum(patint * qint)/denom)
}






