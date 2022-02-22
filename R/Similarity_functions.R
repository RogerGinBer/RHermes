# Column 1 -> RT, column 2 -> intensity
cosineSim <- function(pattern, query, nscans = 5) {
    #Only scans found on both pattern and query are selected
    sametime <- which(query[,1] %in% pattern[,1])
    if (length(sametime) < nscans) {return(0)}
    query <- query[sametime, ]

    sametime <- which(pattern[,1] %in% query[,1])
    if (length(sametime) < nscans) {return(0)}
    pattern <- pattern[sametime, ]

    pattern <- dplyr::distinct(pattern[order(pattern[,1]), ])
    query <- dplyr::distinct(query[order(query[,1]), ])
    dotprod <- sum(pattern[,2] * query[,2])
    scaling <- sqrt(sum(pattern[,2]^2)) * sqrt(sum(query[,2]^2))
    return(dotprod/scaling)
}


pearsonSim <- function(pattern, query, nscans = 5) {
    #Only scans found on both pattern and query are selected
    sametime <- which(query[,1] %in% pattern[,1])
    if (length(sametime) < nscans) {return(0)}
    query <- query[sametime, ]

    sametime <- which(pattern[,1] %in% query[,1])
    if (length(sametime) < nscans) {return(0)}
    pattern <- pattern[sametime, ]

    pattern <- dplyr::distinct(pattern[order(pattern[,1]), ])
    query <- dplyr::distinct(query[order(query[,1]), ])

    num <- sum((pattern[,2] - mean(pattern[,2])) * (query[,2] -
        mean(query[,2])))
    denom <- sqrt(sum((pattern[,2] - mean(pattern[,2]))^2) *
        sum((query[,2] - mean(query[,2]))^2))
    return(num/denom)
}

MSMScosineSim <-
    function (pattern, query, minhits = 0, mzdiff = 0.02, minint = 0.1,
                do.sqrt=FALSE) {
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

#' @import philentropy
calculate_MS2_distance <- function(P, Q, method = "cosine", minint = 0.1,
                                    minhits = 1){
    if(method == "cosine"){
        MSMScosineSim(P, Q, minint = minint, minhits = minhits)
    } else if (method %in% c("fidelity", "squared_chord", "topsoe",
                                "jaccard", "canberra")) {
        l <- match_spec(P, Q, minint = minint, minhits = minhits)
        if(length(l) == 0) return(0)
        p <- rbind(l[[1]]/sum(l[[1]]), l[[2]]/sum(l[[2]]))
        #small delta added to avoid 0log0 error in topsoe
        suppressMessages(philentropy::distance(p + 1e-12, method = method))
    } else {
        stop("Unvalid method provided")
    }
}




