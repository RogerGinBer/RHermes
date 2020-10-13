IonicForm <- function(F_DB, Ad_DB, BiocParallelParam = SerialParam()) {
    suppressWarnings({
        RES <- bplapply(seq_len(nrow(F_DB)), RHermes:::calculate_ionic_forms,
                        BPPARAM = BiocParallelParam, F_DB = F_DB, Ad_DB = Ad_DB)
    })
    db <- do.call(rbind, lapply(RES, function(x) {x[[1]]}))
    db <- db[!duplicated(db[, 1]), ]
    connections <- do.call(rbind, lapply(RES, function(x) {x[[2]]}))
    return(list(db, connections))
}

calculate_ionic_forms <- function(i, F_DB, Ad_DB){
    f <- as.character(F_DB$fms[i])
    j <- apply(Ad_DB, 1, function(x) {
        current_f <- f
        if (x[3] != 1) current_f <- RHermes:::multform(current_f,
                                                       as.numeric(x[3]))
        if (x[6] != "FALSE") current_f <- RHermes:::sumform(current_f, x[6])
        if (x[7] != "FALSE") current_f <- RHermes:::subform(current_f, x[7])
        if (is.na(current_f)) return(NA)

        ch <- ifelse(x[5] == "positive", "+", "-")
        current_f <- paste0("[", current_f, "]",
                            ifelse(abs(as.numeric(x[2])) == 1,
                                   ch,
                                   ifelse(as.numeric(x[2]) > 0,
                                          c(x[2], ch),
                                          c(strsplit(x[2], "-")[[1]][2], ch))))
        return(current_f)
    })
    good <- which(!is.na(j))
    j <- j[!is.na(j)]
    envi <- strsplit(j, split = "[", fixed = TRUE) %>%
        vapply(function(x) {x[[2]]}, FUN.VALUE = character(1)) %>%
        strsplit(j, split = "]", fixed = TRUE) %>%
        vapply(function(x) {x[[1]]}, FUN.VALUE = character(1))

    db <- data.table(f = j,
                     m = round((F_DB$m[i] * as.numeric(Ad_DB[good,3])/
                                    abs(as.numeric(Ad_DB[good, 2])))
                               + as.numeric(Ad_DB[good,4]), digits = 5),
                     ch = as.numeric(Ad_DB[good, 2]),
                     envi = envi)
    con <- data.table(f = rep(f, length(j)), ion = j, an = Ad_DB[good, 1])
    return(list(db, con))
    }

sumform <- function(f1, f2) {
    f1 <- tryCatch({
        CHNOSZ:::count.elements(f1)
    }, error = function(cond){return(NA)})
    if(is.na(f1)) return(NA)
    n1 <- names(f1)
    f2 <- CHNOSZ:::count.elements(f2)
    n2 <- names(f2)
    interAB <- n1[which(n1 %in% n2)]
    if (length(interAB) != 0) {
        s1 <- vapply(interAB, function(x) {
            f1[[x]] + f2[[x]]
        }, numeric(1))
        names(s1) <- interAB
    } else {
        s1 <- integer()
    }
    uniqueA <- n1[which(!(n1 %in% n2))]
    if (length(uniqueA) != 0) {
        s1 <- c(s1, f1[uniqueA])
    }
    uniqueB <- n2[which(!(n2 %in% n1))]
    if (length(uniqueB) != 0) {
        s1 <- c(s1, f2[uniqueB])
    }
    alln <- unique(c(n1, n2))
    alln <- alln[!(alln %in% c("C", "H"))]
    ord <- c()
    if ("C" %in% unique(c(n1, n2))) {
        ord <- c("C")
    }
    if ("H" %in% unique(c(n1, n2))) {
        ord <- c(ord, "H")
    }
    ord <- c(ord, sort(alln))
    s1 <- s1[ord]

    paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(s1), s1), collapse = "")
}

subform <- function(f1, f2) {
    f1 <- tryCatch({
        CHNOSZ:::count.elements(f1)
        }, error = function(cond){return(NA)})
    if(is.na(f1)){return(NA)}
    n1 <- names(f1)
    f2 <- CHNOSZ:::count.elements(f2)
    n2 <- names(f2)
    interAB <- n1[which(n1 %in% n2)]
    if (length(interAB) < length(n2)) {return(NA)}
    if (length(interAB) != 0) {
        s1 <- vapply(interAB, function(x) {
            res <- ifelse(f1[[x]] - f2[[x]] > 0, f1[[x]] - f2[[x]], NA)
        }, numeric(1))
        if (any(is.na(s1))) {return(NA)}
        if (any(s1 == 0)) {
            n1 <- n1[which(n1 %in% n2)[s1 != 0]]
            n2 <- n2[which(n2 %in% n1)[s1 != 0]]
            s1 <- s1[s1 != 0]
        }
        names(s1) <- interAB
    } else {
        s1 <- integer()
    }
    uniqueA <- n1[which(!(n1 %in% n2))]
    if (length(uniqueA) != 0) {
        s1 <- c(s1, f1[uniqueA])
    }
    uniqueB <- n2[which(!(n2 %in% n1))]
    if (length(uniqueB) != 0) {
        s1 <- c(s1, f2[uniqueB])
    }
    alln <- unique(c(n1, n2))
    alln <- alln[!(alln %in% c("C", "H"))]
    ord <- c()
    if ("C" %in% unique(c(n1, n2))) {
        ord <- c("C")
    }
    if ("H" %in% unique(c(n1, n2))) {
        ord <- c(ord, "H")
    }
    ord <- c(ord, sort(alln))
    s1 <- s1[ord]

    paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(s1), s1), collapse = "")
}

multform <- function(f, k) {
    f <- CHNOSZ:::count.elements(f)
    n <- names(f)
    f <- f * k
    paste0(mapply(function(x, y) {
        if (y == 1) {
            return(x)
        }
        return(paste0(x, y))
    }, names(f), f), collapse = "")
}

