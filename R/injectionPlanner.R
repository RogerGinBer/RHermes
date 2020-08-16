injectionPlanner <- function(IL, injections, maxover, byMaxInt = TRUE,
                             returnAll = FALSE) {
    if (returnAll) {
        injections <- 1e4
    }
    if (byMaxInt) {
        IL <- IL[order(-IL$MaxInt), ]
    }

    idx <- which(is.na(IL$start) | is.na(IL$end))  #NA depuration
    if (length(idx) != 0) {
        IL <- IL[-idx, ]
    }
    initial_entries <- nrow(IL)
    plan <- list()
    while (nrow(IL) != 0 & injections > 0) {
        TimeInt <- seq(min(IL$start, na.rm = TRUE), max(IL$end,
            na.rm = TRUE), by = 0.5)
        OL <- rep(0, length(TimeInt))
        ok_entries <- c()
        for (i in seq_len(nrow(IL))) {
            timeidx <- which(TimeInt >= IL$start[i] & TimeInt <=
                IL$end[i])
            if (any(OL[timeidx] >= maxover)) {
                next
            } else {
                ok_entries <- c(ok_entries, i)
                OL[timeidx] <- OL[timeidx] + 1
            }
        }
        curinj <- IL[ok_entries, ]
        IL <- IL[-ok_entries, ]
        plan <- c(plan, list(curinj))
        injections <- injections - 1
    }
    if (nrow(IL) != 0) {
        warning(paste0("Some SOI haven't been added to the injection plan due",
            " to lack of space. Try again with more injections, more maxover",
            " or returnAll = T"))
    }
    return(plan)
}
