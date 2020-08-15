#' @export
injectionPlanner <- function(IL, injections, maxover, byMaxInt = TRUE, 
    stats = TRUE, returnAll = FALSE) {
    # IL <- IL@IL #Select only IL slot
    if (returnAll) {
        injections <- 99999
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
    if (stats) {
        print(paste("Total entry coverage:", round((initial_entries - 
            nrow(IL))/initial_entries * 100, digits = 3), "%"))
        num <- vapply(plan, function(x) {
            return(nrow(x))
        })
        df <- data.frame(inj = seq_along(num), entries = num, 
            ID = "Entries per injection")
        df <- rbind(df, data.frame(inj = seq_along(num), entries = cumsum(num), 
            ID = "Total entries analyzed"))
        print(ggplot(df) + geom_line(mapping = aes(x = inj, y = entries, 
            color = ID), size = 1) + geom_point(mapping = aes(x = inj, 
            y = entries, color = ID), size = 2) + geom_hline(yintercept = initial_entries, 
            linetype = "dotted", size = 0.5) + xlab("Injection") + 
            ylab("Number of entries") + ggtitle("Evolution of IL total entry coverage and entries per injection") + 
            scale_color_discrete(breaks = c("darkred", "blue")))
        avsnr <- data.frame(snr = log10(vapply(plan, function(x) {
            mean(x$SNR)
        })))
        print(ggplot(avsnr) + geom_bar(aes(x = seq_along(nrow(avsnr)), 
            y = snr), stat = "identity") + xlab("Injection") + 
            ylab("Mean log10(SNR)") + ggtitle("SNR evolution over injections"))
    }
    return(plan)
}
