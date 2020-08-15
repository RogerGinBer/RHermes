parallelInterpreter <- function(x, uf, cutoff, BinRes, id){
    do.call(rbind,
            lapply((seq_along(uf[[x]])) + id[x] + 2,
                function(i) {
                    bin <- BinRes[[i]]
                    interlist <- RHermes:::densityInterpreter(bin, cutoff)
                    if (length(interlist[[1]]) == 0) {
                      return()
                    }
                    df1 <- data.frame(start = interlist[[1]],
                                    end = interlist[[2]],
                                    formula = rep(i - 2,
                                                times = length(interlist[[1]]))
                    )
                    return(df1)
                }))}

densityInterpreter <- function(list, cutoff) {
    ls <- rbind(list, cutoff)
    bool <- apply(ls, 2, function(x) {
        return(x[1] >= x[2])
    })
    if (!any(bool)) {
        return(list(c(), c()))
    }

    good <- which(bool)
    if (length(good) == 1) {
        return(list(good, good + 1))
    }

    diff <- diff(good)
    start <- good[c(1, which(diff != 1) + 1)]
    end <- good[c(which(diff != 1), length(good))] + 1
    return(list(start, end))
}
