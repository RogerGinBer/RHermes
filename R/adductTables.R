#'@export
adductTables <- function(ch_max = 1, mult_max = 1) {
    data(adducts, package = "enviPat")
    adducts$Mass[49] <- adducts$Mass[49] * (-1)  #Fixed wrong one
    negative.envi <- adducts[which(adducts$Ion_mode == "negative"),
        ]
    negative.envi <- negative.envi[which(negative.envi$Charge %in%
        as.character(seq(from = -1, to = -ch_max)) & negative.envi$Mult %in%
        seq(mult_max)), ]
    positive.envi <- adducts[which(adducts$Ion_mode == "positive"),
        ]
    positive.envi <- positive.envi[which(positive.envi$Charge %in%
        as.character(seq(ch_max)) & positive.envi$Mult %in% seq(mult_max)),
        ]
    negative.ad <- negative.envi[, -c(2, 9)]
    colnames(negative.ad)[c(1, 4)] <- c("adduct", "massdiff")
    positive.ad <- positive.envi[, -c(2, 9)]
    colnames(positive.ad)[c(1, 4)] <- c("adduct", "massdiff")
    return(list(negative.ad, positive.ad))
}
