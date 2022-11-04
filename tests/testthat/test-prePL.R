context("Formula and adduct table functionality")
test_that("Database importer works", {
    hmdb <- RHermes:::database_importer(template = "hmdb")
    expect_equal(nrow(hmdb), 59)

    norman <- RHermes:::database_importer(template = "norman")
    expect_equal(nrow(norman), 118)

    custom <- RHermes:::database_importer(
    template = "custom", filename = system.file("extdata", "hmdb.csv",
                                                package = "RHermes"))
    expect_equal(nrow(custom), 59)

    custom2 <- RHermes:::database_importer(
    template = "custom", filename = system.file("extdata", "norman.xlsx",
                                                package = "RHermes"))
    expect_equal(nrow(custom2), 118)

    skip_if_offline()
    kegg <- RHermes:::database_importer(template = "kegg_p",
                                      keggpath = "hsa00010")
    expect_equal(nrow(kegg), 26)
})

test_that("Adduct tables generate successfully",{
    ad <- RHermes:::adductTables(1,1)
    expect_equal(nrow(ad[[1]]), 11)
    expect_equal(nrow(ad[[2]]), 21)
    expect_equal(ad[[2]][1,"adduct"], "[M+H]+")
})


context("Ionic formulas")
test_that("Ionic formulas generate correctly", {
    BiocParallel::register(BiocParallel::SerialParam())
    hmdb <- RHermes:::database_importer(template = "hmdb")
    ad <- RHermes:::adductTables(1,1)
    colnames(hmdb)[c(2,3)] <- c("m","fms")
    ionf <- RHermes:::IonicForm(hmdb[1:5,], ad[[2]][1:5,])
    expect_length(ionf, 2)
    expect_equal(nrow(ionf[[1]]), 24)
})


context("Isotopic distribution calculation")
test_that("Isotopes generate nicely", {
    BiocParallel::register(BiocParallel::SerialParam())
    ppm <- 2
    minmass <- 80
    maxmass <- 1050
    noiselevel <- 1e3
    FWHM <- 120000
    ion <- "+"
    par <- ExpParam(ppm = ppm, res = FWHM, nthr = noiselevel,
                  minmz = minmass, maxmz = maxmass, ion = ion)
    hmdb <- RHermes:::database_importer(template = "hmdb", minmass = 50,
                                      maxmass = 100)
    ad <- RHermes:::adductTables(1,1)
    colnames(hmdb)[c(2,3)] <- c("m","fms")
    test <- RHermes:::IonicForm(hmdb[1:5,], ad[[2]][1:5, ])
    IC <- RHermes:::IsoCalc(test[[1]], FWHM = par@res, intTHR = 0.2, kTHR = 1)
    expect_length(IC, 2)
    expect_length(IC[[1]], 24)
    expect_equal(nrow(IC[[2]]), 6)
})

