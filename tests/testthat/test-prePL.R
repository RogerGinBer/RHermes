context("Formula and adduct table functionality")
test_that("The abridged HMDB formula df is loaded correctly", {
  hmdb <- RHermes:::database_importer(template = "hmdb")
  expect_equal(nrow(hmdb), 59)
})
test_that("The abridged NORMAN formula df is loaded correctly", {
  norman <- RHermes:::database_importer(template = "norman")
  expect_equal(nrow(norman), 106)
})

test_that("Adduct tables generate successfully",{
  ad <- RHermes:::adductTables(1,1)
  expect_equal(nrow(ad[[1]]), 10)
  expect_equal(nrow(ad[[2]]), 14)
  expect_equal(ad[[2]][1,"adduct"], "M+H")
})


context("Ionic formulas")
test_that("Ionic formulas generate correctly", {
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)

  hmdb <- RHermes:::database_importer(template = "hmdb")
  ad <- RHermes:::adductTables(1,1)
  colnames(hmdb)[c(2,3)] <- c("m","fms")
  ionf <- RHermes:::IonicForm(hmdb[1:5,], ad[[2]])
  expect_length(ionf, 2)
  expect_equal(nrow(ionf[[1]]), 70)
})


context("Isotopic distribution calculation")
test_that("Isotopes generate nicely", {
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  ppm <- 2
  minmass <- 80
  maxmass <- 1050
  noiselevel <- 1e3
  FWHM <- 120000 #Resolving power at full width half maximum at mz = 200 (RP = m/dm)
  ion <- "+"
  par <- ExpParam(ppm = ppm, res = FWHM, nthr = noiselevel,
                  minmz = minmass, maxmz = maxmass, ion = ion)
  hmdb <- RHermes:::database_importer(template = "hmdb", minmass = 50, maxmass = 100)
  ad <- RHermes:::adductTables(1,1)
  colnames(hmdb)[c(2,3)] <- c("m","fms")
  test <- RHermes:::IonicForm(hmdb[1:5,], ad[[2]][1:5, ])
  IC <- RHermes:::IsoCalc(test[[1]], FWHM = par@res, intTHR = 0.2, kTHR = 1)
  expect_length(IC, 2)
  expect_length(IC[[1]], 25)
  expect_equal(nrow(IC[[2]]), 5)
})
