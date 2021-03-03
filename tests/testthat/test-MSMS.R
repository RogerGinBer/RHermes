context("Functional MS2 identification")

test_that("cosineSim works", {
    cosine <- cosineSim(data.frame(rt = seq(1,10), rtiv = seq(1,10)^2),
                    data.frame(rt = seq(1,10), rtiv = seq(1,10)^3))
    expect_equal(cosine, 0.986387, tolerance = 1e-4)
})

test_that("Compounds are identified", {
    #Depends on a large MS2 local files, so we skip it
    skip_on_bioc()
    skip_on_cran()
    skip_if(length(list.files("E:/ABrunner Plasma/MS2data",
                            pattern = ".*pos.*.mzML", full.names = TRUE)) == 0)

    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    MS2files <- list.files("E:/ABrunner Plasma/MS2data",
                           pattern = ".*pos.*.mzML", full.names = TRUE)[1:9]
    myHermes <- processMS2(myHermes, 1, MS2files, sstype = "regular",
                           useDB = FALSE)
    expect_equal(nrow(myHermes@data@MS2Exp[[1]]@Ident[[1]]),  7)
})

test_that("Superspectra can be exported", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    exportMGF(myHermes, 1, "test")
    exportMSP(myHermes, 1, "test")
    exportmzML(myHermes, 1, "test")
    exportSIRIUS(myHermes, 1, "test")
    file.remove(c("./test.mgf", "./test.msp", "./test.mzML", "./test.ms"))
    succeed()
})


test_that("Raw MS2 plot works", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))

    p <- RHermes::plotRawMS2(myHermes, ms2id = 1, entryid = 2)
    succeed()
})
