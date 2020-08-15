context("Functional MS2 identification")

test_that("Compounds are identified", {
    skip_on_bioc() #Depends on local databases and MS2 files, so we skip it
    skip_if(length(list.files("D:/ABrunner Plasma/MS2data",
                            pattern = ".*pos.*.mzML", full.names = TRUE)) == 0)
    myHermes <- readRDS(system.file("extdata", "testIL.rds",
                                    package = "RHermes"))

    MS2files <- list.files("D:/ABrunner Plasma/MS2data",
                           pattern = ".*pos.*.mzML", full.names = TRUE)
    myHermes <- MS2Proc(myHermes, 1, MS2files,
                        referenceDB = "D:/sp_MS2ID_RogerGB.RData")

    expect_equal(nrow(myHermes@data@MS2Exp[[1]]@Ident[[1]]),  39)
})

test_that("Superspectra can be exported", {
    myHermes <- readRDS(system.file("extdata", "withIdent.rds",
                                    package = "RHermes"))
    exportMGF(myHermes, 1, "test")
    exportMSP(myHermes, 1, "test")
    file.remove(c("test.mgf", "test.msp"))
    succeed()
})
