context("Functional MS2 identification")

test_that("cosineSim works", {
    cosine <- cosineSim(data.frame(rt = seq(1,10), rtiv = seq(1,10)^2),
                    data.frame(rt = seq(1,10), rtiv = seq(1,10)^3))
    expect_equal(cosine, 0.986387, tolerance = 1e-4)
})

test_that("Compounds are identified", {
    skip_on_bioc() #Depends on copyrighted databases and local MS2 files, so we skip it
    skip_if(length(list.files("D:/ABrunner Plasma/MS2data",
                            pattern = ".*pos.*.mzML", full.names = TRUE)) == 0)
    myHermes <- readRDS(system.file("extdata", "testIL.rds",
                                    package = "RHermes"))

    MS2files <- list.files("D:/ABrunner Plasma/MS2data",
                           pattern = ".*pos.*.mzML", full.names = TRUE)
    myHermes <- MS2Proc(myHermes, 1, MS2files,
                        referenceDB = "D:/MS2ID_20200824_202808.rds", useDB = TRUE)

    expect_equal(nrow(myHermes@data@MS2Exp[[1]]@Ident[[1]]),  15)
})

test_that("Superspectra can be exported", {
    myHermes <- readRDS(system.file("extdata", "withIdent.rds",
                                    package = "RHermes"))
    exportMGF(myHermes, 1, "test")
    exportMSP(myHermes, 1, "test")
    file.remove(c("test.mgf", "test.msp"))
    succeed()
})

test_that("Mirror plot works", {
    myHermes <- readRDS(system.file("extdata", "withIdent.rds",
                                    package = "RHermes"))
    p <- RHermes::MirrorPlot(myHermes, 1, 2, patform = 1, mode = "versus")
    expect_true(is(p, "plotly"))
})

test_that("Raw MS2 plot works", {
    myHermes <- readRDS(system.file("extdata", "withIdent.rds",
                                    package = "RHermes"))

    p <- RHermes::RawMS2Plot(myHermes, ms2id = 1, entryid = 22,
                                   bymz = TRUE)
    p2 <- RHermes::RawMS2Plot(myHermes, 1, 4, bymz = FALSE)
    expect_true(is(p[[1]], "plotly") & is(p[[2]], "visNetwork"))
    expect_true(is(p2[[1]], "plotly") & is(p2[[2]], "visNetwork"))
})
