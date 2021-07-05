context("Peaklist generation works fine")

test_that("Raw data can be loaded",{
  ms1data <- RHermes:::import_and_filter(lf = system.file("extdata",
                                                          "MS1TestData.mzML",
                                                          package = "RHermes"))
  expect_length(ms1data, 3)
  expect_equal(nrow(ms1data[[1]]), 40040)
})

test_that("ScanSearch works",{
    BiocParallel::register(BiocParallel::SerialParam())
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- processMS1(myHermes, system.file("extdata",
                                                "MS1TestData.mzML",
                                                package = "RHermes"))
    #Same result as precalculated version
    expect_equal(nrow(myHermes@data@PL[[2]]@peaklist),
                 nrow(myHermes@data@PL[[1]]@peaklist))
})

test_that("Labelled proc works",{
    BiocParallel::register(BiocParallel::SerialParam())
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- processMS1(myHermes,
                            system.file("extdata",
                                        "MS1TestData.mzML",
                                        package = "RHermes"),
                            labelled = TRUE)

    expect_equal(nrow(myHermes@data@PL[[2]]@peaklist), 1378)
})

test_that("PL plot works", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))

    p <- RHermes::plotPL(myHermes, 1, "C3H7NO2", rtrange = c(0,1500),
                              dynamicaxis = TRUE, ads = NA)
    expect_true(is(p, "plotly"))
})


test_that("Coverage plot works", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- processMS1(myHermes,
                            system.file("extdata",
                                        "MS1TestData.mzML",
                                        package = "RHermes"),
                            labelled = TRUE)
    p <- RHermes:::plotCoverage(myHermes, 2)
    expect_true(is(p[[1]], "plotly") & is(p[[2]], "plotly"))
})

