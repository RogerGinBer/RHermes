context("Peaklist generation works fine")

test_that("Raw data can be loaded",{
  ms1data <- RHermes:::import_and_filter(lf = system.file("extdata",
                                                          "MS1TestData.mzML",
                                                          package = "RHermes"))
  expect_length(ms1data, 3)
  expect_equal(nrow(ms1data[[1]]), 142619)
})

test_that("ScanSearch works",{
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- RHermesExp()
  myHermes <- setDB(myHermes, db = "hmdb")
  myHermes@metadata@cluster <- BiocParallel::SnowParam(1)
  myHermes <- FileProc(myHermes, system.file("extdata",
                                             "MS1TestData.mzML",
                                             package = "RHermes"))
  expect_equal(length(myHermes@data@PL), 1)
  expect_equal(nrow(myHermes@data@PL[[1]]@peaklist), 13536)
})

test_that("Labelled proc works",{
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  myHermes <- RHermesExp()
  myHermes <- setDB(myHermes, db = "hmdb")
  myHermes@metadata@cluster <- BiocParallel::SnowParam(1)
  myHermes <- FileProc(myHermes, system.file("extdata",
                                             "MS1TestData.mzML",
                                             package = "RHermes"),
                       labelled = TRUE)
  expect_equal(length(myHermes@data@PL), 1)
  expect_equal(nrow(myHermes@data@PL[[1]]@peaklist), 33462)
})

test_that("PL plot works", {
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- RHermesExp()
  myHermes <- setDB(myHermes, db = "hmdb")
  myHermes@metadata@cluster <- BiocParallel::SnowParam(1)
  myHermes <- FileProc(myHermes, system.file("extdata",
                                             "MS1TestData.mzML",
                                             package = "RHermes"))
  p <- RHermes::PLPlot(myHermes, 1, "C6H12O6", rtrange = c(0,1500),
                              dynamicaxis = TRUE, ads = c("M+Na"))
  expect_true(is(p, "plotly"))
})


test_that("Coverage plot works", {
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- RHermesExp()
  myHermes <- setDB(myHermes, db = "hmdb")
  myHermes@metadata@cluster <- BiocParallel::SnowParam(1)
  myHermes <- FileProc(myHermes, system.file("extdata",
                                             "MS1TestData.mzML",
                                             package = "RHermes"))
  p <- RHermes:::coveragePlot(myHermes, 1)
  expect_true(is(p[[1]], "plotly") & is(p[[2]], "plotly"))
})

