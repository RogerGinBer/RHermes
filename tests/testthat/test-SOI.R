context("SOI are generated successfully")

test_that("A SOI param class can be created",{
  s <- RHermes:::SOIParam()
  expect_s4_class(s, "SOIParam")
  d <- getSOIpar()
  expect_s4_class(d, "SOIParam")
})

test_that("SOI generation works",{
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
  myHermes <- SOIfinder(myHermes, getSOIpar(), 1)
  expect_equal(nrow(myHermes@data@SOI[[1]]@SoiList), 165)
})

test_that("Blank substraction is configured and works",{
  skip_on_bioc()
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)
  library(reticulate)
  library(keras)
  skip_if(!py_available(initialize = TRUE))
  expect(py_module_available("keras"), failure_message = "No Keras")
  expect(py_module_available("tensorflow"), failure_message = "No TensorFlow")
  model <- load_model_tf(system.file("extdata", "model",
                                     package = "RHermes"))
  expect(is(model, "python.builtin.object"),
         failure_message = "Model doesn't load")
  set.seed(2)
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds" ,
                                  package = "RHermes"))
  myHermes@data@PL <- c(myHermes@data@PL, myHermes@data@PL)
  myHermes@data@PL[[1]]@peaklist$rtiv <- myHermes@data@PL[[2]]@peaklist$rtiv*
    rnorm(100)^2
  myHermes <- SOIfinder(myHermes, getSOIpar(), 1, 2)

  p <- RHermes::SoiPlot(struct = myHermes, id = 2, formula = "C2H7NO3S",
                              blankid = 2, rtrange = c(0,1500),
                              dynamicaxis = TRUE, ads = "M+H")
  expect_true(is(p, "plotly"))
})

test_that("SOI plot works", {
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds" ,
                                  package = "RHermes"))
  p <- RHermes::SoiPlot(myHermes, 1, "C6H12O6",
                              rtrange = c(0,1500), dynamicaxis = TRUE,
                              ads = "M+Na")
  expect_true(is(p, "plotly"))
})


context("SOI cleanup works")
test_that("SOI are filtered correctly", {
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds" ,
                                  package = "RHermes"))
  myHermes <- SOIcleaner(myHermes, 1, 50000, TRUE)
  expect_equal(nrow(myHermes@data@SOI[[1]]@SoiList), 2)
})


