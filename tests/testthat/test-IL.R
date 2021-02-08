context("Inclusion List generation")

test_that("Regular IL generation works and can be exported",{
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds",
                                  package = "RHermes"))
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- generateIL(myHermes, 1, ILParam())
  expect_equal(nrow(myHermes@data@MS2Exp[[1]]@IL@IL), 16)

  myHermes <- generateIL(myHermes, 1, ILParam(filtermz = 0.1,
                                              priorization = "full"))
  expect_equal(nrow(myHermes@data@MS2Exp[[2]]@IL@IL), 85)
})

test_that("Prioritized IL generation works",{
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds",
                                  package = "RHermes"))
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)
  myHermes <- cleanSOI(myHermes, 1, 50000, TRUE)
  myHermes <- generateIL(myHermes, 1, ILParam(filtermz = 0.1,
                                         priorization = "yes", ad = "M+H"))
  expect_equal(nrow(myHermes@data@MS2Exp[[1]]@IL@IL), 24)

})

test_that("IL can be exported", {
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds",
                                  package = "RHermes"))
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- generateIL(myHermes, 1, ILParam())
  exportIL(myHermes, id = 1, folder = getwd(), maxOver = 5, sepFiles = FALSE)
  exportIL(myHermes, id = 1, folder = getwd(), maxOver = 5, sepFiles = TRUE)
  file.remove(c("./ExportedIL.csv",
                paste0(paste("./Injection", seq(1,5), sep = "_"), ".csv")))
  succeed()
})


