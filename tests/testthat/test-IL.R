context("Inclusion List generation")

test_that("Regular IL generation works",{
  myHermes <- readRDS(system.file("extdata", "afterSOI.rds",
                                  package = "RHermes"))
  library(BiocParallel)
  require(CHNOSZ)
  require(magrittr)
  library(data.table)
  require(tidyverse)

  myHermes <- genIL(myHermes, 1, ILParam())
  myHermes <- genIL(myHermes, 1, ILParam(filtermz = 0.1, priorization = "full"))


})

