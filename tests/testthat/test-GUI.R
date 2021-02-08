context("GUI works")

test_that("GUI server works", {
  library(shiny)
  testServer(system.file("app", package = "RHermes"), {
    stopApp()
  })
  succeed()
})
