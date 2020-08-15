context("GUI works")

test_that("GUI server works", {
  library(shiny)
  testServer(system.file("app", package = "RHermes"), {
    stopApp()
  })
  succeed()
})

# test_that("GUI interface works",{
#   RHermes:::PL_UI("")
#   RHermes:::SOI_UI("")
#   RHermes:::IL_UI("")
#   RHermes:::MS2_UI("")
#   RHermes:::PLPlotUI("")
#   RHermes:::SOIPlotUI("")
#   RHermes:::MS2PlotUI("")
#   RHermes:::Settings_UI("")
#   RHermes:::ExtraInfo_UI("")
#   succeed()
# })
