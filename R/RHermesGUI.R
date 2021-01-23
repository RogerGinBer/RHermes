#' @import shinyWidgets
#' @import shinydashboard
#' @import visNetwork
#' @import KEGGREST
#' @import slickR
#' @import DT
#' @import shinyFiles

#' @export
RHermesGUI <- function() {
    shiny::runApp(system.file("app", package = "RHermes"), launch.browser = TRUE)
}
