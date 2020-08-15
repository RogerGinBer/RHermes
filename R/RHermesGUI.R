#' @export
RHermesGUI <- function() {
    shiny::runApp(system.file("app", package = "RHermes"), launch.browser = TRUE)
}
