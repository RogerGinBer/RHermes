#' @title RHermesGUI
#' @description Calling this function starts the GUI in a new window on your
#' default browser.
#' @author Roger Gine
#' @import KEGGREST
#' @return Nothing, just starts the GUI
#' @examples
#' RHermesGUI()
#' @export
RHermesGUI <- function() {
    shiny::runApp(system.file("app", package = "RHermes"),
                    launch.browser = TRUE)
}
