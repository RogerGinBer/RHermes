#' @export
status_UI <- function(id) {
    ns <- NS(id)
    tagList()
    
}

#' @export
statusServer <- function(id, struct) {
    moduleServer(id, function(input, output, session) {
        
    })
}
