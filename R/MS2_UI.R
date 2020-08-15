#' @export
MS2_UI <- function(id){
  ns <- NS(id)
  tagList(
    useSweetAlert(),
    verticalLayout(
      sidebarPanel(
        h2("MS2 File Input", style = "text-align: center;"),
        hr(style = "border-top : 1px dashed #C9B29B"),
        shinyFilesButton(ns("files"), "Enter the mzML file adresses:", "Select",
                         TRUE, icon = icon("database"), style = "margin-bottom: 10px"),
        verbatimTextOutput(ns("selecteddir"), placeholder = TRUE),
        hr(), width = 13),

      sidebarPanel(
        h2("MS2 processing", style = "text-align: center;"),
        hr(),
        numericInput(ns("costhr"), "Minimum cosine score", value = 0.8, min = 0,
                     max = 1, step = 0.01),
        selectInput(ns("MS2database"), "MS2 reference spectral database",
                    choices = c("MassBank EU"), selected = "MassBank EU"),
        uiOutput(ns("whenIL")),
        width = 13
      )
    )
  )

}

#' @export
MS2Server <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){
      roots <- getVolumes()
      shinyFileChoose(
        input,
        'files',
        roots = getVolumes(),
        filetypes = c("mzml", "mzxml")
      )
      output$selecteddir <- renderPrint({parseFilePaths(roots, input$files)})

      ns <- session$ns
      observeEvent(struct$hasIL,{
        output$whenIL <- renderUI(tagList(
          radioButtons(ns("selectIL"), "Select the IL from which the MS2 files were generated:",
                       choices = seq_along(struct$dataset@data@MS2Exp), selected = "1"),
          actionButton(ns("startMS2"), "Start MS2 data processing")
        ))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

      toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
      observeEvent(input$startMS2, {
        db <- switch(input$MS2database,
               "MassBank EU" = system.file("extdata","sp_MassBankEU_20200316_203615.RData", package = "RHermes")
               )
        if(db == ""){
          sendSweetAlert(
            session = session,
            title = "MS2 data couldn't be processed",
            text = "We couldn't find the specified database in the extdata folder",
            type = "error"
          )
        } else {
          toReturn$dataset <- MS2Proc(struct = struct$dataset, id = as.numeric(input$selectIL),
                                      MS2files = unname(parseFilePaths(roots, input$files)$datapath),
                                      referenceDB = db,
                                      mincos = input$costhr)
          toReturn$trigger <- toReturn$trigger + 1
          sendSweetAlert(
            session = session,
            title = "MS2 data processed",
            text = "Check the identifications in the corresponding tab",
            type = "success"
          )
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)



      return(toReturn)
    }
  )
}
