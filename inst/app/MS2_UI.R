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
        materialSwitch(inputId = ns("useDB"), label = "Perform matching against an MS2 database?", value = FALSE),
        conditionalPanel(condition = "input.useDB",
            shinyFilesButton(ns("MS2database"), "Enter an MS2ID rds file:", "Select",
                    FALSE, icon = icon("database"), style = "margin-bottom: 10px")
                         , ns = ns),
        uiOutput(ns("whenIL")),
        width = 13
      )
    )
  )

}

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

      shinyFileChoose(
        input,
        'MS2database',
        roots = getVolumes(),
        filetypes = c("rds")
      )

      ns <- session$ns
      observeEvent(struct$hasIL,{
        output$whenIL <- renderUI(tagList(
          radioButtons(ns("selectIL"), "Select the IL from which the MS2 files were generated:",
                       choices = seq_along(struct$dataset@data@MS2Exp), selected = "1"),
          uiOutput(ns("startMS2_box"))
        ))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

      observeEvent(input$files , {
        if(nrow(parseFilePaths(roots,input$files)) == 0){
          output$startMS2_box <- renderUI({})
        } else {
          output$startMS2_box <- renderUI({
            tags$div(actionButton(ns("startMS2"), "Start MS2 processing",
                                  style = "background-color: #4d4263; color: #F0F0F0"),
                     style = "text-align: center;")})
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE
      )


      toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
      observeEvent(input$startMS2, {
        if(length(parseFilePaths(roots, input$MS2database)$datapath) == 0 & input$useDB){
          sendSweetAlert(
            session = session,
            title = "MS2 data couldn't be processed",
            text = "No DB selected",
            type = "error"
          )
        } else {
          toReturn$dataset <- processMS2(struct = struct$dataset, id = as.numeric(input$selectIL),
                                      MS2files = unname(parseFilePaths(roots, input$files)$datapath),
                                      referenceDB = parseFilePaths(roots, input$MS2database)$datapath,
                                      useDB = input$useDB, sstype = "regular",
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
