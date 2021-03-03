Settings_UI <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      fluidRow(
        column(selectInput(inputId = ns("partype"),
                           label = "Select Parallelization type:",
                           choices = c("serial", "snow"), selected = "snow"),
               width = 4),
        column(numericInput(inputId = ns("coreN"), label = "Number of cores",
                            value = 1, min = 1,
                            max = BiocParallel::bpworkers(),
                            step = 1), width = 4),
        column(uiOutput(ns("warning")), width = 4)
      ),
      fluidRow(
        column(actionButton(ns("saveSet"), label = "Save Settings",
                            icon = icon("save")), width = 4)
      ),
      width = 12)
    )
}

SettingsServer <- function(id, struct){
  moduleServer(id, function(input, output, session) {

    if(Sys.info()["sysname"] != "Windows"){
      updateSelectInput(session = session, inputId = "partype",
                        choices = c(serial = "serial", snow = "snow",
                                    "multicore (FORK)" = "multicore"))
    }

    observeEvent({
          input$partype
        }, {
          if(input$partype == "serial"){updateNumericInput(session = session,
                                                           inputId = "coreN",
                                                           value = 1, max = 1)
          } else {updateNumericInput(session = session, inputId = "coreN",
                                     max = BiocParallel::bpworkers())}
        })

    toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
    observeEvent({
      input$saveSet
    }, {
      struct$dataset@metadata@cluster <- switch(input$partype,
             serial = BiocParallel::SerialParam(),
             snow = BiocParallel::SnowParam(workers = as.numeric(input$coreN)),
             multicore = BiocParallel::MulticoreParam(workers =
                                                        as.numeric(input$coreN))
      )
      toReturn$dataset <- struct$dataset
      toReturn$trigger <- toReturn$trigger + 1
    },ignoreInit = TRUE, ignoreNULL = TRUE)

    observeEvent({input$coreN}, {
      if((as.numeric(input$coreN)/ BiocParallel::bpworkers()) > 0.75){
        output$warning <- renderUI(
          p(tags$b("Warning:"),
            "Using more than 75% cores may cause undesired effects like system unstability",
            style = "background-color: #4d4263; color: #fff;
            padding: 15px 15px 15px 15px; text-align: center")
        )
      } else {
          output$warning <- renderUI(p(""))
      }})

    return(toReturn)
  })
}
