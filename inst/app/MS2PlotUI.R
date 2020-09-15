MS2PlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$br(),
    radioGroupButtons(inputId = ns("selectplot"), label = "Select a graph :",
                      choices = c(`<i class='fa fa-bar-chart'></i> Identifications` = "ident",
                                  `<i class="fas fa-database"></i> Raw MSMS data` = "raw"),
                      justified = TRUE),
    br(),
    
    conditionalPanel("input.selectplot == 'ident'", ns = ns,
                     fluidRow(
                       column(6,
                              selectizeInput(ns("selectident"),
                                             label = "Select the superspectra entry to check:",
                                             choices = NULL, selected = NULL)),
                       column(6,
                              pickerInput(ns("selecthits"),
                                          label = "Select which hits to check against:",
                                          choices = NULL, selected = NULL, multiple = TRUE,
                                          options = list(
                                            `actions-box` = TRUE))
                       )
                     ),
                     plotlyOutput(ns("mirrorplot"))
    ),
    conditionalPanel("input.selectplot == 'raw'", ns = ns,
                     fluidRow(
                       column(6,
                              radioButtons(ns("id"), "Select MS2Exp:", choices = "",
                                           selected = "", inline = TRUE)
                       ),
                       column(6,
                              checkboxGroupButtons(
                                inputId = ns("rawset"),
                                label = "Select what to plot:",
                                choices = c("Raw data plot", "Network", "Peak table"),
                                checkIcon = list(
                                  yes = tags$i(class = "fa fa-check-square",
                                               style = "color: #4D4263"),
                                  no = tags$i(class = "fa fa-square-o",
                                              style = "color: #4D4263")),
                                selected = c("Raw data plot", "Network", "Peak table")
                              )
                       )
                     ),
                     fluidRow(column(width = 6,
                                     selectizeInput(ns("selectILentry"),
                                                    label = "Select the IL entry to check:", choices = NULL,
                                                    selected = NULL
                                     )
                     ),
                     column(width = 6, tags$b("Select plot mode:"),
                            switchInput(ns("bymz"), onLabel = "By m/z",
                                        offLabel = "By group", value = TRUE,
                                        labelWidth = "AUTO"))),
                     fluidRow(
                       conditionalPanel("input.rawset.includes('Raw data plot')",
                                        ns = ns,
                                        conditionalPanel("input.bymz", ns = ns,
                                                         plotlyOutput(ns("rawMSMS_bymz"))),
                                        conditionalPanel("!input.bymz", ns = ns,
                                                         plotlyOutput(ns("rawMSMS_bygroup"))),
                       ),
                       conditionalPanel("input.rawset.includes('Network')",
                                        ns = ns, visNetworkOutput(ns("rawss"))),
                       conditionalPanel("input.rawset.includes('Peak table')",
                                        ns = ns,
                                        column(12, align="center",
                                               tableOutput(ns('rawpks'))))
                     )
    )
  )
  
}

MS2PlotServer <- function(id, struct) {
  moduleServer(id, function(input, output, session) {
    # observeEvent({},{}, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    ####Updates to buttons, selections, etc.####
    observeEvent({
      struct$hasMS2
      struct$dataset@data@MS2Exp
    }, {
      if (struct$hasMS2) {
        #Determine which have ms2 data
        whichMS2 <- vapply(struct$dataset@data@MS2Exp,
                           function(ms2) {
                             return(length(ms2@MS2Data) != 0)
                           }, logical(1))
        whichMS2 <- which(whichMS2)
        
        #Update accordingly
        updateRadioButtons(session, "id", choices = whichMS2,
                           selected = whichMS2[1])
        
      } else {
        #Hide panels again
        updateRadioButtons(session, "id", choices = NULL,
                           selected = NULL)
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 100)
    
    observeEvent({
      input$id
    }, {
      ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
      sslist <- ms2@Ident[["MS2Features"]]
      with_hits <- which(vapply(sslist$results,
                                function(x){is.data.frame(x)},
                                logical(1)))
      original_IL <- unique(sslist$ILentry)
      updateSelectizeInput(session, "selectident",
                           choices = as.character(with_hits), server = TRUE,
                           selected = as.character(with_hits[1]),
                           options = list(maxOptions = 50000))
      updateSelectizeInput(session, "selectILentry",
                           choices = as.character(original_IL),
                           selected = as.character(original_IL[1]),
                           server = TRUE, options = list(maxOptions = 50000))
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)
    
    
    #### Plots ####
    observeEvent({
      input$selectILentry
      input$id
    }, {
      if (input$selectILentry != "" & !is.na(input$selectILentry)) {
        rawplots <- RawMS2Plot(struct$dataset,
                               as.numeric(input$id), as.numeric(input$selectILentry),
                               bymz = input$bymz)
        output$rawMSMS_bymz <- renderPlotly(rawplots[["p_bymz"]])
        output$rawMSMS_bygroup <- renderPlotly(rawplots[["p_bygroup"]])
        output$rawpks <- renderTable(rawplots[["pks"]])
        if (!is.na(rawplots[[2]][1])) {
          output$rawss <- renderVisNetwork(rawplots[["net"]])
        } else {
          # output$rawss <- renderVisNetwork({})
        }
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -1)
    
    observeEvent({
      input$selectident
      input$id
    },{
      ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
      sslist <- ms2@Ident[["MS2Features"]]
      if(is.data.frame(sslist) & !is.na(as.numeric(input$selectident))){
        if(as.numeric(input$selectident) <= nrow(sslist)){
          cur <- sslist$results[[as.numeric(input$selectident)]]
          if(is.data.frame(cur)){
            updatePickerInput(session, "selecthits", choices = cur$formula,
                              selected = cur$formula[[1]])
          }
        }
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)
    
    observeEvent({
      input$selectident
      input$selecthits
      input$id
    }, {
      if (input$selectident != "" & !is.na(input$selectident) &
          length(input$selecthits) != 0) {
        ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
        sslist <- ms2@Ident[["MS2Features"]]
        with_hits <- which(vapply(sslist$results,
                                  function(x){is.data.frame(x)},
                                  logical(1)))
        if(as.numeric(input$selectident) %in% with_hits){
          identplots <- MirrorPlot(struct$dataset,
                                   as.numeric(input$id), as.numeric(input$selectident),
                                   input$selecthits)
          output$mirrorplot <- renderPlotly(identplots)
        }
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -2)
  })
}


