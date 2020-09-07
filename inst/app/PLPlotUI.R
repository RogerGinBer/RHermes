PLPlotUI <- function(id){
  ns <- NS(id)
  tagList(
    tags$br(),
    dropdown(circle = FALSE, status = "info",
             label = HTML('<span style = "margin-left:10px; margin-right:10px">Settings</span>'),
                   tooltip = tooltipOptions(title = "Plotting settings"),
                   icon = icon("align-justify"),
                   sidebarPanel(fluidRow(
      column(3, offset = 0,
             radioButtons(ns("PLplotmode"), label = "Select mode:",
                          c("By compound name", "By formula"), selected = "By compound name"),
             conditionalPanel("input.PLcompselect != 'No file detected'", ns = ns,
                              radioButtons(ns("PLfiles"), "Select which PL to use: ", choices = "", 1, inline = TRUE))
      ),
      column(9, offset = 0,
             fluidRow(
               column(8,
                      conditionalPanel(condition = "input.PLplotmode == 'By compound name'", ns = ns,
                                       selectizeInput(inputId = ns("PLcompselect"), label ="Select compound",
                                                      choices = "No file detected", selected = NULL, multiple = FALSE),
                      ),
                      conditionalPanel(condition = "input.PLplotmode == 'By formula'", ns = ns,
                                       selectizeInput(inputId = ns("PLformselect"), label = "Select formula",
                                                      choices = "No file detected", selected = NULL, multiple = FALSE),
                      )),

               column(3, offset = 1, tags$b("Dynamic axis"), switchInput(
                 inputId = ns("dynamicaxis"),
                 onStatus = "success",
                 offStatus = "danger",
                 value = TRUE, size = "small"
               ))
             ),
             conditionalPanel("input.PLcompselect != 'No file detected'", ns = ns,
                              sliderInput(ns("RTinterval"), label = "Select an RT interval:",
                                          min = 0, max = 1800, value = c(0,1800)),
                              pickerInput(
                                inputId = ns("ads"),
                                label = "Select adducts to plot",
                                choices = NULL,
                                options = list(
                                  `actions-box` = TRUE),
                                multiple = TRUE
                              ),
                              verbatimTextOutput(outputId = ns("othercomp")),
                              tags$head(tags$style("#PLPlotUI-othercomp{overflow-y:scroll; max-height: 250px; background: ghostwhite;}")))
      )
    ),width = 12), width = "100%"),
    plotlyOutput(outputId = ns("PLPlot"), height = "800px")
  )
}
PLPlotServer <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){
      ns <- session$ns
      observeEvent({
        struct$hasPL
        struct$dataset@data@PL
      },{
        if(struct$hasPL){
          peakLists <- seq_along(struct$dataset@data@PL)
          nameList <- strsplit(struct$dataset@metadata@filenames, split = "/")
          nameList <- lapply(nameList, function(x){
            return(x[length(x)])
          }) %>% unlist()
          names(peakLists) <- nameList

          updateSelectizeInput(session, "PLcompselect", choices = struct$dataset@metadata@ExpParam@DB$Name,
                               selected = struct$dataset@metadata@ExpParam@DB$Name[1],
                               server = TRUE)
          updateSelectizeInput(session, "PLformselect", choices = struct$dataset@metadata@ExpParam@DB$MolecularFormula,
                               selected = struct$dataset@metadata@ExpParam@DB$MolecularFormula[1],
                               server = TRUE)
          updateRadioButtons(session, "PLfiles", choices = peakLists)
          updateSliderInput(session, "RTinterval", max = ceiling(max(struct$dataset@data@PL[[1]]@raw$rt)),
                            value = c(0, ceiling(max(struct$dataset@data@PL[[1]]@raw$rt))))
          updatePickerInput(session, "ads", choices = struct$dataset@metadata@ExpParam@adlist$adduct,
                            selected = struct$dataset@metadata@ExpParam@adlist$adduct)
        }else{
          updateSelectizeInput(session, "PLcompselect", choices = "No file detected")
          updateSelectizeInput(session, "PLformselect", choices = "No file detected")
          updateRadioButtons(session, "PLfiles", choices = "")
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 100)

      observeEvent({
        input$PLcompselect
        input$PLformselect
        input$PLfiles
        input$PLplotmode
        input$RTinterval
        input$dynamicaxis
        input$ads
      },{
        if(struct$hasPL){
          if(input$PLplotmode == "By compound name"){
            if(input$PLcompselect == ""){return()}
            f <- with(struct$dataset@metadata@ExpParam@DB,{MolecularFormula[which(Name == input$PLcompselect)[1]]})[1]
          } else {
            if(input$PLformselect == ""){return()}
            f <- input$PLformselect
          }
          other <- with(struct$dataset@metadata@ExpParam@DB, {Name[MolecularFormula == f]})
          output$othercomp <- renderText(paste(other, collapse = "\n"))
          output$PLPlot <- renderPlotly(RHermes::PLPlot(struct = struct$dataset,
                                                                id = as.numeric(input$PLfiles),
                                                                formula = f, rtrange = as.numeric(input$RTinterval),
                                                                dynamicaxis = as.logical(input$dynamicaxis),
                                                                ads = as.character(input$ads)
                                                              ))
        }


      }, ignoreNULL = TRUE, ignoreInit = TRUE)
      return()
    }
  )
}
