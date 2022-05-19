PLPlotUI <- function(id){
  ns <- NS(id)
  tagList(
    tags$br(),
    radioGroupButtons(inputId = ns("selectclass"), label = "Select a graph :",
                      choices = c(`<i class='fa fa-bar-chart'></i> Compounds` = "comp",
                                  `<i class="fas fa-database"></i> Raw MS1 data` = "raw"),
                      justified = TRUE),
    conditionalPanel("input.selectclass == 'comp'", ns = ns,

        fluidRow(
            column(8, {plotlyOutput(outputId = ns("PLPlot"), height = "800px")}),
            column(4,
                sidebarPanel(
                         radioButtons(ns("PLplotmode"), label = "Select mode:",
                                      c("By compound name", "By formula"), selected = "By compound name"),
                         conditionalPanel("input.PLcompselect != 'No file detected'", ns = ns,
                                          radioButtons(ns("PLfiles"), "Select which PL to use: ", choices = "", 1, inline = TRUE)),
                                  conditionalPanel(condition = "input.PLplotmode == 'By compound name'", ns = ns,
                                                   selectizeInput(inputId = ns("PLcompselect"), label ="Select compound",
                                                                  choices = "No file detected", selected = NULL, multiple = FALSE),
                                  ),
                                  conditionalPanel(condition = "input.PLplotmode == 'By formula'", ns = ns,
                                                   selectizeInput(inputId = ns("PLformselect"), label = "Select formula",
                                                                  choices = "No file detected", selected = NULL, multiple = FALSE),
                                         tags$b("Dynamic axis"),
                                         switchInput(
                                             inputId = ns("dynamicaxis"),
                                             onStatus = "success",
                                             offStatus = "danger",
                                             value = TRUE, size = "small"
                                           )),
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
                                          tags$head(tags$style("#PLPlotUI-othercomp{overflow-y:scroll; max-height: 250px; background: ghostwhite;}"))
                                  )

                         , width = 12),
                    )
                )

        ),
    conditionalPanel("input.selectclass == 'raw'", ns = ns,
        dropdown(circle = FALSE, status = "info",
            label = HTML('<span style = "margin-left:10px; margin-right:10px">Settings</span>'),
            tooltip = tooltipOptions(title = "Plotting settings"),
            icon = icon("align-justify"),
            sidebarPanel(fluidRow(
                column(4, offset = 0,
                        radioButtons(ns("PLfiles_raw"), "Select which PL to use: ", choices = "", 1, inline = TRUE)
                  ),
                column(8,
                    numericInput(ns("targetmz"), "Target mz:", value = 118.08626, min = 0, max = 5000),
                    column(4, switchInput(ns("mzmetric"), "Use ppm?", value = TRUE, onLabel = "ppm", offLabel = "absolute")),
                    column(8, numericInput(ns("mztol"), "mz tolerance", value = 3, min = 0, max = 5000, step = 0.001)),
                    sliderInput(ns("RTinterval_raw"), label = "Select an RT interval:",
                                min = 0, max = 1800, value = c(0,1800))
                )
            ), width = 12),
            width = "100%"),
        hr(),
        plotlyOutput(outputId = ns("RawMS1Plot"), height = "1200px"),

    )
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
          updateRadioButtons(session, "PLfiles_raw", choices = peakLists)
          updateSliderInput(session, "RTinterval",
                            max = ceiling(max(struct$dataset@data@PL[[1]]@raw$rt, 1500)),
                            value = c(0, ceiling(max(struct$dataset@data@PL[[1]]@raw$rt, 1500))))
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
          output$PLPlot <- renderPlotly(RHermes::plotPL(struct = struct$dataset,
                                                                id = as.numeric(input$PLfiles),
                                                                formula = f, rtrange = as.numeric(input$RTinterval),
                                                                dynamicaxis = as.logical(input$dynamicaxis),
                                                                ads = as.character(input$ads)
                                                              ))
        }


      }, ignoreNULL = TRUE, ignoreInit = TRUE)


      observeEvent({
          input$RTinterval_raw
          input$targetmz
          input$mzmetric
          input$mztol
          input$PLfiles_raw
      },{
          if(struct$hasPL){
              d <- struct$dataset@data@PL[[as.numeric(input$PLfiles_raw)]]@raw
              if(is.data.frame(d) & nrow(d) > 0){
                    target <- as.numeric(input$targetmz)
                    if (input$mzmetric){
                        tol <- as.numeric(input$mztol) * target * 1e-6
                    } else {
                        tol <- as.numeric(input$mztol)
                    }
                    d <- filter(d, between(mz, target - tol, target + tol) &
                                  between(rt, as.numeric(input$RTinterval_raw[[1]]),
                                          as.numeric(input$RTinterval_raw[[2]])))
                    if(nrow(d) < 1e5){
                        p1 <- ggplotly(ggplot(d) +
                                           geom_point(aes(x=rt, y=mz, color = log10(rtiv)))+
                                           labs(color = "log10(Intensity)") +
                                           scale_colour_gradient(low = "#CFB7E5", high = "#15004D") +
                                           xlab("Retention time (s)") +
                                           ylab("m/z")
                        )
                        
                        p2 <- ggplotly(ggplot(d) +
                                         geom_point(aes(x=rt, y=log10(rtiv), color = mz)) +
                                         labs(color = "mz") +
                                         scale_colour_gradient(low = "#CFB7E5", high = "#15004D") +
                                         xlab("Retention time (s)") +
                                         ylab("Intensity")
                        )
                                       
                        output$RawMS1Plot <- renderPlotly(
                            subplot(p1, p2, nrows=2, shareX = TRUE)
                        )

                    } else {
                        warning("Too many datapoints in PlotRawMS1. ",
                                "Try a lower mz width")
                    } 
              }
          }


      }, ignoreNULL = TRUE, ignoreInit = TRUE)
      return()
    }
  )
}
