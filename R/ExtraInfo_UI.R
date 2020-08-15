#' @export
ExtraInfo_UI <- function(id){
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("History",
                  sidebarPanel(
                    verbatimTextOutput(ns("histext"), placeholder = TRUE), width = 13
                    )
               ),
      tabPanel("Peak List",
                div(radioButtons(inputId = ns("PLfiles"), label = "Select PL:",
                             choices = "", selected = ""), style = "margin-left: 5%; margin-top: 3%"),
                plotlyOutput(outputId = ns("coveragePlot1")),
                plotlyOutput(outputId = ns("coveragePlot2"))
               ),
      tabPanel("SOI List",
               sidebarPanel(div(radioButtons(inputId = ns("SOIfiles"), label = "Select SOI List:",
                                choices = "", selected = ""), style = "margin-left: 5%; margin-top: 3%"),width = "AUTO"),
               sidebarPanel(DT::dataTableOutput(ns("SOItable"), width = "auto"), width = "auto")
               ),
      tabPanel("MS2 Data")

      )
    )

}

#' @export
ExtraInfoServer <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){
      observeEvent({
        struct$dataset
        struct$dataset@metadata@timestamps
      },
        {
          if(length(struct$dataset@metadata@timestamps) != 0){
              output$histext <- renderText(paste(struct$dataset@metadata@timestamps, collapse = "\n"))
          } else{
              output$histext <- renderText("There have been no modifications")
          }
          },ignoreInit = FALSE, ignoreNULL = TRUE
      )

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
          updateRadioButtons(session, "PLfiles", choices = peakLists, selected = peakLists[1])
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent({input$PLfiles},{
        covPlots <- RHermes:::coveragePlot(struct$dataset,
                                           entry = as.numeric(input$PLfiles))
        output$coveragePlot1 <- renderPlotly(covPlots[[1]])
        output$coveragePlot2 <- renderPlotly(covPlots[[2]])

      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent({
        struct$hasSOI
        struct$dataset@data@SOI
      },{
        if(struct$hasSOI){
        soiLists <- seq_along(struct$dataset@data@SOI)
        nameList <- lapply(struct$dataset@data@SOI, function(x){
           spl <- strsplit(x@filename[[1]], split = "/")[[1]]
          return(spl[length(spl)])
        }) %>% unlist()
        names(soiLists) <- nameList
        updateRadioButtons(session, "SOIfiles", choices = soiLists, selected = soiLists[1])
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent({input$SOIfiles},{
        data <- struct$dataset@data@SOI[[as.numeric(input$SOIfiles)]]@SoiList[,-c("peaks", "MLdata")]
      data$anot <- paste(data$anot, sep = "@")
      output$SOItable <- DT::renderDataTable(data, options = list(scrollX = TRUE, autoWidth = TRUE))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    }
  )
}
