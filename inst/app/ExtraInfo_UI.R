ExtraInfo_UI <- function(id){
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("History",
               sidebarPanel(
                 verbatimTextOutput(ns("histext"), placeholder = TRUE), width = 13
               )
      ),
      tabPanel("Compound database",
               sidebarPanel(
                   DT::dataTableOutput(ns("compound_df")), width = 13
               )
      ),
      tabPanel("Peak List",
               sidebarPanel(
                 div(radioButtons(inputId = ns("PLfiles"), label = "Select PL:",
                                  choices = "", selected = ""), style = "margin-left: 5%; margin-top: 3%"), width = "AUTO"),
               sidebarPanel(
                 plotlyOutput(outputId = ns("coveragePlot1")),
                 plotlyOutput(outputId = ns("coveragePlot2")), width = "AUTO")
      ),
      tabPanel("SOI List",
               sidebarPanel(div(radioButtons(inputId = ns("SOIfiles"), label = "Select SOI List:",
                                             choices = "", selected = ""), style = "margin-left: 5%; margin-top: 3%"),width = "AUTO"),
               sidebarPanel(DT::dataTableOutput(ns("SOItable"), width = "auto"), width = "AUTO")
      ),
      tabPanel("IL Data",
               sidebarPanel(
                 div(radioButtons(inputId = ns("ILfiles"), label = "Select inclusion list:",
                                  choices = "", selected = ""), style = "margin-left: 5%; margin-top: 3%"), width = "AUTO"),
               sidebarPanel(
                   plotlyOutput(outputId = ns("plotIL"),height = "800px"),
                   DT::dataTableOutput(ns("ILtable"), width = "auto"),
                   width = "AUTO")
      ),
      tabPanel("Acquired MS2 data",
               sidebarPanel(
                 div(radioButtons(inputId = ns("MS2files"), label = "Select inclusion list with MS2 data:",
                                  choices = "", selected = ""),
                     selectInput(ns("MS2ILentry"), "Select IL entry to check:", choices = "", selected = ""),
                     style = "margin-left: 5%; margin-top: 3%"), width = "AUTO"),
               sidebarPanel(
                 plotlyOutput(ns("MS2raw")),
                 plotlyOutput(ns("MS2selectedscan")),
                 dataTableOutput(ns("MS2header")), width = "AUTO"
               )
      ))
  )

}

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

      observeEvent({struct$dataset@metadata@ExpParam@DB},{
        if(is.data.frame(struct$dataset@metadata@ExpParam@DB)){
            output$compound_df <- DT::renderDataTable(struct$dataset@metadata@ExpParam@DB,
                                                      options = list(scrollX = TRUE,
                                                      autoWidth = TRUE))
        }
      })


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
        covPlots <- RHermes:::plotCoverage(struct$dataset,
                                           id = as.numeric(input$PLfiles))
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
        data <- struct$dataset@data@SOI[[as.numeric(input$SOIfiles)]]@SOIList[,-c("peaks")]
        data$anot <- paste(data$anot, sep = "@")
        output$SOItable <- DT::renderDataTable(data, options = list(scrollX = TRUE, autoWidth = TRUE))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)


      observeEvent({
        struct$hasIL
        struct$dataset@data@MS2Exp
      },{
        if(struct$hasIL){
          numIL <- seq_along(struct$dataset@data@MS2Exp)
          updateRadioButtons(session, "ILfiles", choices = numIL, selected = numIL[1])
        }
      })
      observeEvent({input$ILfiles},{
        data <- struct$dataset@data@MS2Exp[[as.numeric(input$ILfiles)]]@IL@IL
        if(any(colnames(data) == "XIC")){
            data <- data[,-c("XIC")]
        }
        output$ILtable <- DT::renderDataTable(data, options = list(scrollX = TRUE, autoWidth = TRUE))
        output$plotIL <- renderPlotly(plotIL(struct$dataset, as.numeric(input$ILfiles)))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)


      observeEvent({
        struct$hasMS2
        struct$dataset@data@MS2Exp
      },{
        if(struct$hasMS2){
          numMS2 <- which(vapply(struct$dataset@data@MS2Exp, function(x){length(x@MS2Data) != 0}, logical(1)))
          updateRadioButtons(session, "MS2files", choices = numMS2, selected = numMS2[1])
        }
      })

      observeEvent({
        input$MS2files
      },{
        if(struct$hasMS2){
          numMS2entry <- seq_along(struct$dataset@data@MS2Exp[[as.numeric(input$MS2files)]]@MS2Data)
          updateSelectInput(session, "MS2ILentry", choices = numMS2entry, selected = numMS2entry[1])
        }
      })

      observeEvent({
        input$MS2files
        input$MS2ILentry
      },{
        data <- struct$dataset@data@MS2Exp[[as.numeric(input$MS2files)]]@MS2Data[[as.numeric(input$MS2ILentry)]]
        if(length(data) != 0){
          output$MS2header <- renderDataTable(data[[1]], options = list(scrollX = TRUE, autoWidth = TRUE))
          output$MS2raw <- renderPlotly(ggplotly(ggplot(data[[2]]) +
                                                   geom_point(aes(x=rt, y=mz, color = int, shape = as.factor(ID)))))
        } else {
          output$MS2header <- renderDataTable(data.frame())
          output$MS2raw <- renderPlotly(ggplotly(ggplot()))
        }
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

      observeEvent({input$MS2header_rows_selected},{
        s <- input$MS2header_rows_selected
        data <- struct$dataset@data@MS2Exp[[as.numeric(input$MS2files)]]@MS2Data[[as.numeric(input$MS2ILentry)]]
        if(length(data) & length(s)){
          rts <- unlist(data[[1]][s, "retentionTime"])
          data[[2]] <- data[[2]][data[[2]]$rt %in% rts, ]
          output$MS2selectedscan <- renderPlotly(ggplotly(ggplot(data[[2]]) +
                                                   geom_segment(aes(x=mz, xend = mz, y=0, yend = int)) +
                                                   facet_grid(rt ~ .) +
                                                   theme_minimal()))
        } else {
          output$MS2selectedscan <- renderPlotly(ggplotly(ggplot()))
        }
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

    }

  )
}
