Ident_UI <- function(id) {
    ns <- NS(id)
    tagList(
        
        sidebarPanel(width = "AUTO", h2("Identification Table",
                                        style = "text-align: center;"),
                     hr(style = "border-top : 1px dashed #C9B29B"),
                     uiOutput(ns("identTable"))
        )
    )
}

IdentServer <- function(id, struct) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        ####Updates to buttons, selections, etc.####
        observeEvent({
            struct$hasMS2
            struct$dataset@data@MS2Exp
        }, {
            whichMS2 <- vapply(struct$dataset@data@MS2Exp,
                               function(ms2) {
                                   return(length(ms2@MS2Data) != 0)
                               }, logical(1))
            whichMS2 <- which(whichMS2)
            
            if(length(whichMS2) != 0){
                output$identTable <- renderUI({
                    tagList(fluidRow(
                        column(width = 4,
                               radioButtons(ns("id"), "Select MS2Exp:", choices = whichMS2,
                                            selected = whichMS2[1], inline = TRUE),
                               switchInput(ns("onlyhits"), label = "Only hits", value = TRUE)
                        ),
                        column(width = 4,
                               tags$b("Export identifications:"),
                               div(downloadButton(outputId = ns("savecsv"),
                                                  label = "Export as CSV", style = "margin-top: 10px"))),
                        column(width = 4, radioButtons(ns("format"),
                                                       label = "Superspectra export format:",
                                                       choices = c("msp", "mgf"), selected = "msp"),
                               shinySaveButton(id = ns("saveselector"),
                                               label = "Select folder",
                                               title = "Save superspectra", style = "margin-bottom: 10px"),
                               verbatimTextOutput(ns("savepath"), placeholder = TRUE),
                               actionButton(ns("savebutton"), "Save your super spectra",
                                            icon("save"),
                                            style = "display: block; margin: 0 auto;")
                        )),
                        
                        div(dataTableOutput(ns("ms2table")), style = "margin-top: 30px;")
                    )
                    
                })
            }}, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 10)
        
        savedf <- reactiveValues(data = NULL)
        
        observeEvent({
            input$id
            input$onlyhits
        }, {
            if(input$id != "nothing"){
                ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]@Ident[[1]]
                ms2 <- cbind(data.frame(ID = seq_len(nrow(ms2))), ms2) 
                
                if(input$onlyhits){
                    ms2 <- ms2[vapply(ms2$results, is.data.frame, logical(1)), ]
                }
                
                ms2$hits <- lapply(ms2$results, function(hits){
                    if(!is.data.frame(hits)){return(hits)}
                    hits$formula
                })
                ms2$bestscore <- lapply(ms2$results, function(hits){
                    if(!is.data.frame(hits)){return(NA)}
                    round(max(hits$cos), digits = 3)
                })
                ms2 <- dplyr::select(ms2, -c("ssdata", "results"))
                ms2$start <- as.numeric(ms2$start)
                ms2$end <- as.numeric(ms2$end)
                ms2$apex <- as.numeric(ms2$apex)
                ms2$precmass <- as.numeric(ms2$precmass)
                ms2$bestscore <- as.numeric(ms2$bestscore)
                
                output$ms2table <- renderDataTable(ms2, plugins = 'natural', server = FALSE,
                                                   options = list(columnDefs = list(list(type = "natural", targets = "_all")),
                                                                  scrollX = TRUE))
                ms2$start <- format(round(ms2$start,2),nsmall = 2)
                ms2$end <- format(round(ms2$end,2),nsmall = 2)
                ms2$apex <- format(round(ms2$apex,2),nsmall = 2)
                ms2$bestscore <- format(round(ms2$bestscore,4),nsmall = 4)
                ms2$precmass <- format(round(ms2$precmass,4),nsmall = 4)
                ms2$anot <- lapply(ms2$anot, function(x){paste(x, collapse = " ")})
                ms2$hits <- lapply(ms2$hits, function(x){paste(x, collapse = " ")})
                ms2$hits <- lapply(ms2$hits, function(x){gsub(pattern = "\n", replacement = "", x = x,)})
                ms2$hits <- lapply(ms2$hits, function(x){gsub(pattern = "\t", replacement = "", x = x,)})                
                savedf$data <- ms2
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE)
        
        output$savecsv <- downloadHandler(filename = "ms2data.csv",
                                          content = function(file) {
                                              data.table::fwrite(savedf$data, file)
                                          }, contentType = "text/csv")
        
        shinyFileSave(input, "saveselector", roots = getVolumes(),
                      filetypes = c(""))
        savepath <- reactive(as.character(parseSavePath(getVolumes(),
                                                        input$saveselector)$datapath))
        output$savepath <- renderText(savepath())
        observeEvent(input$savebutton, {
            if (length(savepath()) != 0) {
                
                switch(input$format,
                       msp = {exportMSP(struct$dataset, as.numeric(input$id), savepath())},
                       mgf = {exportMGF(struct$dataset, as.numeric(input$id), savepath())}
                )
                sendSweetAlert(session = session, title = "Saved",
                               text = paste("The file", savepath(), "has been saved successfully"),
                               type = "success")
            } else {
                sendSweetAlert(session = session, title = "Error",
                               text = paste("Please select a valid path"),
                               type = "warning")
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE)
    })
}

