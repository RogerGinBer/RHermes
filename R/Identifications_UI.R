#' @export
Ident_UI <- function(id) {
    ns <- NS(id)
    tagList(
        conditionalPanel(condition = "input.id != 'nothing'",
        ns = ns,
        sidebarPanel(width = "AUTO",
            fluidRow(
                column(width = 4,
                    radioButtons(ns("id"), "Select MS2Exp:", choices = "nothing",
                    selected = "nothing", inline = TRUE)),
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
                ))),
        sidebarPanel(width = "AUTO", dataTableOutput(ns("ms2table")))),
                    conditionalPanel(condition = "input.id == 'nothing'",
                                    ns = ns,
                    sidebarPanel(width = "AUTO", h2("Load your MS2 data first",
                                style = "text-align: center;"))))
}

#' @export
IdentServer <- function(id, struct) {
    moduleServer(id, function(input, output, session) {
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

        savedf <- reactiveValues(data = NULL)

        observeEvent({
            input$id
        }, {
            ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
            identdf <- ms2@Ident[[1]][, -c("qMSMS", "patMSMS")]
            identdf$anot <- vapply(identdf$anot, function(x) {
                paste(x, collapse = " ")
            }, character(1))
            identdf$score <- round(identdf$score, digits = 3)
            output$ms2table <- renderDataTable(identdf, options = list(scrollX = TRUE))
            savedf$data <- identdf

        }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)

        output$savecsv <- downloadHandler(filename = "ms2data.csv",
            content = function(file) {
                write.csv(savedf$data, file)
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
