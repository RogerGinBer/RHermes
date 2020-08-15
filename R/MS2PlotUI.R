#' @export
MS2PlotUI <- function(id) {
    ns <- NS(id)
    tagList(tags$br(), radioGroupButtons(inputId = ns("selectplot"),
        label = "Select a graph :", choices = c(`<i class='fa fa-bar-chart'></i> Identifications` = "ident",
            `<i class="fas fa-database"></i> Raw MSMS data` = "raw"),
        justified = TRUE), br(), radioButtons(ns("id"), "Select MS2Exp:",
        choices = "", selected = "", inline = TRUE), conditionalPanel("input.selectplot == 'ident'",
        ns = ns, selectizeInput(ns("selectident"), label = "Select the identification entry to check:",
            choices = NULL, selected = NULL), plotlyOutput(ns("mirrorplot"))),
        conditionalPanel("input.selectplot == 'raw'", ns = ns,
            fluidRow(column(width = 6, selectizeInput(ns("selectILentry"),
                label = "Select the IL entry to check:", choices = NULL,
                selected = NULL)), column(width = 6, tags$b("Select plot mode:"),
                switchInput(ns("bymz"), onLabel = "By m/z", offLabel = "By group",
                  value = TRUE, labelWidth = "AUTO"))), fluidRow(plotlyOutput(ns("rawMSMS")),
                visNetworkOutput(ns("rawss")))))

}

#' @export
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
            entrieswithscans <- ms2@Ident[[4]]
            identdf <- ms2@Ident[[1]]
            identnames <- paste(identdf$compound, identdf$IL_ID,
                sep = "@")
            updateSelectizeInput(session, "selectident", choices = as.character(identnames),
                selected = as.character(identnames[1]), server = TRUE, options = list(maxOptions = 50000))
            updateSelectizeInput(session, "selectILentry", choices = as.character(entrieswithscans),
                selected = as.character(entrieswithscans[1]),
                server = TRUE, options = list(maxOptions = 50000))

        }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)


        #### Plots ####
        observeEvent({
            input$selectILentry
            input$id
            input$bymz
        }, {
            if (input$selectILentry != "" & !is.na(input$selectILentry)) {
                rawplots <- PlotlyRawMS2Plot(struct$dataset,
                  as.numeric(input$id), as.numeric(input$selectILentry),
                  bymz = input$bymz)
                output$rawMSMS <- renderPlotly(rawplots[[1]])
                if (!is.na(rawplots[[2]][1])) {
                  output$rawss <- renderVisNetwork(rawplots[[2]])
                } else {
                  output$rawss <- renderVisNetwork({})
                }
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -1)


        observeEvent({
            input$selectident
            input$id
        }, {
            if (input$selectident != "" & !is.na(input$selectident)) {
                ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
                identdf <- ms2@Ident[[1]]
                identnames <- paste(identdf$compound, identdf$IL_ID,
                  sep = "@")
                row <- which(identnames == input$selectident)
                identplots <- PlotlyMirrorPlot(struct$dataset,
                  as.numeric(input$id), row)
                output$mirrorplot <- renderPlotly(identplots)
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -2)
    })
}
