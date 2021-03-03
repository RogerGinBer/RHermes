IO_UI <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarPanel(
            div(
                tags$span(h2("Save an RHermesExp RDS"), style = "display: inline-block; width: 90%"),
                tags$span(actionButton(class = "help", label ="", inputId = ns("param_help"),
                                       icon = icon("info-circle"))),
                style = "text-align: center"
            ),
            hr(),
            tags$label("Select which folder to use to save all resulting files:"),
            br(),
            fluidRow(style = "margin-top : 10px; margin-bottom : 10px",
                column(width = 3, shinySaveButton(id = ns("saveselector"),
                            label = "Select folder",
                title = "Select which folder to use to save all resulting files:")),
                column(width = 9, verbatimTextOutput(ns("savepath"), placeholder = TRUE))),
                actionButton(ns("savebutton"), "Save your RHermesExp object as an .RDS",
                         icon("save"), style = "display: block; margin: 0 auto;"),
            width = 13),

        sidebarPanel(
            h2("Load an RHermesExp RDS",
            style = "text-align: center;"),
            h4("(In case the files were processed in another session or using the RHermes functions directly)",
            style = "text-align: center;"),
            hr(),
            div(
                span(shinyFilesButton(ns("loadselector"),
                       "Select .RDS", "Select", multiple = FALSE), style = "margin-right: 10%"),
                span(actionButton(ns("loadbutton"), "Load .RDS")),
                style = "text-align: center"
            ),
            br(),
            verbatimTextOutput(ns("loadpath"), placeholder = TRUE),
            width = 13))
}

IOServer <- function(id, struct) {
    moduleServer(id, function(input, output, session) {

        ns <- session$ns
        #Save logic
        shinyFileSave(input, "saveselector", roots = getVolumes(),
            filetypes = c(""))
        savepath <- reactive(as.character(parseSavePath(getVolumes(),
                                                        input$saveselector)$datapath))
        output$savepath <- renderText(savepath())

        observeEvent(input$savebutton, {
            if (length(savepath()) != 0) {
                saveRDS(struct$dataset, paste0(savepath(), ".rds"))
                output$savepath <- renderText(paste("The file",
                  savepath(), "has been saved successfully"))
                sendSweetAlert(session = session, title = "Saved",
                               text = paste("The file", savepath(), "has been saved successfully"),
                               type = "success")
            } else{
                sendSweetAlert(session = session, title = "Error",
                               text = paste("Please select a valid path"),
                               type = "warning")
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

        #Load logic
        shinyFileChoose(input, "loadselector", roots = getVolumes(),
            filetypes = c("rds"))
        toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
        loadpath <- reactive(as.character(parseFilePaths(getVolumes(),
            input$loadselector)$datapath))
        observeEvent(input$loadbutton, {
            if (length(loadpath()) != 0) {
                toReturn$dataset <- readRDS(loadpath())
                output$loadpath <- renderText(paste("The file",
                  loadpath(), "has been loaded successfully"))
                toReturn$trigger <- toReturn$trigger + 1
                sendSweetAlert(session = session, title = "Loaded",
                               text = paste("The file", loadpath(), "has been loaded successfully"),
                               type = "success")
            } else{
                sendSweetAlert(session = session, title = "Error",
                               text = "Please select a valid path",
                               type = "warning")
            }
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

        observe({
            invalidateLater(1000 * 60)  #Refresh loadpath after 60 seconds independently if it changes
            output$loadpath <- renderText(loadpath())
        }, priority = -10)

        #Help Modals
        observeEvent(input$param_help, {
            showModal(ui = modalDialog(
                tagList(
                    tags$p("Since RHermes requires 2 separate MS experiments, it is sometimes necessary to save
                           your work and load it sometime later. Also, people who prefer to run RHermes using
                           regular R commands may want to visualize their data and results using this app."),
                    tags$p("For this purposes, here you have this Loading/Saving tab. To save, choose your desired
                           folder and filename using the interface. To load an RHermesExp object you have to select
                           the .rds file that contains it using the interface and then click \"Load .RDS\". In both
                           cases, an alert will notice you when the saving/loading process is finished.")
                    ),
                easyClose = TRUE, title = HTML('<h4 class = "help-header">Saving and loading your RHermesExp object</h4>')),
                session = session)
        })



        return(toReturn)
    })
}
