PL_UI <- function(id){
    ns <- NS(id)
    tagList(
        br(),
        useSweetAlert(),
        verticalLayout(
            sidebarPanel(
                HTML("<h2 style = 'text-align: center;'>MS1 File Input</h2>"),
                hr(),
                fluidRow(
                    column(width = 2,
                           shinyFilesButton(ns("files"), "Enter the mzML file adresses:",
                                "Select", TRUE, icon = icon("database"),
                                style = "margin-bottom: 10px"),),
                    column(width = 2,
                           materialSwitch(label ="Load an example mzML",
                                          status = "info",
                                          value = FALSE,
                                          inputId = ns("loadexample")),
                           offset = 7)
                ),
                tableOutput(ns("selecteddir")),
                width = 13),
            sidebarPanel(
                div(
                  tags$span(h2("Peaklist Parameters"),
                            style = "display: inline-block; width: 90%"),
                  tags$span(actionButton(class = "help", label ="",
                                        inputId = ns("param_help"),
                                        icon = icon("info-circle"))),
                  style = "text-align: center; align-items: center;"
                ),
                hr(),
                selectInput(ns("instrumentalmode"),
                            "You can select a preset here: ",
                            choices = c("orbi-pos", "orbi-neg",
                                        "qtof-pos", "qtof-neg"),
                            selected = "orbi-pos"),
                fluidRow(
                    column(4 ,
                    numericInput(ns("ppm"),
                                 "Select ppm error of the instrument: ",
                                 2, min = 0, max = 100, step = 0.1),
                    numericInput(ns("minMz"),"Minimum m/z:", 50),
                    numericInput(ns("maxMz"),"Maximum m/z", 1200)),
                    column(4,
                    numericInput(ns("noise"),"Instrument noise threshold: ",
                                 1000),
                    numericInput(ns("resolution"),
                                 "Instrument resolution at m/z = 200", 120000,
                                 step = 1000),
                    selectInput(ns("ion"),"Ionization mode",
                                choices = c("+", "-"), selected = "+")),
                    column(4,
                    selectInput(ns("instrument"), "Instrument type",
                                c("Orbitrap", "QTOF"), "Orbitrap"),
                    radioButtons(ns("labelled"), "Working with labelled data",
                                 choices = c("T","F"), selected = "F")
                    )
                ),
        hr(),
        fluidRow(column(4, selectInput(ns("db"), "Select formula database",
                                       choices = c(hmdb_demo = "hmdb",
                                                   norman_demo = "norman",
                                                   custom = "custom",
                                                   KEGG_Pathways = "kegg_p"),
                                       selected = "hmdb")),
                    column(8, conditionalPanel(
                        condition = "input.db == 'custom'", ns = ns,
                        {
                        tagList(
                        shinyFilesButton(ns("loadselector"),
                                         "Select custom database",
                                         "Select", multiple = FALSE,
                                         style = "margin-bottom: 10px;"),
                        verbatimTextOutput(ns("loadpath"), placeholder = TRUE))
                        }), conditionalPanel(
                        condition = "input.db == 'kegg_p'", ns = ns,
                        {
                          tagList(
                          pickerInput(ns("sel_organism"), choices = "", selected = "",
                                      options = list(`live-search` = TRUE)),
                          pickerInput(ns("sel_pathway"), choices = "", selected = "",
                                      options = list(`live-search` = TRUE,
                                                     `actions-box` = TRUE),
                                      multiple = TRUE)
                          )
                        })
                    )
                ),
        fluidRow(column(4,selectInput(ns("adcharge"), "Maximum adduct charge",
                                        choices = c(1,2,3), selected = 1)),
                column(4,selectInput(ns("admult"),
                                     "Maximum adduct multiplicity",
                                     choices = c(1,2,3), selected = 1)),
                column(4, pickerInput(
                    inputId = ns("ads"),
                    label = "Select the adducts to search",
                    choices = NULL,
                    options = list(`actions-box` = TRUE),
                    multiple = TRUE))),
        hr(),
        uiOutput(ns("startPL_box")),
        width = 13))
  )
}
PLServer <- function(id, struct){
    moduleServer(
        id,
        function(input, output, session){
            ns <- session$ns
            roots <- getVolumes()
            stat <- reactiveValues(stat = "Awaiting orders")
            #Update shown input values when template changes
            observeEvent(input$instrumentalmode,{
                df <- read.csv2(system.file("extdata",
                                            "InstrumentTemplates.csv",
                                            package = "RHermes"),
                                            stringsAsFactors = FALSE)
                selected <- df[df$name == input$instrumentalmode, ]
                updateNumericInput(session, "ppm", value = selected$ppm)
                updateNumericInput(session, "minMz", value = selected$minmz)
                updateNumericInput(session, "maxMz", value = selected$maxmz)
                updateNumericInput(session, "noise", value = selected$nthr)
                updateNumericInput(session, "resolution", value = selected$res)
                updateSelectInput(session, "ion",
                                  selected = as.character(selected$ion))
                updateSelectInput(session, "instrument",
                                  selected = as.character(selected$instr))
            }, ignoreInit = TRUE)

        observeEvent({
            input$adcharge
            input$admult
            input$ion
        },{
            adtable <- RHermes:::adductTables(as.numeric(input$adcharge),
                                    as.numeric(input$admult))
            if(input$ion == "-"){
                adductlist <- adtable[[1]]$adduct
            } else {
                adductlist <- adtable[[2]]$adduct
            }
            updatePickerInput(session, "ads", choices = adductlist,
                                selected = adductlist)

        }, ignoreNULL = TRUE)

        #mzML file selection logic
        shinyFileChoose(
            input,
            'files',
            roots = getVolumes(),
            filetypes = c("mzml", "mzxml")
        )
        #Database selection logic
        shinyFileChoose(
            input,
            'loadselector',
            roots = getVolumes(),
            filetypes = c("csv", "xlsx", "xls")
        )

        loadpath <- reactive(as.character(
            parseFilePaths(getVolumes(), input$loadselector)$datapath)
        )
        output$loadpath <- renderText(loadpath())
        output$selecteddir <- renderTable({parseFilePaths(roots, input$files)})



        observeEvent({
            input$files
            input$loadexample
        }, {
          if(nrow(parseFilePaths(roots, input$files)) == 0 & !input$loadexample){
            output$startPL_box <- renderUI({})
          } else {
            output$startPL_box <- renderUI({
              tags$div(actionButton(ns("startPL"), "Start PeakList generation",
                                    style = "background-color: #4d4263; color: #F0F0F0"),
                       style = "text-align: center;")})
          }
        }, ignoreNULL = TRUE, ignoreInit = TRUE
      )
        #KEGG pathway data loading
        orgnames <- read.csv(system.file("extdata", "keggOrganismNames.csv",
                                         package = "RHermes"))
        updatePickerInput(session, "sel_organism",
                          choices = paste(as.character(orgnames[,1]),
                                          as.character(orgnames[,2]),
                                          as.character(orgnames[,3]),
                                          sep = "-"),
                          selected = paste(as.character(orgnames[,1]),
                                           as.character(orgnames[,2]),
                                           as.character(orgnames[,3]),
                                           sep = "-")[1]
                          )

        observeEvent(input$sel_organism, {
          genome <- strsplit(x = input$sel_organism, split = "-")[[1]]
          if(length(genome) != 0){
            tryCatch({
            path <- KEGGREST::keggList("pathway", genome[2])
            updatePickerInput(session, "sel_pathway",
                              choices = paste(names(path), path, sep="-"),
                              selected = paste(names(path), path, sep="-")[1])
            }, error = function(e){message("No internet connection, cannot query KEGG")})

          }
        }, ignoreInit = TRUE, ignoreNULL = TRUE)

        #Output allocation
        toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)

        #Start data processing to PL
        observeEvent(input$startPL, {
            par <- ExpParam(ppm = input$ppm, res = input$resolution,
                            nthr = input$noise, minmz = input$minMz,
                            maxmz = input$maxMz,ion = input$ion,
                            instr = input$instrument)
            struct$dataset <- setExpParam(struct$dataset, par)

            keggpath <- lapply(strsplit(input$sel_pathway, "-"), function(x){
              x[[1]]
            })
            struct$dataset <- setDB(struct$dataset, db = input$db,
                                    adcharge = as.numeric(input$adcharge),
                                    admult = as.numeric(input$admult),
                                    filename = loadpath(), keggpath = keggpath)
            sel <- struct$dataset@metadata@ExpParam@adlist$adduct %in% input$ads
            ads <- struct$dataset@metadata@ExpParam@adlist[sel,]
            struct$dataset@metadata@ExpParam@adlist <- ads
            pth <- parseFilePaths(roots, input$files)$datapath
            if(input$loadexample){pth <- system.file("extdata", "MS1TestData.mzML",
                                                              package = "RHermes")}
            struct$dataset <- processMS1(struct$dataset, pth,
                                       labelled = as.logical(input$labelled))
            sendSweetAlert(
                session = session,
                title = "All files have been processed",
                text = "Now you can visualize the results and/or start SOI
                generation",
                type = "success"
            )
            toReturn$dataset <- struct$dataset
            toReturn$trigger <- toReturn$trigger + 1
        })


        #Help Modals
        observeEvent(input$param_help, {
          showModal(ui = modalDialog(
              tagList(
                tags$p("Here you have all parameters explained: "),
                tags$ul(
                  tags$li(tags$b("Presets"),": we have designed some ready-to-use presets so that
                          you don't have to mess around with the parameters. Feel free to tune them
                          to your liking."),
                  tags$li(tags$b("ppm"),": instrumental mass error in parts per million.
                        In case of doubt, set a value a bit (1-2 ppm) higher
                        than usual to avoid losing peaks."),
                  tags$li(tags$b("m/z",), ": set both the minimum and maximum mz range you are
                                 interested with. This interval is used to filter the formulas
                                 in your selected database."),
                  tags$li(tags$b("Instrument noise"),": the minimum intensity of the signals to detect.
                          This parameter is especially useful for Q-TOF instruments which tend to have
                          a lot of noisy signals (Eg: \"grass\")."),
                  tags$li(tags$b("Resolution"),": the instrument resolution power at m/z = 200. It is
                          used to calculate distinguishable isotopic peaks. In Q-TOF instruments it is
                          assumed to remain constant along the mzs while in Orbitrap it grows proportionally
                          with the square root of the mz."),
                  tags$li(tags$b("Ionization mode"), ": either positive or negative (but not both).
                          We currently don't support switched polarity experiments. In case you want
                          to process such data, you could filter the positive and negative scans of each file
                          (eg. with MSConvert) and process them in separate RHermesExp experiments
                          (that is, running this app twice)."),
                  tags$li(tags$b("Instrument type"), ": either Orbitrap or Q-TOF. The main difference is in the
                          distinguishable isotopic peak calculation (see Resolution)."),
                  tags$li(tags$b("Labelled data"), ": activating this will make RHermes significantly slower,
                          since it checks for ALL 13C carbon signals, meaning that, for a compound like C6H12O6 it would try
                          to find all M0, M1, ..., M6 peaks. Don't process many files at once (not that you should
                          in the regular processing) and, if possible, avoid high carbon metabolites in your database.")
                ),
                tags$p("To use RHermes you need to provide a suitable
                       formula database. You can try the \"hmdb_demo\" or the \"norman_demo\" which just have a handful of
                       compounds to get a taste of the RHermes processing. For real usage, you must either provide a .csv or
                       .xls/.xlsx that list all formulas to search in the data (select custom option) or use our KEGG Pathways
                       interface to select the pathways whose compounds you're interested in."),
                tags$p("The csv/xls/xlsx file needs to have 2 columns that match the following names: Name and MolecularFormula.
                       This is crucial for the correct usage of the program. Other columns are fine and, although they will be saved in
                       your RHermesExp object, do not play any role in the data processing.")
              ),
          easyClose = TRUE, title = HTML('<h4 class = "help-header">Peaklist Parameters and Databases</h4>')),
          session = session)
        })
        return(toReturn)
    }
  )
}
