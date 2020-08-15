#' @export
SOI_UI <- function(id){
  ns <- NS(id)
  tagList(
    useSweetAlert(),
    tabsetPanel(
      tabPanel("SOI generation",
               verticalLayout(
                 sidebarPanel(
                   div(
                     tags$span(h2("SOI Detection Parameters"), style = "display: inline-block; width: 90%"),
                     tags$span(actionButton(class = "help", label ="", inputId = ns("param_help"),
                                            icon = icon("info-circle"))),
                     style = "text-align: center"
                   ),
                   hr(),
                   fluidRow(
                     column(4,
                            selectInput(inputId = ns("paramcsv"),
                                        label = "Input the filter preset you want to use:",
                                        choices = "", selected = "",
                                        multiple = FALSE),
                     ),
                     column(6,offset = 2,
                            tableOutput(outputId = ns("paramtable"))
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(6,
                            numericInput(ns("minint"),
                                         "Minimum SOI Intensity", 1000),
                            numericInput(ns("maxlen"),
                                         "Maximum SOI RT length", 30),
                     ),
                     column(6,
                            radioButtons(ns("blanksub"),"Blank Substraction",
                                         c("T", "F"), selected = "T")
                     )
                   ),
                   selectInput(ns("sampleID"),
                               "Which of the files (index) corresponds to the sample:",
                               choices = "", selected = ""),
                   conditionalPanel("input.blanksub == 'T'", ns = ns,{
                     selectInput(ns("blankID"),
                     "If using Blank Substraction, which of the files (index)
                          corresponds to the blank:", choices = "",
                     selected = "")
                   }),

                   uiOutput(ns("soibuttons")),
                   tags$p("Here will appear the selected entries to process into SOI:"),
                   verbatimTextOutput(ns("description"), TRUE),
                   width = 13
                 )
               )),
      tabPanel("SOI filtering - Cleanup",
               verticalLayout(
                 sidebarPanel(
                   div(
                     tags$span(h2("SOI Cleanup Parameters"), style = "display: inline-block; width: 90%"),
                     tags$span(actionButton(class = "help", label ="", inputId = ns("param_help_clean"),
                                            icon = icon("info-circle"))),
                     style = "text-align: center"
                   ),
                   hr(),
                   uiOutput(ns("SOIcleanup")),

                  width = 12)
                )
              )
    ),

  )
}

#' @export
SOIServer <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){
      #SOI parameter auto-update
      presetdf <- read.csv2(system.file("extdata", "SOITemplates.csv",
                                        package = "RHermes"),
                            stringsAsFactors = FALSE)
      updateSelectInput(session, "paramcsv", choices = unique(presetdf$name),
                        selected = unique(presetdf$name)[1])
      observeEvent(input$paramcsv,{
        selected <- presetdf[presetdf$name == input$paramcsv, ]
        output$paramtable <- renderTable(as.matrix(selected[,c("binwidth",
                                                               "minscan",
                                                               "shift")]))
        updateNumericInput(session, "minint", value = selected$minint[1])
        updateNumericInput(session, "maxlen", value = selected$maxlen[1])
      }, ignoreInit = TRUE)

      #Show crucial buttons only when we know there are PL to process
      ns <- session$ns #To refer to the correct namespace
      observeEvent(struct$hasPL, {
        output$soibuttons <- renderUI({
          tagList(
            fluidRow(
              hr(style = "border-top : 1px dashed #BEBBEB"),
              span(
                actionButton(inputId = ns("addfile"), "Add selected config"),
                actionButton(inputId = ns("remfile"), "Remove selected config")
              ),
              uiOutput(ns("startbut"), inline = TRUE,
                       style = "margin-left: 5%"),
              br(),
              hr(style = "border-top : 1px dashed #BEBBEB"),
              style = "margin: auto; text-align:center;"
            )
          )
        })
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      #Modify file names as they are processed into PLs
      fnames <- reactive(unname(struct$dataset@metadata@filenames))
      observeEvent(fnames(), {
        updateSelectInput(session, "sampleID", choices = fnames(),
                          selected = fnames()[1])
        updateSelectInput(session, "blankID", choices = fnames(),
                          selected = fnames()[1])
      })



      #Configuration reactivity
      configdf <- reactiveValues(DF = data.frame(file = character(0),
                                                 blank = character(0),
                                                 template = character(0),
                                                 minint = numeric(0),
                                                 maxlen = numeric(0)))
      observeEvent(input$addfile,{
        if(as.logical(input$blanksub)){
            curdf <- isolate(configdf$DF)
            configdf$DF <- rbind(curdf, data.frame(file = input$sampleID,
                                                   blank = input$blankID,
                                                   template = input$paramcsv,
                                                   minint = input$minint,
                                                   maxlen = input$maxlen,
                                                   stringsAsFactors = FALSE))
        } else {
            curdf <- isolate(configdf$DF)
            configdf$DF <- rbind(curdf, data.frame(file = input$sampleID,
                                                   blank = "NO",
                                                   template = input$paramcsv,
                                                   minint = input$minint,
                                                   maxlen = input$maxlen,
                                                   stringsAsFactors = FALSE))
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent(input$remfile,{
        if(nrow(configdf$DF) != 0){
            curdf <- configdf$DF
            configdf$DF <- curdf[-nrow(curdf), ]
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      output$description <- renderPrint(configdf$DF)

      observeEvent(nrow(configdf$DF),{
        if(nrow(configdf$DF) == 0){
          output$startbut <- renderUI(expr = {})
        } else {
          output$startbut <- renderUI(
            actionButton(ns("startSOIdetect"), "Start SOI detection",
                       icon = icon("barcode"),
                       style = "background-color: #4d4263; color: #F0F0F0"))
        }
      })

      toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
      observeEvent(input$startSOIdetect, {
        #Parse the user input configdf to get the right input for SOIfinder
        fileidx <- c()
        blankidx <- c()
        paramlist <- c()
        for(i in seq_len(nrow(configdf$DF))){
          currow <- configdf$DF[i,]
          fileidx <- c(fileidx, which(fnames() == currow$file))
          blankidx <- c(blankidx,
                        ifelse(currow$blank == "NO", 0,
                               which(fnames() == currow$blank)))
          paramlist <- c(paramlist, getSOIpar(currow$template))
          paramlist[[i]]@minint <- currow$minint
          paramlist[[i]]@maxlen <- currow$maxlen
        }

        #Process and output
        toReturn$dataset <- SOIfinder(struct$dataset, paramlist, fileidx,
                                      blankidx)
        toReturn$trigger <- toReturn$trigger + 1

        sendSweetAlert(
          session = session,
          title = "All SOI lists have been generated",
          text = "Now you can visualize the results and/or generate an IL",
          type = "success"
        )

      }, ignoreNULL = TRUE, ignoreInit = TRUE)


      observeEvent({
        struct$hasSOI
        struct$dataset@data@SOI
        },{

        output$SOIcleanup <- renderUI({
          tagList(
            fluidRow(
                column(6, radioButtons(inputId = ns("soiChoices"),
                           label =  "Select which SOI list to clean: ",
                           choices = seq_along(struct$dataset@data@SOI),
                           selected = 1)),
                column(6,
                      tags$b("Enable isotopic fidelity filtering"),
                      div(switchInput(
                        inputId = ns("isofiltering"),
                        value = TRUE,
                        onLabel = "Yes",
                        offLabel = "No"),
                      style = "margin-top: 5px"
                      ))
            ),
            numericInput(ns("numericval"), label = "Minimum SOI intensity", value = 0,
                         min = 0, max = 1e5, step = 1000),
            fluidRow(
              actionButton(inputId = ns("soiClean"), label = "Start SOI Cleaning",
                           style = "text-align: center; background-color: #4d4263; color: #F0F0F0"),
              style = "margin: auto; text-align:center;"
            )
          )})

      }, ignoreNULL = TRUE,ignoreInit = TRUE, priority = 100)

      observeEvent(input$minfilter,{

        updateNumericInput(session, "numericval", value = input$minfilter)

      })

      observeEvent(input$soiClean,{


        SOIcleaner(struct$dataset, as.numeric(input$soiChoices),
                   as.numeric(input$minfilter), input$isofiltering)
        sendSweetAlert(session = session, title = "Finished",
                       text = paste("The SOI list", input$soiChoices,
                                    "has been cleaned"),
                       type = "success")

      })

      #Help Modals
      observeEvent(input$param_help, {
        showModal(ui = modalDialog(
          tagList(
            tags$p("Here you have all parameters explained: "),
            tags$ul(
              tags$li(p(tags$b("Filter Presets"),": These are different kinds of filters you can run over the data points
                      in order to find SOIs (Scans of Interest). The main difference between them is their robustness
                      aganist peak splitting. Double and triple generally give better results but take a bit longer."),
                      p("All the \"-x\" filters have increased RT window sizes and are designed with long chromatography
                      experiments in mind (ie. >20-30 min), where compound peaks are generally longer and small windows
                      would take longer."),
                      p("Our advice is to try out \"simple\" or \"double\" and adjust if necessary (but we have good
                        results with these).")),
              tags$li(tags$b("Minimum SOI intensity"), ": The minimum intensity of scans for the algorithm to consider.
                      Is best to leave it as default, but if you're looking only for intense signals, go on and increase it to
                      5-10-20K and the algorithm running time will be greatly reduced. Remember that you can always filter
                      SOIs in the optional SOI cleanup stage."),
              tags$li(tags$b("Maximum SOI RT"), ": The algorithm will split long SOIs (larger than the selected value),
                      into smaller ones. Longer values take shorter (since there aren't as many SOI entries), but blank
                      substraction problems can arise and it's not really worth it. Better left as default."),
              tags$li(tags$b("Blank Substraction"), ": Whether to perform Blank Substraction (Optional but recommended).",
                      tags$b("Important:"), "The blank substraction step requires setting up Keras and Tensorflow first.
                      If you haven't done it please check the package Readme file to see how to do it.")
            ),
            tags$p("To generate the SOI list (or lists) you have to set up the desired parameters and click \"Add selected config\".
            When you do so, you will see that the input parameters appear in the box below the button. You can either add more configs
            (for instance, try different sample/blank pairs for comparison, try filter presets, etc.) or remove them. Once you're done,
            click \"Start SOI detection\" to start the processing. Check the R terminal for details on the process. You will recieve an
            alert on the app once the processing is finished.")
          ),
          easyClose = TRUE, title = HTML('<h4 class = "help-header">SOI list generation</h4>')),
          session = session)
      })

      observeEvent(input$param_help_clean, {
        showModal(ui = modalDialog(
          tagList(
            tags$p("Here you have all parameters explained: "),
            tags$ul(
              tags$li(tags$b("Isotopic Fidelity Filtering"),": Whether to filter intense SOIs based on their
                      isotopic ratios. The idea is to compare the isotopic intensities of the SOI to their
                      theoretical abundances. We've experimentally checked that all identified compounds
                      pass this filter. The only caveat is when the ppm parameter is set too tight and certain
                      isotopologue peaks are outside the error margin, which leads to some false negatives.
                      It is recommended to manually check the PL plots of some known compounds if you want to
                      make sure this isn't your case."),
              tags$li(tags$b("Minimum SOI intensity"), ": Intensity threshold required at the most intense point of the SOI
                      (analogous to maxo in XCMS). It filters all SOIs whose intensity is lower than it. Set it to 0 if
                      you don't want to filter anything. Keep in mind that more stringent filters will result in fewer SOIs
                      but a \"richer\" inclusion list that requires less injections and has most intense metabolites.")
            ),
            tags$p("Choose the SOI list to clean, select your desired parameters and click \"Start SOI cleaning\" to begin.
                   The cleaning process is also required if you want to prioritize certain adduct annotations over others when
                   generating the inclusion list.")
          ),
          easyClose = TRUE, title = HTML('<h4 class = "help-header">SOI list refining</h4>')),
          session = session)
      })



      return(toReturn)
    }
  )
}
