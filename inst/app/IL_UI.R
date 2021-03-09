IL_UI <- function(id){
  ns <- NS(id)
  tagList(
    useSweetAlert(),
    verticalLayout(
      sidebarPanel(
        div(
          tags$span(h2("Inclusion List generation"), style = "display: inline-block; width: 90%"),
          tags$span(actionButton(class = "help", label ="", inputId = ns("param_help"),
                                 icon = icon("info-circle"))),
          style = "text-align: center"
        ),
        hr(),
        uiOutput(ns("whenSOI")),
        width = 13),
      sidebarPanel(
        h2("Export your inclusion list", style = "text-align: center;"),
        hr(),
        fluidRow(
          column(6,
            radioButtons(ns("ILexportmode"), "Export mode:",
                         choices = c("As a single csv file",
                                     "In multiple csv files"),
                         selected = "As a single csv file"),
            shinySaveButton(ns("singleexport"),
                            "Select your base filename:",
                            title = "Select your base filename:"),

          offset = 1),
          column(3,
                 uiOutput(ns("whenIL")),
                 offset = 0)
        ),
        verbatimTextOutput(ns("completedtext")),
        width = 13
      )
    )
  )
}

ILServer <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){
      ns <- session$ns #To refer to the correct namespace
      adlist <- reactive(struct$dataset@metadata@ExpParam@adlist)
      observeEvent(struct$hasSOI, {
        if(struct$hasSOI){
          output$whenSOI <- renderUI({tagList(
            numericInput(ns("filtermz"), label = "Select the precursor ion filter mz width (+-x)", value = 0.5, min = 0.1,
                         max = 20, step = 0.1),
            radioButtons(ns("soiid"), "Select which SOI (id) you want to process into an IL:",
                         choices = seq_along(struct$dataset@data@SOI)),
            radioButtons(ns("mode"), "Select the inclusion list generation mode",
                         choices = c("All adducts", "Only a subset of adducts", "Prioritize some adducts over the rest"),
                         selected = "All adducts"),
            conditionalPanel("input.mode != 'All adducts'", ns = ns, {
              selectInput(ns("adlist"), label = "Select which adduct/s (if any) you want to prioritize/select:",
                        selectize = TRUE, multiple = TRUE, choices = adlist()$adduct)
              }),
            actionButton(ns("startgenIL"), label = "Start IL generation")
          )})
        } else {
          output$whenSOI <- renderUI({tagList()}) #Remove option if PL not available
        }
      })
      numsoi <- reactive(seq_along(struct$dataset@data@SOI))
      observeEvent(numsoi(),
                   if(struct$hasSOI){
                     updateRadioButtons(session = session, inputId = ns("soiid"),
                                                choices = numsoi())
                     },
                   ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent(struct$hasIL,{
        if(struct$hasIL){
          output$whenIL <- renderUI({
            tagList(
              radioButtons(ns("ILidx"), "Select which IL you want to export: ", choices = numIL(), selected = numIL()[1]),
              numericInput(ns("maxover"), "Maximum number of entries to monitor at the same time: ",
                           min = 1, max = 50, value = 5, step = 1,),
              tags$div(actionButton(ns("exportIL"), "Export the selected IL",
                                    style = "background-color: #4d4263; color: #F0F0F0"),
                       style = "text-align: center;")
            )
          })
        } else {
          output$whenIL <- renderUI({tagList()})
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)
      numIL <- reactive(seq_along(struct$dataset@data@MS2Exp))
      observeEvent(numIL(),
                   if(struct$hasIL){
                     updateRadioButtons(session = session, inputId = ns("ILidx"),
                                        choices = numIL())
                   },
                   ignoreNULL = TRUE, ignoreInit = TRUE)

      shinyFileSave(
        input,
        'singleexport',
        roots = getVolumes()
      )

      observeEvent(input$exportIL,{
        sepFiles <- switch(input$ILexportmode,
                           "In multiple csv files" = TRUE,
                           "As a single csv file" = FALSE
                           )
        fname <- as.character(parseSavePath(getVolumes(), input$singleexport)$datapath)
        message("Now exporting the IL, it may take a bit")
        exportIL(struct$dataset, as.numeric(input$ILidx), file = fname ,
                 maxOver = input$maxover, sepFiles = sepFiles)
        output$completedtext <- renderText(paste("The selected IL has been exported to:", fname))
      }, ignoreInit = TRUE, ignoreNULL = TRUE)

      toReturn <- reactiveValues(dataset = RHermesExp(), trigger = 0)
      observeEvent(input$startgenIL,{
        filtermz <- input$filtermz
        mode <- switch (input$mode,
                        "All adducts" = "full",
                        "Only a subset of adducts" = "only",
                        "Prioritize some adducts over the rest" = "yes"
                        )
        ad <- ifelse(is.null(input$adlist), character(0), input$adlist)
        par <- ILParam(filtermz = filtermz, filterrt = 10, rtmargin = 5,
                       priorization = mode, ad = ad)
        toReturn$dataset <- generateIL(struct$dataset, id = as.numeric(input$soiid), par = par)
        toReturn$trigger <- toReturn$trigger + 1
        sendSweetAlert(
          session = session,
          title = "IL generated",
          text = "Now you can export it below",
          type = "success"
        )
      })

      #Help Modals
      observeEvent(input$param_help, {
        showModal(ui = modalDialog(
          tagList(
            tags$p("Here you have all parameters explained: "),
            tags$ul(
              tags$li(tags$b("mz window"),": the precursor ion isolation window. It is dependent on your instrument
                      capabilities and the same value must be used when performing the continuous DDA experiments.
                      Higher mz windows allow for shorter inclusion lists, but are harder to deconvolute. Low mz windows
                      are very selective (easier deconvolution) but are associated with sensitivity loss. The value set
                      here is the semiwidth, which is the maximum error from the targeted mz"),
              tags$li(tags$b("IL mode"),": there are 3 differents modes to choose from. You can include all SOIs independently
                      of their adduct annotation, you can select only SOIs that have a subset of adduct annotations or you
                      can prioritize some annotations over others (which means that, if a compound elutes as M+H, M+Na and M+NH4
                      and you choose M+H, only that SOI will be included in the IL)."),
              tags$li(tags$b("Overlap limit"), ": the maximum number of IL entries the instrument can monitor at any given time.
                      Higher limits will result in less injections but lower scan rate for each IL entry, which can negatively
                      affect the MS2 analysis. It really depends on your instrument, but values between 5-10 have been tested
                      with good results. A duty cycle of less that 1 second is desired to catch all eluting compounds but,
                      if the chromatography is long enough, one can raise the overlap limit.")
            ),
            tags$p("To generate the IL you have to choose which SOI list to use, set up the parameters and click \"Generate Inclusion
            List\". Once finished, you can export it below as a single csv or multiple csv (one for each planned injection)."),
            tags$p("Perform the continuous DDA importing the corresponding csv experiment and setting up your instrument
                   accordingly. When you have the MS2 files, proceed with the next step. In the meantime, consider saving your
                   progress in the Loading/Saving tab.")
          ),
          easyClose = TRUE, title = HTML('<h4 class = "help-header">IL Parameters</h4>')),
          session = session)
      })

      return(toReturn)
    }
  )
}
