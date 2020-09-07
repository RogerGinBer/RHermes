SOIPlotUI <- function(id){
  ns <- NS(id)
  tagList(
    tags$br(),
    radioGroupButtons(
      inputId = ns("selectplot"),
      label = "Select a graph :",
      choices = c(`<i class="fas fa-mountain"></i> SOI plots` = "soiplot",
                  `<i class='fa fa-bar-chart'></i> Isotopic fidelity` = "fidelityplot"),
      justified = TRUE
    ),
    conditionalPanel("input.selectplot == 'soiplot'", ns = ns,
       fluidRow(
         column(4, offset = 1,
                radioButtons(ns("SOIplotmode"), label = "Select mode:",
                             c("By compound name", "By formula"), selected = "By compound name"),
                conditionalPanel("input.SOIcompselect != 'No file detected'", ns = ns,
                                   radioButtons(ns("soifiles"), "Select which SOI list to use: ", choices = "", 1, inline = TRUE)),
                tags$b("Blank Substraction file:"),
                textOutput(ns("blanktext"))
         ),

         column(6, offset = 0,
                fluidRow(
                  column(8,
                         conditionalPanel(condition = "input.SOIplotmode == 'By compound name'", ns = ns,
                                          selectizeInput(inputId = ns("SOIcompselect"), label ="Select compound",
                                                         choices = "No file detected", selected = NULL, multiple = FALSE),
                         ),
                         conditionalPanel(condition = "input.SOIplotmode == 'By formula'", ns = ns,
                                          selectizeInput(inputId = ns("SOIformselect"), label = "Select formula",
                                                         choices = "No file detected", selected = NULL, multiple = FALSE),
                         )),

                  column(4, tags$b("Dynamic axis"), switchInput(
                    inputId = ns("dynamicaxis"),
                    onStatus = "success",
                    offStatus = "danger",
                    value = TRUE, size = "small"
                  ))

                ),
                conditionalPanel("input.SOIcompselect != 'No file detected'", ns = ns,
                                 sliderInput(ns("RTinterval"), label = "You can select an RT interval here:",
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
                                 tags$head(tags$style("#SOIPlotUI-othercomp{overflow-y:scroll; max-height: 250px; background: ghostwhite;}")))
         )
       ),
     plotlyOutput(outputId = ns("SoiPlot"), height = "800px")
    ),
    conditionalPanel("input.selectplot == 'fidelityplot'", ns = ns,
                     selectizeInput(ns("choicesfidelity"), choices = NULL, selected = NULL,
                                    label = "Select SOI to check", width = "600px"),
                     verbatimTextOutput(ns("valoration")),
                     plotlyOutput(ns("fplot"), height = "600px")
                     )
  )
}
SOIPlotServer <- function(id, struct){
  moduleServer(
    id,
    function(input, output, session){

    ns <- session$ns
    observeEvent({
      struct$hasSOI
      struct$dataset@data@SOI
      },{
     if(struct$hasSOI){


       soiNames <- lapply(struct$dataset@data@SOI, function(x){
          parsed <-  strsplit(x@filename, "/")[[1]]
          return(parsed[length(parsed)])
       })

       soiIndex <- seq_along(struct$dataset@data@SOI)
       names(soiIndex) <- soiNames

       updateSelectizeInput(session, "SOIcompselect", choices = unique(struct$dataset@metadata@ExpParam@DB$Name),
                            selected = unique(struct$dataset@metadata@ExpParam@DB$Name[1]),
                              server = TRUE)
       updateSelectizeInput(session, "SOIformselect", choices = unique(struct$dataset@metadata@ExpParam@DB$MolecularFormula),
                            selected = unique(struct$dataset@metadata@ExpParam@DB$MolecularFormula[1]),
                            server = TRUE)
       updateRadioButtons(session, "soifiles", choices = soiIndex)
       updateSliderInput(session, "RTinterval", max = ceiling(max(struct$dataset@data@SOI[[1]]@PlotDF$rt)),
                         value = c(0, ceiling(max(struct$dataset@data@SOI[[1]]@PlotDF$rt))))
       updatePickerInput(session, "ads", choices = struct$dataset@metadata@ExpParam@adlist$adduct,
                         selected = struct$dataset@metadata@ExpParam@adlist$adduct)
     }else{
       updateSelectizeInput(session, "SOIcompselect", choices = "No file detected")
       updateSelectizeInput(session, "SOIformselect", choices = "No file detected")
       updateRadioButtons(session, "soifiles", choices = "")
     }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 100)

    observeEvent({
      input$SOIcompselect
      input$SOIformselect
      input$soifiles
      input$SOIplotmode
      input$RTinterval
      input$dynamicaxis
      input$ads
    },{
      if(struct$hasSOI){
        if(input$SOIplotmode == "By compound name"){
          if(input$SOIcompselect == ""){return()}
          f <- with(struct$dataset@metadata@ExpParam@DB,{MolecularFormula[Name == input$SOIcompselect]})[1]
        } else {
          if(input$SOIformselect == ""){return()}
          f <- input$SOIformselect[1]
        }
        blankNames <- lapply(struct$dataset@data@SOI, function(x){
          if(x@SoiParam@blanksub){
              parsed <- strsplit(x@SoiParam@blankname, "/")[[1]]
              return(parsed[length(parsed)])
          } else {
              return("No blank substraction")
          }
        })
        blankID <- lapply(struct$dataset@data@SOI, function(x){
          if(x@SoiParam@blanksub){
            id <- which(struct$dataset@metadata@filenames ==
                          x@SoiParam@blankname)
            return(id)
          } else {
            return(NA)
          }
        })
        output$blanktext <- renderText(blankNames[[as.numeric(input$soifiles)]])
        other <- with(struct$dataset@metadata@ExpParam@DB, {Name[MolecularFormula == f]})
        output$othercomp <- renderText(paste(other, collapse = "\n"))
        output$SoiPlot <- renderPlotly(SoiPlot(struct = struct$dataset,
                                                              id = as.numeric(input$soifiles),
                                                              formula = f, rtrange = as.numeric(input$RTinterval),
                                                              dynamicaxis = as.logical(input$dynamicaxis),
                                                              ads = as.character(input$ads),
                                                              blankid = blankID[[as.numeric(input$soifiles)]]))
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    observeEvent({
      input$SOIcompselect
      input$SOIformselect
      input$soifiles
      input$SOIplotmode
      input$RTinterval
      input$ads
    },{
      if(struct$hasSOI){
        if(input$SOIplotmode == "By compound name"){
          if(input$SOIcompselect == ""){return()}
          formula <- with(struct$dataset@metadata@ExpParam@DB,{MolecularFormula[Name == input$SOIcompselect]})
        } else {
          if(input$SOIformselect == ""){return()}
          formula <- input$SOIformselect
        }
        # other <- with(struct$dataset@metadata@ExpParam@DB, {Name[MolecularFormula == f]})
        # output$othercomp <- renderText(paste(other, collapse = "\n"))
        #

        equivalences <- struct$dataset@metadata@ExpParam@ionF[[2]]
        equivalences <- equivalences[equivalences$f == formula,]
        equivalences <- equivalences[equivalences$an %in% input$ads]
        ions <- equivalences$ion

        SOIlist <- struct$dataset@data@SOI[[as.numeric(input$soifiles)]]@SoiList
        ids <- which(SOIlist$formula %in% ions)
        fa_names <- vapply(SOIlist$formula[ids], function(x){
          paste(unlist(equivalences[equivalences$ion == x, c("f", "an")]), collapse = "~")
        }, character(1)
        )
        rows <- paste(ids, fa_names, rep(" Intensity:", length(ids)),round(SOIlist$MaxInt[ids]))
        updateSelectizeInput(session, "choicesfidelity", choices = rows)


      }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)



    observeEvent({input$choicesfidelity},{
      if(!is.null(input$choicesfidelity) & input$choicesfidelity != ""){
        selected <- strsplit(input$choicesfidelity, " ")[[1]][[1]] %>% as.numeric(.)
        isoresults <- IsoFidelity(struct$dataset, as.numeric(input$soifiles), selected)
        output$fplot <- renderPlotly(isoresults[[1]])

        heteroatoms <- isoresults[[5]]
        heteroatoms <- heteroatoms[!is.null(heteroatoms)]
        if(length(heteroatoms) != 0){
          het_text <- paste0(names(heteroatoms), heteroatoms, collapse = "\n")
        } else {
          het_text <- ""
        }
        output$valoration <- renderText(paste0("Calculated number of carbon atoms: ", isoresults[[2]], "\n",
                                              "Veredict: ", isoresults[[3]], "\n",
                                            "Isotopic cosine score: ", isoresults[[4]], "\n",
                                            het_text))
      }
    })

    return()
  }
  )
}


