MS2PlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$br(),

    radioGroupButtons(inputId = ns("selectclass"), label = "Select a graph :",
                      choices = c(`<i class='fa fa-bar-chart'></i> Identifications` = "ident",
                                  `<i class="fas fa-database"></i> Raw MSMS data` = "raw"),
                      justified = TRUE),
    radioGroupButtons(ns("id"), "Select MS2Exp:", choices = "",
                      selected = ""),
    br(),

    conditionalPanel("input.selectclass == 'ident'", ns = ns,
                     radioGroupButtons(ns("selectplot"),
                                       choices = c(`Superspectra Plot` = "sup",
                                                   `Superspectra Comparison` = "comp",
                                                   `Compound matches` = "hits"), justified = TRUE
                     ),
                     conditionalPanel("input.selectplot == 'sup'", ns = ns,
                                      fluidRow(column(6, selectizeInput(ns("selectss"),
                                                                        label = "Select superspectrum to plot",
                                                                        choices = NULL, selected = NULL))),
                                      plotlyOutput(ns("ssplot"))),
                     conditionalPanel("input.selectplot == 'hits'", ns = ns,
                                      fluidRow(
                                        column(6,
                                               selectizeInput(ns("selectident"),
                                                              label = "Select the superspectra entry to check:",
                                                              choices = NULL, selected = NULL)),
                                        column(6,
                                               pickerInput(ns("selecthits"),
                                                           label = "Select which hits to check against:",
                                                           choices = NULL, selected = NULL, multiple = TRUE,
                                                           options = list(
                                                             `actions-box` = TRUE))
                                        )
                                      ),
                                      plotlyOutput(ns("mirrorplot"))
                     ),

                     conditionalPanel("input.selectplot == 'comp'", ns = ns,
                                      fluidRow(
                                        column(6,
                                               selectizeInput(ns("selectident2"),
                                                              label = "Select the superspectra entry to check:",
                                                              choices = NULL, selected = NULL)),
                                        column(6,
                                               pickerInput(inputId = ns("otherss"),
                                                           label = "Select which superspectra/s to compare with:",
                                                           choices = NULL, selected = NULL, multiple = TRUE,
                                                           options = list(`max-options` = 10))
                                        )
                                      ),
                                      plotlyOutput(ns("versusplot"))
                     )
    ),
    conditionalPanel("input.selectclass == 'raw'", ns = ns,
                     fluidRow(
                       column(6,
                              checkboxGroupButtons(
                                inputId = ns("rawset"),
                                label = "Select what to plot:",
                                choices = c("Raw data plot", "Network", "Peak table"),
                                checkIcon = list(
                                  yes = tags$i(class = "fa fa-check-square",
                                               style = "color: #4D4263"),
                                  no = tags$i(class = "fa fa-square-o",
                                              style = "color: #4D4263")),
                                selected = c("Raw data plot", "Network", "Peak table")
                              )
                       )
                     ),
                     fluidRow(column(width = 6,
                                     selectizeInput(ns("selectILentry"),
                                                    label = "Select the IL entry to check:", choices = NULL,
                                                    selected = NULL
                                     )
                     ),
                     column(width = 6, tags$b("Select plot mode:"),
                            switchInput(ns("bymz"), onLabel = "By m/z",
                                        offLabel = "By group", value = TRUE,
                                        labelWidth = "AUTO"))),
                     fluidRow(
                       conditionalPanel("input.rawset.includes('Raw data plot')",
                                        ns = ns,
                                        conditionalPanel("input.bymz", ns = ns,
                                                         plotlyOutput(ns("rawMSMS_bymz"), height = "800")),
                                        conditionalPanel("!input.bymz", ns = ns,
                                                         plotlyOutput(ns("rawMSMS_bygroup"), height = "800")),
                       ),
                       conditionalPanel("input.rawset.includes('Network')",
                                        ns = ns, visNetworkOutput(ns("rawss"))),
                       conditionalPanel("input.rawset.includes('Peak table')",
                                        ns = ns,
                                        column(12, align="center",
                                               tableOutput(ns('rawpks'))))
                     )
    ),

  )

}

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
        updateRadioGroupButtons(session, "id", choices = whichMS2,
                                selected = whichMS2[1])

      } else {
        #Hide panels again
        updateRadioGroupButtons(session, "id", choices = NULL,
                                selected = NULL)
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 100)

    observeEvent({
      input$id
    },
    {
      ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
      sslist <- ms2@Ident[["MS2Features"]]
      with_hits <- which(vapply(sslist$results,
                                function(x){is.data.frame(x)},
                                logical(1)))
      original_IL <- unique(sslist$ILentry)
      updateSelectizeInput(session, "selectident",
                           choices = as.character(with_hits), server = TRUE,
                           selected = as.character(with_hits[1]),
                           options = list(maxOptions = 50000))
      updateSelectizeInput(session, "selectILentry",
                           choices = as.character(original_IL),
                           selected = as.character(original_IL[1]),
                           server = TRUE, options = list(maxOptions = 50000))
      updateSelectizeInput(session, "selectident2",
                           choices = as.character(seq_len(nrow((sslist)[1]))), server = TRUE,
                           selected = as.character(1),
                           options = list(maxOptions = 50000))
      updateSelectizeInput(session, "selectss",
                           choices = as.character(seq_len(nrow((sslist)[1]))), server = TRUE,
                           selected = as.character(1),
                           options = list(maxOptions = 50000))

    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)

    observeEvent({input$selectident2},{
      ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
      sslist <- ms2@Ident[["MS2Features"]]

      otherss <- seq_len(nrow(sslist))[-as.numeric(input$selectident2)]

      updatePickerInput(session, "otherss",
                        choices = as.character(otherss),
                        selected = as.character(otherss[1]))
    }, ignoreNULL = TRUE, ignoreInit = TRUE
    )

    #### Plots ####
    observeEvent({
      input$selectILentry
      input$id
    }, {
      if (input$selectILentry != "" & !is.na(input$selectILentry)) {
        tryCatch({
          rawplots <- plotRawMS2(struct$dataset,
                                 as.numeric(input$id), as.numeric(input$selectILentry))
          rtrange <- range(rawplots$p_bymz$data$rt)

          filtermz <- IL(struct$dataset, as.numeric(input$id))@ILParam@filtermz
          mzprec <- IL(struct$dataset, as.numeric(input$id))@IL$mass[as.numeric(input$selectILentry)]
          mzrange <- c(mzprec - filtermz, mzprec + filtermz)
          
          original_file <- SOI(struct$dataset,
              IL(struct$dataset, as.numeric(input$id))@SOInum)@filename
          original_PL <- which(vapply(struct$dataset@data@PL, function(x) {
              return(x@filename == original_file)
          }, logical(1)))[1]
          original_PL <- PL(struct$dataset, original_PL)
          
          annotated_scans <- filter(original_PL@peaklist,
                                    between(rt, rtrange[1], rtrange[2]) &
                                    between(mz, mzrange[1], mzrange[2]))
          raw_scans <- filter(original_PL@raw,
                              between(rt, rtrange[1], rtrange[2]) &
                              between(mz, mzrange[1], mzrange[2]))
          suppressWarnings({
              ms1plot <- ggplot() + 
                  geom_point(aes(x=rt, y=rtiv, mz = mz), data = raw_scans, color = "grey") + 
                  geom_point(aes(x=rt, y=rtiv, color = formv, mz = mz), data = annotated_scans) +
                  labs(color = "Ion. Formula") + 
                  ggtitle("Above: continuous MS2 scans<br>Below: MS1 data within isolation window")+
                  theme_minimal()
          })
          
          output$rawMSMS_bymz <- renderPlotly(
              subplot(rawplots[["p_bymz"]], ggplotly(ms1plot), nrows = 2,
                      shareX = TRUE)
          )
          output$rawMSMS_bygroup <- renderPlotly(
              subplot(rawplots[["p_bygroup"]], ggplotly(ms1plot), nrows = 2,
                      shareX = TRUE)
          )
          
          output$rawpks <- renderTable(rawplots[["pks"]])
          if (is(rawplots[["net"]], "visNetwork")) {
            output$rawss <- renderVisNetwork(rawplots[["net"]])
          }
        }, error = function(e){message("Raw MS2 plot failed")})

      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -1)

    observeEvent({
      input$selectident
      input$id
    },{
    tryCatch({
      ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
      sslist <- ms2@Ident[["MS2Features"]]
      if(is.data.frame(sslist) & !is.na(as.numeric(input$selectident))){
        if(as.numeric(input$selectident) <= nrow(sslist)){
          cur <- sslist$results[[as.numeric(input$selectident)]]
          if(is.data.frame(cur)){
            updatePickerInput(session, "selecthits", choices = cur$formula,
                              selected = cur$formula[[1]])
          }
        }
      }
    }, error = function(cond){})
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 50)

    #Identification plot (against MS2 DB hits)
    observeEvent({
      input$selectident
      input$selecthits
      input$id
    }, {
      if (input$selectident != "" & !is.na(input$selectident) &
          length(input$selecthits) != 0) {
        tryCatch({
            ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
            sslist <- ms2@Ident[["MS2Features"]]
            with_hits <- which(vapply(sslist$results,
                                      function(x){is.data.frame(x)},
                                      logical(1)))
            if(as.numeric(input$selectident) %in% with_hits){
              identplots <- plotMirror(struct$dataset,
                                       as.numeric(input$id), as.numeric(input$selectident),
                                       input$selecthits, mode = "hits")
              output$mirrorplot <- renderPlotly(identplots)
            }
        }, error = function(cond){warning("Identity plot failed")})
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -2)

    #Pairwise MS2 spectrum plot
    observeEvent({
      input$selectident2
      input$otherss
      input$id
    }, {
      if (input$selectident2 != "" & !is.na(input$selectident2) &
          length(input$otherss) != 0) {
        tryCatch({
          ms2 <- struct$dataset@data@MS2Exp[[as.numeric(input$id)]]
          sslist <- ms2@Ident[["MS2Features"]]
          identplots <- plotMirror(struct$dataset,
                                   as.numeric(input$id), as.numeric(input$selectident2),
                                   as.numeric(input$otherss), mode = "versus")
          output$versusplot <- renderPlotly(identplots)
        }, error = function(e){message("Mirrorplot failed")})
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -2)

    #Single MS2 spectrum plot
    observeEvent({
      input$selectss
      input$id
    }, {
      if (!is.na(input$selectss) & input$selectss != "") {
        identplots <- plotSS(struct$dataset, as.numeric(input$id),
                             as.numeric(input$selectss))
        output$ssplot <- renderPlotly(identplots)
      }
    }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = -2)
  })
}

plotMirror <- function(struct, id, ssnumber, patform, mode = "hits") {
    entry <- struct@data@MS2Exp[[id]]@Ident$MS2Features[ssnumber,]
    query <- entry$ssdata[[1]]
    maxint <- max(query$int)
    query$int <- query$int/max(query$int) * 100
    bestdf <- query[query$int > 10, ]
    bestdf$mz <- round(bestdf$mz, 4)
    molecmass <- entry$precmass
    baseline <- 1000
    subtitle <- ""
    title <- ""
    baseline <- baseline/maxint * 100

    f <- list(
        family = "Open Sans",
        size = 16,
        color = "black"
    )

    if(mode == "hits"){
        ref <- struct@data@MS2Exp[[id]]@Ident$MS2_correspondance[[ssnumber]]
        pattern <- struct@data@MS2Exp[[id]]@Ident$DatabaseSpectra[ref]
        pattern <- unname(pattern) #Avoids "name.subname" when unlisting
        pattern <- unlist(pattern, recursive = FALSE, use.names = TRUE)
    }

    mirrplot <- lapply(patform, function(x){
        pl <- ggplot()
        if(mode == "hits"){
            index <- entry$results[[1]]$id[entry$results[[1]]$formula == x]
            if(length(index) == 0){
                refspec <- NULL
            } else {
                refspec <- pattern[[x]][[index]]
                spec_energy <- names(pattern[[x]])[index]
                comp_name <- strsplit(patform, split = "#")[[which(patform == x)[1]]][3]
                a <- list(
                    text = paste(comp_name, spec_energy),
                    font = f,
                    xref = "paper",
                    yref = "paper",
                    yanchor = "bottom",
                    xanchor = "center",
                    align = "center",
                    x = 0.5,
                    y = 1,
                    showarrow = FALSE
                )
            }
            if(is.null(refspec)){return(ggplotly(pl))}
            refspec <- as.data.frame(t(refspec))
            pl <- pl + scale_x_continuous(limits = c(min(c(refspec$mz,
                                                    query$mz, molecmass)) - 20,
                                                    max(c(refspec$mz, query$mz,
                                                        molecmass)) + 20))
        }
        if(mode == "versus"){
            refmass<-struct@data@MS2Exp[[id]]@Ident$MS2Features$precmass[[x]]
            refspec <- struct@data@MS2Exp[[id]]@Ident$MS2Features$ssdata[[x]]
            if(is.null(refspec)){return(ggplotly(pl))}

            a <- list(
                text = paste("SS comparison between", ssnumber, "and", x),
                font = f,
                xref = "paper",
                yref = "paper",
                yanchor = "bottom",
                xanchor = "center",
                align = "center",
                x = 0.5,
                y = 1,
                showarrow = FALSE
            )

            refdf <- data.frame(mz = refmass)
            pl <- pl +
                geom_point(data = refdf, aes(x = .data$mz, y = 0), shape = 25,
                            size = 4, color = "red", fill = "red") +
                scale_x_continuous(limits = c(min(c(refspec$mz, query$mz,
                                                    molecmass, refmass)) - 20,
                                    max(c(refspec$mz, query$mz, molecmass,
                                            refmass)) + 20))
        }
            colnames(refspec) <- c("mz", "int")
            refspec$int <- refspec$int/max(refspec$int) * 100

            moldf <- data.frame(mz = molecmass)
            bldf <- data.frame(xmin = min(c(refspec$mz, query$mz,
                                        molecmass)) - 5,
                                xmax = max(c(refspec$mz,query$mz,
                                        molecmass) + 5),
                                y = baseline)

            pl <- pl + geom_segment(data = query, aes(x = .data$mz,
                        xend = .data$mz, y = 0, yend = .data$int),
                        color = "black") +
                    geom_segment(data = refspec, aes(x = .data$mz,
                        xend = .data$mz, y = 0, yend = -.data$int),
                        color = "red") +
                    geom_segment(data = bldf, aes(x = .data$xmin,
                        xend = .data$xmax, y = .data$y, yend = .data$y),
                        linetype = "dashed", color = "black", alpha = 0.3) +
                    geom_segment(data = bldf, aes(x = .data$xmin,
                        xend = .data$xmax, y = -.data$y, yend = -.data$y),
                        linetype = "dashed", color = "red", alpha = 0.3) +
                    geom_point(data = moldf, aes(x = .data$mz, y = 0),
                                shape = 17, size = 2) +
                    theme_minimal() + ylab("% Intensity") +
                    theme(plot.margin = unit(c(1, 0.7, 1, 0.8), "cm"),
                            text = element_text(size = 11,
                                                family = "Segoe UI Light"),
                            plot.title = element_text(hjust = 0.5)) +
                    geom_text(data = bestdf, aes(x = .data$mz,
                                                    y = .data$int + 5,
                                                    label = .data$mz),
                                family = "Segoe UI Light", check_overlap = TRUE)

            base_height <- ifelse(length(patform)<5, 850/length(patform), 200)
            ggplotly(pl, height = base_height * length(patform)) %>%
                layout(annotations = a)
        })

    return(subplot(mirrplot, nrows = length(mirrplot), shareX = TRUE,
                    which_layout = 1))
}

