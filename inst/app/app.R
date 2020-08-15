library(shiny)
library(shinyFiles)
library(shinydashboard)
library(shinyWidgets)
library(RHermes)
library(data.table, quietly = TRUE)
library(enviPat)
library(doSNOW)
library(foreach)
library(parallel)
library(tidyverse)
require(enviPat)
library(igraph)
require(mzR)
require(magrittr)
require(doParallel)
require(CHNOSZ)
require(ggplot2)
library(keras)
library(plotly)
library(visNetwork)
library(KEGGREST)
library(slickR)

header <- dashboardHeader(
  title = p("RHermes", style = "font-family: Open Sans, -apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Arial; font-size: 22px"),
  dropdownMenuOutput(outputId = "dropdown_m")
  )

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Intro", tabName = "intro", icon = icon("info")),
    menuItem("From mzML to PL", tabName = "mz2PL", icon = icon("bolt")),
    menuItem("SOI detection", tabName = "SOIdet", icon = icon("barcode")),
    menuItem("Inclusion list generation", tabName = "IL", icon = icon("list")),
    menuItem("MSMS data processing", tabName = "MSMS", icon = icon("atom")),
    menuItem("Identifications", tabName = "ident", icon = icon("fingerprint")),
    menuItem("Plots", tabName = "Plots", icon = icon("chart-bar"),
             menuSubItem("PL Exploration", tabName = "PLplot"),
             menuSubItem("SOI Exploration", tabName = "SOIplot"),
             menuSubItem("MS2 Exploration", tabName = "MS2plot")
    ),
    menuItem("Tables", tabName = "extra", icon = icon("table")),
    menuItem("Settings", tabName = "sett", icon = icon("file-alt")),
    menuItem("Loading/Saving", tabName = "IO", icon = icon("save")),
    div(actionButton("shutdown", label="Shutdown app", icon = icon("power-off"),
                 style ="background-color: #222d32; color: #b8c7ce; border-color: #b8c7ce; border-width:1px"),
        style = "margin-top:100%; width: 220px; display: inline-grid; position:fixed; bottom:0")
  )
)

body <- dashboardBody(
  tags$head(tags$link(rel="stylesheet", href="https://cdnjs.cloudflare.com/ajax/libs/animate.css/4.0.0/animate.min.css")), ##Adds animate.css animations
  tags$script(HTML("$('body').addClass('fixed');")),
  includeCSS("./www/customStyle.css"),
  tabItems(
    # First tab content
    tabItem(tabName = "intro",
            titlePanel(
            div(class="animate__animated animate__pulse", h1(tags$b("RHermes"), align = "center"), style = "text-align: center;")),
            mainPanel(width = 20,
              p("RHermes is a bioinformatics tool designed to metabolically characterize biological samples of any origin. It works by throwing a database of possible m/z associated to different formulas and adducts onto LC-MS spectra and determine which formulas are present and where they appear on the chromatogram. These results are then condensed into an inclusion list of ions to analyze with MS/MS fragmentation."),
              p("Combined with this fragmentation, RHermes enables LC-MS spectra \"deconvolution\" and characterization of virtually all ionizable species in a sample."),
              hr(),
              h2("General outline of the workflow:"),
              br(),
              div(
                slickR(c("./www/Step1.png", "./www/Step2.png", "./www/Step3.png"),
                       height = "500px") +
                settings(dots = TRUE, accessibility = TRUE, infinite = TRUE),
                style = "height: 600px"
              ),
              p("This GUI is split into the main steps of the RHermes workflow."),
              h2("Troubleshooting"),
              style = "padding : 20px 50px 50px 50px")
    ),

    # Second tab content
    tabItem(tabName = "mz2PL",
      PL_UI("PL_UI")
    ),
    tabItem(tabName = "SOIdet",
      SOI_UI("SOI_UI"),
    ),
    tabItem(tabName = "IL",
      IL_UI("IL_UI"),
    ),
    tabItem(tabName = "MSMS",
      MS2_UI("MS2_UI")
            ),
    tabItem(tabName = "Plots",
            verticalLayout(
              )),
    tabItem(tabName = "PLplot",  RHermes:::PLPlotUI("PLPlotUI")),
    tabItem(tabName = "SOIplot",  RHermes:::SOIPlotUI("SOIPlotUI")),
    tabItem(tabName = "MS2plot",  RHermes:::MS2PlotUI("MS2PlotUI")),
    tabItem(tabName = "ident", RHermes:::Ident_UI("Identifications")),
    tabItem(tabName = "extra", RHermes:::ExtraInfo_UI("ExtraInfo_UI")),
    tabItem(tabName = "sett", RHermes:::Settings_UI("Settings_UI")),
    tabItem(tabName = "IO", RHermes:::IO_UI("IO_UI"))
  )
)

ui <- dashboardPage(header, sidebar, body, title = "RHermes")

server <- function(input, output, session){

  struct <- reactiveValues(dataset=RHermesExp(), hasPL = FALSE, hasSOI = FALSE,
                           hasIL = FALSE, hasMS2 = FALSE)

  observe({
    struct$hasPL <- ifelse(length(struct$dataset@data@PL)!=0, TRUE, FALSE)
    struct$hasSOI <- ifelse(length(struct$dataset@data@SOI)!=0, TRUE, FALSE)
    struct$hasIL <- ifelse(length(struct$dataset@data@MS2Exp)!=0, TRUE, FALSE)
    if(struct$hasIL){
      struct$hasMS2 <- any(vapply(struct$dataset@data@MS2Exp, function(x){
        length(x@Ident)!=0
      }, logical(1)))    }
  })

  output$dropdown_m <- renderMenu({
    dropdownMenu(type = "tasks", badgeStatus = "primary",
                 taskItem(value = ifelse(struct$hasPL, 100, 0), color = "green",
                          "The object has PL/s"
                          ),
                 taskItem(value = ifelse(struct$hasSOI, 100, 0), color = "aqua",
                          "The object has SOI list/s"
                          ),
                 taskItem(value = ifelse(struct$hasIL, 100, 0), color = "yellow",
                          "The object has IL/s"
                          ),
                 taskItem(value = ifelse(struct$hasMS2, 100, 0), color = "red",
                          "The object has MS2 data"
                          )
  )})

  PLresults <- RHermes:::PLServer("PL_UI", struct = struct)
  observeEvent(PLresults$trigger, {
    struct$dataset <- PLresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  SOIresults <- RHermes:::SOIServer("SOI_UI", struct = struct)
  observeEvent(SOIresults$trigger, {
    struct$dataset <- SOIresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  ILresults <- RHermes:::ILServer("IL_UI", struct = struct)
  observeEvent(ILresults$trigger, {
    struct$dataset <- ILresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  MS2results <- RHermes:::MS2Server("MS2_UI", struct)
  observeEvent(MS2results$trigger, {
    struct$dataset <- MS2results$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  IOresults <- RHermes:::IOServer("IO_UI", struct = struct)
  observeEvent(IOresults$trigger, {
    struct$dataset <- IOresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  RHermes:::PLPlotServer("PLPlotUI", struct = struct)
  RHermes:::SOIPlotServer("SOIPlotUI", struct = struct)
  RHermes:::MS2PlotServer("MS2PlotUI", struct = struct)
  RHermes:::IdentServer("Identifications", struct = struct)
  RHermes:::ExtraInfoServer("ExtraInfo_UI", struct = struct)

  setResults <- RHermes:::SettingsServer("Settings_UI", struct = struct)
  observeEvent(setResults$trigger, {
    struct$dataset <- setResults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$shutdown, {shiny::stopApp()})
}

# setmessage <- function(mobj, m){
#   mobj <- rbind(mobj, data.frame(from = "The Developer", message = m, date = as.character(date()),
#                                   stringsAsFactors = FALSE))
#   return(mobj)
# }
# popmessage <- function(mobj){return(mobj[-1,])}


shinyApp(ui, server)
