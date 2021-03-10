for(i in list.files(system.file("app", package = "RHermes"), pattern = "UI",
                    full.names = TRUE)){
  source(i, local = TRUE)
}

library(shiny)
library(shinyFiles)
library(shinydashboard)
library(shinyWidgets)
library(RHermes)
library(data.table, quietly = TRUE)
library(enviPat)
library(igraph)
require(mzR)
require(magrittr)
require(CHNOSZ)
require(ggplot2)
library(keras)
library(plotly)
library(visNetwork)
library(KEGGREST)
library(slickR)
library(BiocParallel)
library(DT)

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
              HTML("<p style = font-size:1.1em>RHermes is a broad-scoped targeted metabolomics software designed to identify compounds in biological and environmental samples.
                We invert the classical <i>peak-first</i> metabolite detection workflows by first annotating the LC-MS1 raw data points themselves
                using a database of plasusible formulas and adducts. This ultimately results in a very sample-specific inclusion list of ions to monitor in a MS2 experiment</p>
                <p style = font-size:1.1em> RHermes has shown that <b>up to 90%</b> of a typical data-dependent acquisition (DDA) MS2 experiment is wasted in redundant and non-biological ions, leading to fake positive
                identifications while missing the true compounds in the sample.</p>"),
              HTML("<blockquote style = text-align:left><p> We believe there's a better way to characterize the metabolites in our samples</p> <p>- RHermes</p> </blockquote>"),
              HTML("<p style = font-size:1.1em>You can find the molecular formula databases and sample mzML files in the paper in our <a href = 'https://zenodo.org/record/4581662'>Zenodo dataset</a></p>"),
              hr(),
              HTML("<h2><b>The Workflow</b></h2>"),
              br(),
              tags$script(src = 'scroll.js'),
              HTML("<p style = font-size:1.1em> </p>"),
              img(src = "step1.svg", width = "100%", class = "steps"),
              HTML("<p style = font-size:1.1em> </p>"),
              img(src = "step2.svg", width = "100%", class = "steps"),
              HTML("<p style = font-size:1.1em> </p>"),
              img(src = "step3.svg", width = "100%", class = "steps"),
              hr(),
              h2("Bug reports"),
              HTML("<p style = font-size:1.1em>In case you find any abnormal behaviour (eg: app unexpectedly crashes, etc.) feel free to email
                   the mantainer at <a href = 'mailto:roger.gine@estudiants.urv'>roger.gine@estudiants.urv</a> or open a <a href = 'https://github.com/RogerGinBer/RHermes/issues'>Github issue</a>.
                   You should provide a detailed description about what happened, screenshots of the error and, <i>desirably</i>, a <a href = 'https://stackoverflow.com/help/minimal-reproducible-example'>minimally reproducible example</a></p>"),
              h2("Citation"),
              HTML("<p style = font-size:1.1em>Please cite this software as:</p>"),
              HTML("<blockquote style = text-align:left><p> HERMES: a molecular formula-oriented method to target the metabolome </p><p>
                    Roger Giné, Jordi Capellades, Josep M. Badia, Dennis Vughs, Michaela Schwaiger-Haber, Maria Vinaixa, Andrea M. Brunner, Gary J. Patti, Oscar Yanes</p><p>
                    bioRxiv 2021.03.08.434466; doi: https://doi.org/10.1101/2021.03.08.434466</p> </blockquote>"),
              style = "padding : 20px 50px 50px 50px")
    ),

    tabItem(tabName = "mz2PL", PL_UI("PL_UI")),
    tabItem(tabName = "SOIdet", SOI_UI("SOI_UI")),
    tabItem(tabName = "IL", IL_UI("IL_UI")),
    tabItem(tabName = "MSMS", MS2_UI("MS2_UI")),
    tabItem(tabName = "Plots", verticalLayout()),
    tabItem(tabName = "PLplot",  PLPlotUI("PLPlotUI")),
    tabItem(tabName = "SOIplot",  SOIPlotUI("SOIPlotUI")),
    tabItem(tabName = "MS2plot", MS2PlotUI("MS2PlotUI")),
    tabItem(tabName = "ident", Ident_UI("Identifications")),
    tabItem(tabName = "extra", ExtraInfo_UI("ExtraInfo_UI")),
    tabItem(tabName = "sett", Settings_UI("Settings_UI")),
    tabItem(tabName = "IO", IO_UI("IO_UI"))
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
                 taskItem(value = ifelse(struct$hasIL, 100, 0),color = "yellow",
                          "The object has IL/s"
                          ),
                 taskItem(value = ifelse(struct$hasMS2, 100, 0), color = "red",
                          "The object has MS2 data"
                          )
  )})

  PLresults <- PLServer("PL_UI", struct = struct)
  observeEvent(PLresults$trigger, {
      struct$dataset <- PLresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  SOIresults <- SOIServer("SOI_UI", struct = struct)
  observeEvent(SOIresults$trigger, {
    struct$dataset <- SOIresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  ILresults <- ILServer("IL_UI", struct = struct)
  observeEvent(ILresults$trigger, {
    struct$dataset <- ILresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  MS2results <- MS2Server("MS2_UI", struct)
  observeEvent(MS2results$trigger, {
    struct$dataset <- MS2results$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  IOresults <- IOServer("IO_UI", struct = struct)
  observeEvent(IOresults$trigger, {
    struct$dataset <- IOresults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  PLPlotServer("PLPlotUI", struct = struct)
  SOIPlotServer("SOIPlotUI", struct = struct)
  MS2PlotServer("MS2PlotUI", struct = struct)
  IdentServer("Identifications", struct = struct)
  ExtraInfoServer("ExtraInfo_UI", struct = struct)

  setResults <- SettingsServer("Settings_UI", struct = struct)
  observeEvent(setResults$trigger, {
    struct$dataset <- setResults$dataset
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$shutdown, {shiny::stopApp()})
}

shinyApp(ui, server)
