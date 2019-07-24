library(shiny)
library(devtools)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
load("~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/srm.rna.averaged.clean.RData")

#mtx.srm <- get(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/SRMforPublication20190614.RData"))
#mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
max.time.points <- 13
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
   tags$head(tags$style("#d3{height:90vh !important;}")),
   titlePanel("mRNA vs protein counts"),
   sidebarLayout(
      sidebarPanel(
         selectInput("geneSelector", "", head(rownames(mtx.rna)), selected=rownames(mtx.rna)[1],  multiple=FALSE),
         actionButton("forwardTimeStepButton", "+", style="margin-bottom: 20px; margin-left: 100px; font-size:200%"),
         actionButton("backwardTimeStepButton", "-", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),
         verbatimTextOutput("timeStepDisplay"),
         #sliderInput("dayNumberSlider", label = "Day:", min = 0, max = 12, value = 0, step = 1, round=TRUE),
         width=2
         ),
      mainPanel(
         d3Output("d3", height="800px"),
         width=10
         )
      )
   ) # fluidPage

server <- function(input, output) {

  reactiveState <- reactiveValues(timeStep=1, genes=list())

  observeEvent(input$geneSelector, ignoreInit=TRUE, {
     tf <- input$geneSelector
     printf("new tf: %s", tf)
     })

  output$timeStepDisplay <- renderText({
      reactiveState$timeStep
      })

  observeEvent(input$forwardTimeStepButton, ignoreInit=FALSE, {
     currentValue <- reactiveState$timeStep
     if(currentValue == max.time.points) return()
     reactiveState$timeStep <- reactiveState$timeStep + 1
     })

  observeEvent(input$backwardTimeStepButton, ignoreInit=FALSE, {
     currentValue <- reactiveState$timeStep
     if(currentValue == 1) return()
     reactiveState$timeStep <- reactiveState$timeStep - 1
     })

  output$d3 <- renderD3({
     currentDay <- reactiveState$timeStep
     if(currentDay <= 0) return();
     if(currentDay > max.time.points) return();
     xValues <- as.numeric(sub("d_", "", colnames(mtx.srm)))
     xMax <- max(xValues)
     yMax <- 1.0
     tf <- input$geneSelector[1]
     rna <- as.numeric(mtx.rna[tf, ]) # c(1,3,5,7,9,11,13,15,17,19,21,23,25)])
     #rna <- rna/max(rna)
     srm <- as.numeric(mtx.srm[tf,])
     #srm <- srm/max(srm)
     timePoints <- as.numeric(sub("d_", "", colnames(mtx.rna)))
     rna.values <- as.numeric(mtx.rna[tf,])
     printf("---- rna.values")
     srm.values <- as.numeric(mtx.srm[tf,])
     xMin <- min(timePoints)
     xMax <- max(timePoints)
     yMin <- 0
     yMax <- max(c(rna.values, srm.values))

     rna.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=rna.values[i])))
     srm.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=srm.values[i])))

     data <- list(rna=rna.xy, srm=srm.xy, xMax=xMax, yMax=yMax)
     browser()
     r2d3(data, script = "linePlot.js")
    })
}

app <- shinyApp(ui = ui, server = server)
