library(shiny)
library(devtools)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
load("~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/srm.rna.averaged.clean.RData")
# mtx.rna <- asinh(mtx.rna)
# mtx.srm <- asinh(mtx.srm)

#mtx.srm <- get(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/SRMforPublication20190614.RData"))
#mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
max.time.points <- 13
#goi <- c("GATA1", "SPI1", head(rownames(mtx.rna), n=20))
goi <- rownames(mtx.rna)
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
   tags$head(
      #tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/d3-tip/0.9.1/d3-tip.min.js"),
      tags$style("#d3{height:90vh !important;}")),
   titlePanel("srm timecourses"),
   sidebarLayout(
      sidebarPanel(
         radioButtons("transformChoice", "Data Transform",
                      c("None", "Normalized", "Arcsinh")),
         verbatimTextOutput(outputId="currentVectorDisplay"),
         #selectInput("geneSelector", "Single TF: rna + srm", goi, selected=goi[1],  multiple=FALSE),
         selectInput("srmSelector", "", goi, selected=goi[1],  multiple=TRUE, size=20, selectize=FALSE),
         actionButton("plotTFsButton", "Plot", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),

         sliderInput("correlationThresholdSlider", label = "abs(pearson)", min = 0, max = 1, value = 0.9, step = 0.01),
         #actionButton("forwardTimeStepButton", "+", style="margin-bottom: 20px; margin-left: 20px; font-size:200%"),
         #actionButton("backwardTimeStepButton", "-", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),
         #verbatimTextOutput("timeStepDisplay"),
         #actionButton("clearPlotButton", "Clear", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),
         #sliderInput("dayNumberSlider", label = "Day:", min = 0, max = 12, value = 0, step = 1, round=TRUE),
         width=2
         ),
      mainPanel(
         d3Output("d3", height="800px"),
         width=10
         )
      )
   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output) {

  reactiveState <- reactiveValues(timeStep=1, genes=list())

   observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
     newValue <- input$currentlySelectedVector
     printf("newValue: %s", newValue)
     if(nchar(newValue) == 0)
        newValue <- "   "
     output$currentVectorDisplay <- renderText({newValue})
     #output$currentVectorDisplay <- renderText({newValue})
     })

  observeEvent(input$srmSelector, ignoreInit=TRUE, {
     tf <- input$srmSelector
     printf("new tf: %s", paste(tf, collapse=","))
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
     transform <- isolate(input$transformChoice)
     r2d3.command <- "plotBoth"
     currentDay <- reactiveState$timeStep
     if(currentDay <= 0) return();
     if(currentDay > max.time.points) return();
     xValues <- as.numeric(sub("d_", "", colnames(mtx.srm)))
     xMax <- max(xValues)
     yMax <- 1.0
     proceed <- input$plotTFsButton
     tfs <- isolate(input$srmSelector)
     timePoints <- as.numeric(sub("d_", "", colnames(mtx.srm)))
     srm.vectors <- lapply(tfs, function(tf) as.numeric(mtx.srm[tf,]))
     names(srm.vectors) <- tfs
     # vectors <- transformData(rna.values, srm.values, transform)
     # rna.values <- vectors[["rna"]]
     # srm.values <- vectors[["srm"]]

     xMin <- min(timePoints)
     xMax <- max(timePoints)
     yMin <- 0
     yMax <- maxOfVectors(srm.vectors)

     vectorsWithTimes <- vector(mode="list", length(tfs))
     names(vectorsWithTimes) <- tfs

     for(tf in tfs){
        srm.vector <- srm.vectors[[tf]]
        vectorsWithTimes[[tf]] <- lapply(seq_len(length(timePoints)), function(i) return(list(x=timePoints[i], y=srm.vector[i])))
        }

     data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd=r2d3.command)
     r2d3(data, script = "multiPlot.js")
    })

} # server
#------------------------------------------------------------------------------------------------------------------------
maxOfVectors <- function(vectorList)
{
   max <- 0
   for(vector in vectorList){
      vector.max <- max(vector, na.rm=TRUE)
      if(is.na(vector.max)) browser()
      if(vector.max > max)
         max <- vector.max
      } # for vector

   return(max)

} # maxOfVectors
#------------------------------------------------------------------------------------------------------------------------
transformData <- function(rna, srm, transformName)
{
   printf("--- transform by %s", transformName)

   if(transformName == "None"){
      rna.out <- rna;
      srm.out <- srm;
      }

   if(transformName == "Normalized"){
      rna.out <- rna/max(rna)
      srm.out <- srm/max(srm);
      }

   if(transformName == "Arcsinh"){
      rna.out <- asinh(rna)
      srm.out <- asinh(srm)
      }

   return(list(rna=rna.out, srm=srm.out))

} # transformData
#------------------------------------------------------------------------------------------------------------------------

app <- shinyApp(ui = ui, server = server)
