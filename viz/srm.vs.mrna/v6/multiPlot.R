library(shiny)
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
         #selectInput("geneSelector", "Single TF: rna + srm", goi, selected=goi[1],  multiple=FALSE),
         selectInput("srmSelector", "", goi, selected=NULL,  multiple=TRUE, size=20, selectize=FALSE),
         sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
         actionButton("findPositiveCorrelationsButton", "Find Correlated +", style="margin-bottom: 20px; margin-left: 2px; font-size:100%"),
         actionButton("findNegativeCorrelationsButton", "Find Correlated -", style="margin-bottom: 20px; margin-left: 2px; font-size:100%"),
         #actionButton("backwardTimeStepButton", "-", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),
         #verbatimTextOutput("timeStepDisplay"),
         #actionButton("clearPlotButton", "Clear", style="margin-bottom: 20px; margin-left: 10px; font-size:200%"),
         #sliderInput("dayNumberSlider", label = "Day:", min = 0, max = 12, value = 0, step = 1, round=TRUE),
         br(),
         verbatimTextOutput(outputId="currentVectorDisplay"),
         width=2
         ),
      mainPanel(
         d3Output("d3", height="800px"),
         width=10
         )
      )
   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(session, input, output) {

   reactiveState <- reactiveValues(timeStep=1, genes=list())

   observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
     newValue <- input$currentlySelectedVector
     # printf("newValue: %s", newValue)
     if(nchar(newValue) == 0)
        newValue <- "   "
     output$currentVectorDisplay <- renderText({newValue})
     #output$currentVectorDisplay <- renderText({newValue})
     })

  observeEvent(input$transformChoice, ignoreInit=TRUE, {
     tfs <- input$srmSelector
     currentDay <- reactiveState$timeStep
     transform <- input$transformChoice
     output$d3 <- renderD3({
       plotTFs(tfs, input, output, transform)
       })
     })

  observeEvent(input$srmSelector, ignoreInit=TRUE, {
     tfs <- input$srmSelector
     currentDay <- reactiveState$timeStep
     transform <- input$transformChoice
     output$d3 <- renderD3({
       plotTFs(tfs, input, output, transform)
       })
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

   observeEvent(input$findPositiveCorrelationsButton, ignoreInit=TRUE,{
     tfs <- isolate(input$srmSelector)
     threshold <- isolate(input$correlationThresholdSlider)
     tfs.correlated <- findCorrelated(tfs[1], threshold)
     updateSelectInput(session, "srmSelector", selected=tfs.correlated)
     })

   observeEvent(input$findNegativeCorrelationsButton, ignoreInit=TRUE,{
     tfs <- isolate(input$srmSelector)
     threshold <- isolate(input$correlationThresholdSlider)
     tfs.correlated <- findCorrelated(tfs[1], threshold, negative=TRUE)
     updateSelectInput(session, "srmSelector", selected=tfs.correlated)
     })

  #  output$d3 <- renderD3({
  #     transform <- isolate(input$transformChoice)
  #     r2d3.command <- "plotBoth"
  #     currentDay <- reactiveState$timeStep
  #     if(currentDay <= 0) return();
  #     if(currentDay > max.time.points) return();
  #     xValues <- as.numeric(sub("d_", "", colnames(mtx.srm)))
  #     xMax <- max(xValues)
  #     yMax <- 1.0
  #     proceed <- input$plotTFsButton
  #     tfs <- isolate(input$srmSelector)
  #     timePoints <- as.numeric(sub("d_", "", colnames(mtx.srm)))
  #     srm.vectors <- lapply(tfs, function(tf) as.numeric(mtx.srm[tf,]))
  #     names(srm.vectors) <- tfs
  #     # vectors <- transformData(rna.values, srm.values, transform)
  #     # rna.values <- vectors[["rna"]]
  #     # srm.values <- vectors[["srm"]]
  #
  #     xMin <- min(timePoints)
  #     xMax <- max(timePoints)
  #     yMin <- 0
  #     yMax <- maxOfVectors(srm.vectors)
  #
  #     vectorsWithTimes <- vector(mode="list", length(tfs))
  #     names(vectorsWithTimes) <- tfs
  #
  #     for(tf in tfs){
  #        srm.vector <- srm.vectors[[tf]]
  #        vectorsWithTimes[[tf]] <- lapply(seq_len(length(timePoints)), function(i) return(list(x=timePoints[i], y=srm.vector[i])))
  #        }
  #
  #     data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd=r2d3.command)
  #     r2d3(data, script = "multiPlot.js")
  #    })

} # server
#------------------------------------------------------------------------------------------------------------------------
findCorrelated <- function(targetTF, threshold, negative=FALSE)
{
   correlations <- apply(mtx.srm, 1, function(row) cor(mtx.srm[targetTF,], row))
   # browser()

   if(negative)
      result <- names(which(correlations <= (-1 * threshold)))
   else
      result <- names(which(correlations >= threshold))

   return(unique(c(targetTF, result)))

} # findCorrelated
#------------------------------------------------------------------------------------------------------------------------
plotTFs <- function(tfs, input, output, transform)
{
   timePoints <- as.numeric(sub("d_", "", colnames(mtx.srm)))
   srm.vectors <- lapply(tfs, function(tf) as.numeric(mtx.srm[tf,]))
   names(srm.vectors) <- tfs

   srm.vectors <- transformData(srm.vectors, transform)

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

   #browser()

   data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd="plot")

   #printf("calling r2d3")
   #print(data)
   r2d3(data, script = "multiPlot.js")

} # plotTFs
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
transformData <- function(srm, transformName)
{
   # printf("--- transform by %s", transformName)

   if(transformName == "None"){
      srm.out <- srm;
      }

   if(transformName == "Normalized"){
      srm.out <- lapply(srm, function(vec) vec/max(vec))
      }

   if(transformName == "Arcsinh"){
      srm.out <- lapply(srm, asinh)
      }

   return(srm.out)

} # transformData
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)
