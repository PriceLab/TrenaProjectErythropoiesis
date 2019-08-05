library(shiny)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
load("srm.rna.averaged.clean.RData")
max.time.points <- 13
goi <- rownames(mtx.rna)
#------------------------------------------------------------------------------------------------------------------------
srm.rna.tab <- function()
{
   sidebarLayout(
      sidebarPanel(
         radioButtons("srm.rna.transformChoice", "Data Transform",
                      c("None", "Normalized", "Arcsinh")),
         selectInput("geneSelector", "Single TF", goi, selected=goi[1],  multiple=FALSE),
         span(style="color:red", "RNA"),
         span(" + "),
         span(style="color:blue", "SRM"),

         width=2
         ),
      mainPanel(
         d3Output("srm.rna.d3", height="80vh"),
         width=10
         )
      )

} # srm.rna.tab
#------------------------------------------------------------------------------------------------------------------------
srm.coexpression.tab <- function()
{
   sidebarLayout(
      sidebarPanel(
         radioButtons("srm.transformChoice", "Data Transform",
                      c("None", "Normalized", "Arcsinh")),
         selectInput("srmSelector", "", goi, selected=NULL,  multiple=TRUE, size=20, selectize=FALSE),
         sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
         actionButton("findPositiveCorrelationsButton", "Find Correlated +", style="margin-bottom: 20px; margin-left: 2px; font-size:100%"),
         actionButton("findNegativeCorrelationsButton", "Find Correlated -", style="margin-bottom: 20px; margin-left: 2px; font-size:100%"),
         br(),
         verbatimTextOutput(outputId="currentVectorDisplay"),
         width=2
         ),
      mainPanel(
         d3Output("srm.d3", height="80vh"),
         width=10
         )
      )

} # srm.coexpression.tab
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(

   tags$head(tags$style(".tab-pane {margin-top: 20px;}")),

   titlePanel("SRM and RNA-seq in erythropoiesis"),

   tabsetPanel(
       tabPanel("SRM & rna compared", srm.rna.tab()),
       tabPanel("SRM co-expression", srm.coexpression.tab())
       ) # tabsetPanel

   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session) {

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

  output$srm.rna.d3 <- renderD3({
     transform <- input$srm.rna.transformChoice
     r2d3.command <- "plotBoth"
     currentDay <- reactiveState$timeStep
     if(currentDay <= 0) return();
     if(currentDay > max.time.points) return();
     xValues <- as.numeric(sub("d_", "", colnames(mtx.srm)))
     xMax <- max(xValues)
     yMax <- 1.0
     tf <- input$geneSelector[1]
     timepoints <- as.numeric(sub("d_", "", colnames(mtx.rna)))
     rna.values <- as.numeric(mtx.rna[tf,])
     srm.values <- as.numeric(mtx.srm[tf,])

     vectors <- transformData.rna.srm(rna.values, srm.values, transform)
     rna.values <- vectors[["rna"]]
     srm.values <- vectors[["srm"]]

     xMin <- min(timepoints)
     xMax <- max(timepoints)
     yMin <- 0
     yMax <- max(c(rna.values, srm.values))

     rna.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=rna.values[i])))
     srm.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=srm.values[i])))

     data <- list(rna=rna.xy, srm=srm.xy, xMax=xMax, yMax=yMax, cmd=r2d3.command)
     # browser()
     r2d3(data, script = "linePlot.js")
     })


   reactiveState <- reactiveValues(timeStep=1, genes=list())

   observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
     newValue <- input$currentlySelectedVector
     # printf("newValue: %s", newValue)
     if(nchar(newValue) == 0)
        newValue <- "   "
     output$currentVectorDisplay <- renderText({newValue})
     #output$currentVectorDisplay <- renderText({newValue})
     })

  observeEvent(input$srm.transformChoice, ignoreInit=TRUE, {
     tfs <- input$srmSelector
     currentDay <- reactiveState$timeStep
     transform <- input$srm.transformChoice
     output$srm.d3 <- renderD3({
       plotTFs(tfs, input, output, transform)
       })
     })

  observeEvent(input$srmSelector, ignoreInit=TRUE, {
     tfs <- input$srmSelector
     currentDay <- reactiveState$timeStep
     transform <- input$srm.transformChoice
     output$srm.d3 <- renderD3({
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




} # server
#------------------------------------------------------------------------------------------------------------------------
transformData.rna.srm <- function(rna, srm, transformName)
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

} # transformData.rna.srm
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
   printf("plotTFs (%s): %s", transform, paste(tfs, collapse=", "))

   timePoints <- as.numeric(sub("d_", "", colnames(mtx.srm)))
   srm.vectors <- lapply(tfs, function(tf) as.numeric(mtx.srm[tf,]))
   names(srm.vectors) <- tfs

   srm.vectors <- transformData.srm(srm.vectors, transform)

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

   data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd="plot")

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
transformData.srm <- function(srm, transformName)
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

} # transformData.srm
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)

