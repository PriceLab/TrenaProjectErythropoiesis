library(shiny)
library(r2d3)
library(rsconnect)
addResourcePath("www", "www");
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
load("srm.rna.averaged.clean.RData")
max.time.points <- 13
goi <- rownames(mtx.rna)
# these three proteins have spotty srm data, though good rna-seq.
# marjorie and jeff ask that they be eliminated in the protein/rna list
bad.proteins <- c("ETO2", "MLL1", "SPT16")
goi.for.comparison <- goi # setdiff(goi, bad.proteins)
# apply(mtx.srm[bad.proteins,], 1, function(row) length(which(is.na(row))))
#       ETO2  MLL1 SPT16
#          7     5     7
#------------------------------------------------------------------------------------------------------------------------
srm.rna.tab <- function()
{
   sidebarLayout(
      sidebarPanel(
         selectInput("geneSelector", "Plot Protein and mRNA",
                     goi.for.comparison,
                     selected=goi.for.comparison[1],  multiple=FALSE),
         radioButtons("srm.rna.lineTypeSelector", "Smoothing", c("No", "Yes")),
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
         radioButtons("srm.transformChoice", "Data Transform", c("None", "Normalized")), # , "Arcsinh")),
         radioButtons("srm.lineTypeSelector", "Smoothing", c("No", "Yes")),
         selectInput("srmSelector", "Plot Protein", goi, selected=NULL,  multiple=FALSE),
         sliderInput("correlationThresholdSlider", label = "Pearson", min = 0, max = 1, value = 0.9, step = 0.01),
         radioButtons("correlationDirectionChooser", "Find Correlations", c("None", "+", "-")),
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

   tags$head(
        tags$style(".tab-pane {margin-top: 20px;}"),
        tags$link(rel = "stylesheet", type = "text/css", href = "www/app.css")
        ),

   titlePanel("Transcription Factor Protein and RNA Expression Profiles During Erythropoiesis"),

   tabsetPanel(
       tabPanel("Protein & RNA compared", srm.rna.tab()),
       tabPanel("Temporal Protein Abundances", srm.coexpression.tab())
       ) # tabsetPanel

   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session) {

   reactiveState <- reactiveValues(selectedTF=NULL, correlatedTFs=list())

   observeEvent(input$srmSelector, ignoreInit=FALSE, {
      plotCorrelatedProteins(input, output)
      })

   observeEvent(input$geneSelector, ignoreInit=TRUE, {
      tf <- input$geneSelector
      })

   observeEvent(input$srm.transformChoice, ignoreInit=TRUE, {
     plotCorrelatedProteins(input, output)
     })

   observeEvent(input$correlationThresholdSlider, ignoreInit=TRUE, {
      plotCorrelatedProteins(input, output)
      })

   observeEvent(input$correlationDirectionChooser, ignoreInit=TRUE, {
      plotCorrelatedProteins(input, output)
      })

  output$srm.rna.d3 <- renderD3({
     lineSmoothing <- input$srm.rna.lineTypeSelector
     r2d3.command <- "plotBoth"
     xValues <- as.numeric(sub("d_", "", colnames(mtx.srm)))
     xMax <- max(xValues)
     yMax <- 1.0
     tf <- input$geneSelector[1]
     printf("plotting srm+rna for %s (%s)", tf, lineSmoothing)
     timepoints <- as.numeric(sub("d_", "", colnames(mtx.rna)))
     rna.values <- as.numeric(mtx.rna[tf,])
     srm.values <- as.numeric(mtx.srm[tf,])

     vectors <- transformData.rna.srm(rna.values, srm.values, transformName="None")
     rna.values <- vectors[["rna"]]
     srm.values <- vectors[["srm"]]

     xMin <- min(timepoints)
     xMax <- max(timepoints)
     yMin <- 0
     yMax <- max(c(rna.values, srm.values), na.rm=TRUE)
     y2Max <- max(rna.values, na.rm=TRUE)

     rna.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=rna.values[i])))
     srm.xy <- lapply(seq_len(length(timepoints)), function(i) return(list(x=timepoints[i], y=srm.values[i])))

     data <- list(rna=rna.xy, srm=srm.xy, xMax=xMax, yMax=yMax, y2Max=y2Max, cmd=r2d3.command, smoothing=lineSmoothing)
     r2d3(data, script = "linePlot.js")
     })

   observeEvent(input$currentlySelectedVector, ignoreInit=FALSE, {
     newValue <- input$currentlySelectedVector
     if(nchar(newValue) == 0)
        newValue <- "   "
     output$currentVectorDisplay <- renderText({newValue})
     })


} # server
#------------------------------------------------------------------------------------------------------------------------
transformData.rna.srm <- function(rna, srm, transformName)
{
   # printf("--- transform by %s", transformName)

   if(transformName == "None"){
      rna.out <- rna;
      srm.out <- srm;
      }

   if(transformName == "Normalized"){
      rna.out <- rna/max(rna, na.rm=TRUE)
      srm.out <- srm/max(srm, na.rm=TRUE);
      }

   if(transformName == "Arcsinh"){
      rna.out <- asinh(rna)
      srm.out <- asinh(srm)
      }

   return(list(rna=rna.out, srm=srm.out))

} # transformData.rna.srm
#------------------------------------------------------------------------------------------------------------------------
findCorrelated <- function(targetTF, threshold, direction)
{
   if(direction == "None")
      return(targetTF)

   suppressWarnings(
      correlations <- apply(mtx.srm, 1,
                           function(row) cor(mtx.srm[targetTF,], row,  use="complete.obs"))
      )

   if(direction == "-")
      result <- names(which(correlations <= (-1 * threshold)))
   else  # must be "+"
      result <- names(which(correlations >= threshold))

   return(unique(c(targetTF, result)))

} # findCorrelated
#------------------------------------------------------------------------------------------------------------------------
plotCorrelatedProteins <- function(input, output)
{
   tf <- input$srmSelector
   correlationThreshold <- input$correlationThresholdSlider;
   correlationDirection <- isolate(input$correlationDirectionChooser)
   tfs.all <- findCorrelated(tf, correlationThreshold, correlationDirection)
   transform <- input$srm.transformChoice

   output$srm.d3 <- renderD3({
     plotTFs(tfs.all, input, output, transform)
     })

} # plotCorrelatedProteins
#------------------------------------------------------------------------------------------------------------------------
plotTFs <- function(tfs, input, output, transform)
{

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

   lineSmoothing <- input$srm.lineTypeSelector
   printf("plotTFs (%s, %s): %s", transform, lineSmoothing, paste(tfs, collapse=", "))

   data <- list(vectors=vectorsWithTimes, xMax=xMax, yMax=yMax, cmd="plot", smoothing=lineSmoothing)
   # if("ETO2" %in% tfs) browser()
   r2d3(data, script = "multiPlot.js")

} # plotTFs
#------------------------------------------------------------------------------------------------------------------------
maxOfVectors <- function(vectorList)
{
   max <- 0
   for(vector in vectorList){
      vector.max <- max(vector, na.rm=TRUE)
      #if(is.na(vector.max)) browser()
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
      srm.out <- lapply(srm, function(vec) vec/max(vec, na.rm=TRUE))
      }

   if(transformName == "Arcsinh"){
      srm.out <- lapply(srm, asinh)
      }

   return(srm.out)

} # transformData.srm
#------------------------------------------------------------------------------------------------------------------------
deploy <- function(){
    deployApp(account="pshannon", appName="tf-srm-rna",
              appFiles=c("app.R",
                         "linePlot.css",
                         "linePlot.js",
                         "multiPlot.css",
                         "multiPlot.js",
                         "srm.rna.averaged.clean.RData",
                         "www/app.css")

              )
} # deploy
#------------------------------------------------------------------------------------------------------------------------
#runApp(shinyApp(ui = ui, server = server), port=8888)
shinyApp(ui = ui, server = server)


