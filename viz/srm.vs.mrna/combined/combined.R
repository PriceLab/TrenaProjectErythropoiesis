library(shiny)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
load("../docker/srm.rna.averaged.clean.RData")
max.time.points <- 13
goi <- rownames(mtx.rna)
#------------------------------------------------------------------------------------------------------------------------
srm.rna.tab <- function()
{
   sidebarLayout(
      sidebarPanel(
         radioButtons("transformChoice", "Data Transform",
                      c("None", "Normalized", "Arcsinh")),
         selectInput("geneSelector", "Single TF: rna + srm", goi, selected=goi[1],  multiple=FALSE),
         width=2
         ),
      mainPanel(
         d3Output("srm.rna.d3", height="800px"),
         width=10
         )
      )

} # srm.rna.tab
#------------------------------------------------------------------------------------------------------------------------
srm.coexpression.tab <- function()
{
   sidebarLayout(
      sidebarPanel(width = 3,
                   sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 30)
                   ),
      mainPanel(plotOutput("distPlot"))
      ) # sidebarLayout


} # srm.coexpression.tab
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(

   titlePanel("SRM and RNA-seq in erythropoiesis"),

   tabsetPanel(
       tabPanel("SRM & rna compared", srm.rna.tab()),
       tabPanel("SRM co-expression", srm.coexpression.tab())
       ) # tabsetPanel

   ) # fluidPage

#------------------------------------------------------------------------------------------------------------------------
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

  output$srm.rna.d3 <- renderD3({
     transform <- input$transformChoice
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
app <- shinyApp(ui = ui, server = server)

