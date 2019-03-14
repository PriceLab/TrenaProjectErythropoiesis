library(shiny)
library(DT)
library(ggplot2)
library(grid)
#----------------------------------------------------------------------------------------------------
load("KLF1.FLI1.atac.srm.RData") # mtx
rownames(mtx)[1] <- sprintf("%s srm count", rownames(mtx)[1])
rownames(mtx)[2] <- sprintf("%s srm count", rownames(mtx)[2])
mtx[1, ] <- round(mtx[1,])
mtx[2, ] <- round(mtx[2,])

tbl <- as.data.frame(t(mtx))
tbl$time <- seq_len(nrow(tbl))
tbl$tp <- as.factor(rownames(tbl))
colnames(tbl) <- c("KLF1.srm", "FLI1.srm", "KLF1.total.hits", "KLF1.top.hits", "KLF1.total.density",
                   "KLF1.top.density", "FLI1.total.hits", "FLI1.top.hits", "FLI1.total.density",
                   "FLI1.top.density", "time", "tp")

x.axis.labels <- c("4r1", "4r2", "8", "10r1", "10r2", "11r1", "11r2", "12r1", "12r2")

#----------------------------------------------------------------------------------------------------
createPlotCountsTab <- function()
{
   fluidRow(selectInput(inputId = "countsChooser",
                        label = NULL, #"Choose Row of interest:",
                        choices = colnames(tbl)),
            plotOutput("countsPlotter"),
            column=12)

} # createPlotCountsTab
#----------------------------------------------------------------------------------------------------
createCorrelationsTab <- function()
{
   fluidRow(column(3, selectInput(inputId = "variableChooser.1",
                        label = NULL,
                        choices = colnames(tbl))),
            column(3, selectInput(inputId = "variableChooser.2",
                        label = NULL,
                        choices = colnames(tbl))),
            plotOutput("correlationPlotter"),
            column=12)

} # createCorrelationsTab
#----------------------------------------------------------------------------------------------------
createDataTableTab <- function()
{
   fluidRow(selectInput(inputId = "tableRowChooser",
                        label = "Choose Row of interest:",
                        choices = rownames(mtx)),
            DT::dataTableOutput("atacSrmTable"),
            column=12)

} # createDataTableTab
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(
  h5("ATAC-seq vs. SRM: an example using only KLF1 & FLT1"),
  mainPanel(width=12,
            tabsetPanel(type = "tabs",
                        tabPanel("Plot Counts", createPlotCountsTab()),
                        tabPanel("Correlate Counts", createCorrelationsTab()),
                        tabPanel("Table", createDataTableTab()),
                        tabPanel("README", pre(includeHTML("README.html")))
                        )
    ) # mainPanel
  ) # fluidPage

#--------------------------------------------------------------------------------------------------------------
server <- function(input, output) {

  output$atacSrmTable = DT::renderDataTable(
                as.data.frame(mtx),
                options=list(scrollX=TRUE,
                            scrollY="300px",
                            dom='t',
                            ordering=FALSE,
                            paging=FALSE,
                            autowWidth=TRUE,
                            fixedColumns=list(leftColumns=1)
                            ))
   observeEvent(input$plotCorrelationButton, {
                   printf("click!")
                   })

   output$countsPlotter <- renderPlot({
      row.1 <- input$countsChooser; #"FLI1 srm count"
      #plot(mtx[row.1,])
      m <- aes_string(x="time", y=row.1) #FLI1.total.hits)
      p <- ggplot(data=tbl) + geom_point(mapping=m) + geom_smooth(mapping=m) + xlab("Day & replicate") +
              scale_x_continuous(breaks=1:9, labels=x.axis.labels)
      p
      })

   output$correlationPlotter <- renderPlot({
      vec.1 <- input$variableChooser.1
      vec.2 <- input$variableChooser.2
      #plot(tbl[, row.1], tbl[, row.2])
      correlation <- round(cor(tbl[, vec.1], tbl[, vec.2], method="spearman"), 4)
      grob1 <- grobTree(textGrob(paste("Spearman: ", correlation),
                                 x = 0.63, y = 0.97, hjust = 0, gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))

      ggplot(tbl, aes_string(x=vec.1, y=vec.2)) + geom_point() +
        geom_smooth(method=lm, se=FALSE) +
        annotation_custom(grob1) +
        theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),
              axis.line = element_line(color="black"), axis.line.x = element_line(color="black")) +
         theme_bw()
      })


}
#--------------------------------------------------------------------------------------------------------------

runApp(shinyApp(ui, server))
