#library(TrenaProjectErythropoiesis)
library(shiny)
library(shinyModules)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
mtx <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-16173x28.RData"))
mtx <- round(mtx, digits=2)
day0.mean <- apply(mtx[, 1:2], 1, mean)
day2.mean <- apply(mtx[, 3:4], 1, mean)

# tpe <- TrenaProjectErythropoiesis();
# expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")
# checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
# mtx <- getExpressionMatrix(tpe, expected[1])
#------------------------------------------------------------------------------------------------------------------------
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title="Early Erythropoiesis Gene Activation"),
  dashboardSidebar(
      sliderInput(inputId="day0Max", label="Mean Day 0 Max lfc", min=0, max=5, value=0, step=0.1),
      sliderInput(inputId="day2Min", label="Mean Day 2 Min lfc", min=0, max=15, value=3, step=0.1)
      ),
  dashboardBody(
    div(
       dataTableUI("table"),
       style="margin: 20px; padding: 10px; border: 2px solid black; border-radius: 10px;"
       ), # div
   div(messageBoxUI(id="selectResultsDisplay", title="selection", boxWidth=800, boxHeight=50,
                    fontSize=20),
       style="margin-left: 30px")
    ) # dashboardBody
  ) # dashboardPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session)
{

  tbl.current <- reactive({  # lazy, only runs when called
     day0.max <- as.numeric(input$day0Max)
     day2.min <- as.numeric(input$day2Min)
     goi.0 <- names(which(day0.mean <= day0.max))
     goi.2 <- names(which(day2.mean >= day2.min))
     goi <-intersect(goi.0, goi.2)
     if(length(goi) == 0)
         tbl.filtered <- data.frame()
     if(length(goi) > 0)
         tbl.filtered <- as.data.frame(mtx[goi,,drop=FALSE])
     tbl.filtered
     })

  selectedRow <- dataTableServer("table",
                                 tbl=tbl.current,
                                 selectionPolicy=reactive("single"),
                                 wrapLongTextInCells=reactive(FALSE),
                                 )

  observe({
     req(selectedRow())
     printf("selectedRow changed: %s", selectedRow())
     #x <- selectedRow()
     })

  messageBoxServer("selectResultsDisplay", newContent=selectedRow)

} # server
#------------------------------------------------------------------------------------------------------------------------
#shinyApp(ui, server)
runApp(shinyApp(ui, server), port=9999)
