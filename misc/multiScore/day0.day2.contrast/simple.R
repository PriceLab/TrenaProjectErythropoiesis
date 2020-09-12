library(shiny)
library(RUnit)
library(shinyModules)
#------------------------------------------------------------------------------------------------------------------------
mtx <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-16173x28.RData"))
mtx <- round(mtx, digits=2)
#------------------------------------------------------------------------------------------------------------------------
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title="Early Erythropoiesis Gene Activation"),
  dashboardSidebar(
      sliderInput(inputId="day0Max", label="Mean Day 0 Max lfc", min=0, max=5, value=0, step=0.1),
      sliderInput(inputId="day2Min", label="Mean Day 2 Min lfc", min=0, max=15, value=3, step=0.1)
      ),
  dashboardBody(
   div(messageBoxUI(id="selectResultsDisplay", title="selection", boxWidth=800, boxHeight=50, fontSize=20),
       style="margin-left: 30px")
    ) # dashboardBody
  ) # dashboardPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session)
{
    sliderValues <- reactive({
        list(day0=input$day0Max, day2=input$day2Min)
        })

#    msg <- reactive({
#        printf("--- entering reactive");
#        day0.max <- input$day0Max
#        day2.min <- input$day2Min
#        tbl.filtered <- mtcars
#        sprintf("msg: %s-%s", day0.max, day2.min)
        # day0.max <- as.numeric(input$day0Max)
        # day2.min <- 3 #as.numeric(input$day2Min)
        # printf("day0.max: %d", day0.max)
        # printf("day2.min: %d", day2.min)
        # goi.0 <- names(which(day0.mean <= day0.max))
        # goi.2 <- names(which(day2.mean >= day2.min))
        # goi <-intersect(goi.0, goi.2)
        # if(length(goi) == 0)
        #     tbl.filtered <- data.frame()
        # if(length(goi) > 0)
        #     tbl.filtered <- as.data.frame(mtx[goi,,drop=FALSE])
        # printf("filtered line count: %d", nrow(tbl.filtered))
        # return(tbl.filtered)
#        })

    rvs <- reactiveValues(tbl=NULL)

    observe({
        printf("--- entering tbl reactive");
        day0.max <- input$day0Max
        day2.min <- input$day2Min
        printf("msg: %f-%f", day0.max, day2.min)

        rows.to.sample <- sample(1:10, 1)
        mtcars.sub <- mtcars[sample(1:20, rows.to.sample),]
        printf("rows: %d", nrow(mtcars.sub))
        rvs$tbl <- mtcars.sub
        # day0.max <- as.numeric(input$day0Max)
        # day2.min <- 3 #as.numeric(input$day2Min)
        # printf("day0.max: %d", day0.max)
        # printf("day2.min: %d", day2.min)
        # goi.0 <- names(which(day0.mean <= day0.max))
        # goi.2 <- names(which(day2.mean >= day2.min))
        # goi <-intersect(goi.0, goi.2)
        # if(length(goi) == 0)
        #     tbl.filtered <- data.frame()
        # if(length(goi) > 0)
        #     tbl.filtered <- as.data.frame(mtx[goi,,drop=FALSE])
        # printf("filtered line count: %d", nrow(tbl.filtered))
        # return(tbl.filtered)
        })


  messageBoxServer("selectResultsDisplay", newContent=reactive(nrow(rvs$tbl)))
  #messageBoxServer("selectResultsDisplay", newContent=tbl) #reactive("foo"))

} # server
#------------------------------------------------------------------------------------------------------------------------
#shinyApp(ui, server)
runApp(shinyApp(ui, server), port=9999)
