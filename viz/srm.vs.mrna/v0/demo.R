library(shiny)
library(devtools)
library(r2d3)

ui <- fluidPage(
   tags$head(tags$style("#d3{height:90vh !important;}")),
   titlePanel("mRNA vs protein counts"),
   sidebarLayout(
      sidebarPanel(
         sliderInput("dayNumberSlider", label = "Day:", min = 0, max = 16, value = 0, step = 1, round=TRUE),
         width=2
         ),
      mainPanel(
         d3Output("d3", height="800px"),
         width=10
      )
   )
)

server <- function(input, output) {
  output$d3 <- renderD3({
     count <- input$dayNumberSlider
     xMax <- 16
     yMax <- 16^2
     x <- 0
     y <- 0
     if(count > 0){
        x <- seq_len(count)
        y <- x^2
        }
     data <- list(values=lapply(seq_len(count), function(i) return(list(x=x[i], y=y[i]))),
                  xMax=xMax, yMax=yMax)
     r2d3(data, script = "demo-reSize.js")
    })
}

app <- shinyApp(ui = ui, server = server)
