library(shiny)
library(devtools)
library(r2d3)
#------------------------------------------------------------------------------------------------------------------------
mtx.srm <- get(load("../SRMforPublication20190614.RData"))
mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
#------------------------------------------------------------------------------------------------------------------------
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
     day <- as.integer(input$dayNumberSlider)
     xMax <- 1.0
     yMax <- 1.0
     tf <- "GATA1"
     rna <- as.numeric(mtx.rna[tf, c(1,3,5,7,9,11,13,15,17,19,21,23,25)])
     rna <- rna/max(rna)
     srm <- as.numeric(mtx.srm[tf,])
     srm <- srm/max(srm)
     stopifnot(length(rna) == length(srm))
     data <- list(values=lapply(seq_len(day+1), function(i) return(list(x=rna[i], y=srm[i]))),
                  xMax=xMax, yMax=yMax)
     #print(data)
     r2d3(data, script = "demo-reSize.js")
    })
}

app <- shinyApp(ui = ui, server = server)
