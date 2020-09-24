#library(TrenaProjectErythropoiesis)
library(shiny)
library(shinyModules)
library(RUnit)
library(ghdb)
ghdb <- GeneHancerDB()
#------------------------------------------------------------------------------------------------------------------------
mtx <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-16173x28.RData"))
mtx <- round(mtx, digits=2)
day0.mean <- apply(mtx[, 1:2], 1, mean)
day2.mean <- apply(mtx[, 3:4], 1, mean)
tbl.corces <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atac.corces.day0-4.RData"))
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
      sliderInput(inputId="day2Min", label="Mean Day 2 Min lfc", min=0, max=15, value=4.3, step=0.1),
      div(messageBoxUI(id="tableSelection", title="table selection", boxWidth=100, boxHeight=30, fontSize=12),
         style="margin-left: 20px"),
      div(radioButtons("regionOfInterest", "ATAC-seq region",
                       choices=c("TSS +/-5kb", "Enhancer Extent", "Enhancer Intersect"),
                       selected="TSS +/-5kb"),
       style="display: inline-block;vertical-align:top; margin-left: 20px; width: 200px;")

      ),
  dashboardBody(
    div(
       dataTableUI("table"),
       style="margin: 20px; padding: 10px; border: 2px solid black; border-radius: 10px;"
       ), # div
    div(igvUI("igv"),
       style="margin: 10px; margin-bottom: 5px; padding: 10px; border: 3px solid gray; border-radius: 10px; background-color: white;"),
    div(messageBoxUI(id="messageBox.igv", title="igv selection", boxWidth= 600), style="margin-left: 100px;")
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

  geneOfInterest <- dataTableServer("table",
                                    tbl=tbl.current,
                                    selectionPolicy=reactive("single"),
                                    wrapLongTextInCells=reactive(FALSE),
                                    )
  roi <- reactive({
     goi <- geneOfInterest()
     req(goi)
     calculateGenomicRegion(goi, 1000)
     })

  observe({
     goi <- geneOfInterest()
     req(goi)
     printf("---- goi:")
     print(goi)
     #req(goi)
     printf("geneOfInterest: %s", goi)
     gh.roi <- addGeneHancerTracks(session, goi)
     if(length(gh.roi) > 0)
        geneRegion <- gh.roi
     else
        geneRegion <- roi()
     printf("==== roi, using gh if available")
     print(geneRegion)
     chromLoc <- with(geneRegion, sprintf("%s:%d-%d", chrom, start, end))
     printf("chromLoc: %s", chromLoc)
     showGenomicRegion(session, chromLoc)
     addCorcesATACtracks(session, geneRegion)
     })

  messageBoxServer("tableSelection", newContent=geneOfInterest) #reactive(roi()$regionString))

  igvSelection <- igvServer("igv",
                            genome="hg38",
                            geneModelDisplayMode="COLLAPSED",
                            locus="APOE",
                            height=400)

  messageBoxServer("messageBox.igv", newContent=igvSelection)

} # server
#------------------------------------------------------------------------------------------------------------------------
addGeneHancerTracks <- function(session, geneSymbol)
{
  #tbl.gh.cd34 <- retrieveEnhancersFromDatabase(ghdb, geneSymbol, tissues="Common myeloid progenitor CD34+")
  #printf("tbl.gh.cd34: %d rows", nrow(tbl.gh.cd34))

  tbl.gh.all <- retrieveEnhancersFromDatabase(ghdb, geneSymbol, tissues="all")
  printf("tbl.gh.all: %d rows", nrow(tbl.gh.all))

  coi <- c("chrom", "start", "end", "combinedscore")
  #if(nrow(tbl.gh.cd34) > 0)
  #  loadBedGraphTrack(session, "GH.CD34+", tbl.gh.cd34[, coi], color="random", autoscale=FALSE, min=0, max=50)

  if(nrow(tbl.gh.all) == 0)
     return()

  loadBedGraphTrack(session, "GH.all", tbl.gh.all[, coi], color="black", autoscale=FALSE, min=0, max=50)

  chrom <- tbl.gh.all$chrom[1]
  start <- min(tbl.gh.all$start) - 1000
  end <- max(tbl.gh.all$end) + 1000

  return(list(chrom=chrom, start=start, end=end))

} # addGeneHancerTracks
#------------------------------------------------------------------------------------------------------------------------
addCorcesATACtracks <- function(session, roi)
{
   printf("adding corces atac-seq in region")
   print(roi)
   tbl.atac <- subset(tbl.corces, chrom==roi$chrom & start >= roi$start & end <= roi$end)
   print(head(tbl.atac))
   print(dim(tbl.atac))
   for(time in c("day0.1", "day0.2", "day0.3", "day0.4", "day2.1", "day2.2")){
      tbl <- tbl.atac[, c("chrom", "start", "end", time)]
      loadBedGraphTrack(session, time, tbl, color="random", autoscale=FALSE, min=0, max=50)
      } # for time

} # addCorcesATACtracks
#------------------------------------------------------------------------------------------------------------------------
calculateGenomicRegion <- function(geneSymbol, shoulder)
{
   require(jsonlite)
   require(httr)
   uri <- sprintf("http://localhost:8000/geneLoc")
   body.jsonString <- sprintf('%s', toJSON(list(gene=geneSymbol, genome="hg38", shoulder=5000)))
   r <- POST(uri, body=body.jsonString)
   region.hg38 <- fromJSON(content(r)[[1]])
   region.hg38$regionString <- with(region.hg38, sprintf("%s:%d-%d", chrom, start, end))
   printf("=== region.hg38 back from geneLoc service")
   print(region.hg38)
   return(region.hg38)

} # calculateGenomicRegion
#------------------------------------------------------------------------------------------------------------------------
#shinyApp(ui, server)
runApp(shinyApp(ui, server), port=9999)
