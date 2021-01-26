#library(TrenaProjectErythropoiesis)
library(shiny)
library(shinydashboard)
#library(shinyModules)
library(igvWidget)
library(dataTableWidget)
library(msgBoxWidget)
library(RUnit)
library(ghdb)
library(R6)
library(later)
library(httr)
library(jsonlite)
#------------------------------------------------------------------------------------------------------------------------
ghdb <- GeneHancerDB()
mtx <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-16173x28.RData"))
mtx <- round(mtx, digits=2)
day0.mean <- apply(mtx[, 1:2], 1, mean)
day2.mean <- apply(mtx[, 3:4], 1, mean)
tbl.corces <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.corces.hg38.scoredByPaul.RData")) #
tbl.diffbind <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.diffBind.rory.hg38.day0-4.RData"))
  # tbl.atac.corces.day0-4.RData"))
# tpe <- TrenaProjectErythropoiesis();
# expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")
# checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
# mtx <- getExpressionMatrix(tpe, expected[1])
#------------------------------------------------------------------------------------------------------------------------
EarlyActivationApp <- R6Class("app",

   #------------------------------------------------------------
   private = list(tbl = NULL,
                  tbl.filtered = NULL,
                  igv = NULL,
                  dtw = NULL,
                  msgBox = NULL,
                  input = NULL,
                  output = NULL,
                  session = NULL
                  ),
   #------------------------------------------------------------
   public = list(

      #---------------------------------------------------------------------------------------
      initialize = function(){
         private$tbl = as.data.frame(mtx)
         private$igv = igvWidget$new("igv01",
                                     genome="hg38",
                                     locus="APOE",
                                     width="100%",
                                     height="800",
                                     border="1px solid purple; border-radius: 5px;")
         private$dtw = dataTableWidget$new(id="dtw", private$tbl,
                                           width="100%", height="1000px",
                                           border="1px blue solid; border-radius: 5px;",
                                           pageLength=25,
                                           lengthMenu=c(5,15,25,50))
         private$msgBox =  msgBoxWidget$new(id="box1", title="", boxWidth=200)
         },

      #---------------------------------------------------------------------------------------
      ui = function(){
        tabsetPanel(type="tabs",id="mainTabs", selected="geneExpressionTab",
           tabPanel(title="Introduction", value="introductionTab", includeHTML("intro.html")),
           tabPanel(title="Gene Expression", value="geneExpressionTab",
                  dashboardPage(
                      dashboardHeader(title="Early Erythropoiesis Gene Activation"),
                      dashboardSidebar(
                          sliderInput(inputId="day0Max", label="Mean Day 0 Max lfc", min=0, max=5, value=0, step=0.1),
                          sliderInput(inputId="day2Min", label="Mean Day 2 Min lfc", min=0, max=15, value=2.8, step=0.1)
                      ),
                      dashboardBody(private$dtw$ui())
                  ) # dashboardPage
               ),
           tabPanel(title="IGV", value="igvTab",
                   fluidPage(
                      div(private$msgBox$ui(), style="margin: 10px;"),
                      private$igv$ui()
                      )
                    )
           ) # tabsetPanel
         }, # ui

      #---------------------------------------------------------------------------------------
      server = function(input, output, session){
         private$dtw$server(input, output, session)
         private$msgBox$server(input, output, session)
         private$igv$server(input, output, session)
         private$input <- input;
         private$output <- output;
         private$session <- session;

         observe({
            day0.max <- as.numeric(input$day0Max)
            day2.min <- as.numeric(input$day2Min)
            goi.0 <- names(which(day0.mean <= day0.max))
            goi.2 <- names(which(day2.mean >= day2.min))
            goi <-intersect(goi.0, goi.2)
            if(length(goi) == 0)
                private$tbl.filtered <- data.frame()
            if(length(goi) > 0)
                private$tbl.filtered <- as.data.frame(mtx[goi,,drop=FALSE])
            private$dtw$setTable(private$tbl.filtered)
            })

         observe({
           indices <- private$dtw$tableSelection()
           req(indices)
           selected.genes <- rownames(private$tbl.filtered)[indices]
           gene.of.interest <- selected.genes[1]
           updateTabsetPanel(session, "mainTabs", selected="igvTab")
           later(function(){
               private$igv$setLocus(gene.of.interest)
               loc.list <- self$addGeneHancerTracks(gene.of.interest)
               printf("back from gh")
               if(length(loc.list) == 0){
                  uri <- sprintf("http://localhost:8000/geneLoc")
                  body.jsonString <- sprintf('%s', toJSON(list(gene=gene.of.interest,
                                                               genome="hg38", shoulder=100000)))
                  r <- POST(uri, body=body.jsonString)
                  loc.list <- fromJSON(content(r)[[1]])
                  loc.string <- with(loc.list, sprintf("%s:%d-%d", chrom, start, end))
                  showModal(modalDialog(
                     title="GeneHancer",
                     sprintf("No GeneHancer regions for %s, using tss +/- 100k: %s",
                             gene.of.interest, loc.string)),
                     session=private$session)
                  } # if no GeneHancer
               with(loc.list, printf("full loc, from gh or service: %dk", as.integer((end-start)/1000)))
               print(loc.list)
               with(loc.list, private$igv$setLocus(sprintf("%s:%d-%d", chrom, start, end)))
               private$msgBox$setText(paste(selected.genes, collapse=", "))
               self$addCorcesATACtracks(loc.list)
               }, 1)
           }) # observe tableSelection()
         }, # server

      #---------------------------------------------------------------------------------------
      addGeneHancerTracks=function(geneSymbol){
         tbl.gh.all <- retrieveEnhancersFromDatabase(ghdb, geneSymbol, tissues="all")
         printf("tbl.gh.all: %d rows", nrow(tbl.gh.all))
         coi <- c("chrom", "start", "end", "combinedscore")
         if(nrow(tbl.gh.all) == 0)
            return(list())
         private$igv$displayBedGraphTrack("GH", tbl.gh.all[, c("chrom", "start", "end", "combinedscore")],
                                          color="darkRed", autoscale=FALSE,
                                          min=0, max=50, trackHeight=30)
         chrom <- tbl.gh.all$chrom[1]
         start <- min(tbl.gh.all$start) - 10000
         end <- max(tbl.gh.all$end) + 10000
         return(list(chrom=chrom, start=start, end=end))
         }, # addGeneHancerTracks

      #---------------------------------------------------------------------------------------
      addCorcesATACtracks=function(loc){
        printf("======== entering addCorcesATACtracks")
        print(loc)
        tbl.corces.sub <- subset(tbl.corces, chrom==loc$chrom & start >= loc$start & end <= loc$end & score > 2)
        printf("   tbl.corces.sub: %d rows", nrow(tbl.corces.sub))
        tbl.diffbind.sub <-   subset(tbl.diffbind, chrom==loc$chrom & start >= loc$start & end <= loc$end)
        printf("   tbl.diffbind.sub: %d rows", nrow(tbl.diffbind.sub))
        if(nrow(tbl.corces.sub) > 0){
           printf("=== loose")
           tbl.track <- tbl.corces.sub[, c("chrom", "start", "end", "signed.score")]
           print(head(tbl.track))
           private$igv$displayBedGraphTrack("corces.loose", tbl.track,
                                            color="darkblue",
                                            autoscale=TRUE)
           save(tbl.track, file="tbl.corces.loose.RData")
           printf("======== after loose track displayed")
           } # tbl.corces.sub
        if(nrow(tbl.diffbind.sub) > 0){
           print(head(tbl.diffbind.sub))
           print(dim(tbl.diffbind.sub))
           tbl.track <- tbl.diffbind.sub[, c("chrom", "start", "end", "Fold")]
           printf("=== strict")
           save(tbl.track, file="tbl.corces.diffbind.RData")
           print(head(tbl.track))
           #tbl.track$score <- -log10(tbl.track$p.value)
           private$igv$displayBedGraphTrack("corces.strict", tbl.track[, c("chrom", "start", "end", "Fold")],
                                            color="darkgreen", autoscale=TRUE)
           } # tbl.diffbind.sub
        } # addCorcesATACtracks

     ) # public
  ) # app

#------------------------------------------------------------------------------------------------------------------------
# server <- function(input, output, session)
# {
#
#   igv$server(input, output, session)
#
#   tbl.current <- reactive({  # lazy, only runs when called
#      day0.max <- as.numeric(input$day0Max)
#      day2.min <- as.numeric(input$day2Min)
#      goi.0 <- names(which(day0.mean <= day0.max))
#      goi.2 <- names(which(day2.mean >= day2.min))
#      goi <-intersect(goi.0, goi.2)
#      if(length(goi) == 0)
#          tbl.filtered <- data.frame()
#      if(length(goi) > 0)
#          tbl.filtered <- as.data.frame(mtx[goi,,drop=FALSE])
#      tbl.filtered
#      })
#   geneOfInterest <- dataTableServer("table",
#                                     tbl=tbl.current,
#                                     selectionPolicy=reactive("single"),
#                                     wrapLongTextInCells=reactive(FALSE),
#                                     )
#    messageBoxServer("messageBox.igv", newContent=private$igv$trackClickResult)
#
# } # server

#      geneOfInterest <- dataTableServer("table",
#                                    tbl=tbl.current,
#                                    selectionPolicy=reactive("single"),
#                                    wrapLongTextInCells=reactive(FALSE),
#                                    )
#  roi <- reactive({
#     goi <- geneOfInterest()
#     req(goi)
#     calculateGenomicRegion(goi, 1000)
#     })
#
#  observe({
#     goi <- geneOfInterest()
#     req(goi)
#     printf("---- goi:")
#     print(goi)
#     #req(goi)
#     printf("geneOfInterest: %s", goi)
#     gh.roi <- addGeneHancerTracks(session, goi)
#     if(length(gh.roi) > 0)
#        geneRegion <- gh.roi
#     else
#        geneRegion <- roi()
#     printf("==== roi, using gh if available")
#     print(geneRegion)
#     chromLoc <- with(geneRegion, sprintf("%s:%d-%d", chrom, start, end))
#     printf("chromLoc: %s", geneRegion$regionString)
#     igv$setLocus(geneRegion$regionString)
#     #addCorcesATACtracks(session, geneRegion)
#     })
#
#  messageBoxServer("tableSelection", newContent=geneOfInterest)

# } # server
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

  printf("--- leaving addGeneHancerTracks")
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
      loadBedGraphTrack(session, time, tbl.atac, color="random", autoscale=FALSE, min=0, max=50)
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
app <- EarlyActivationApp$new()
#shinyApp(ui, server)
runApp(shinyApp(app$ui, app$server), port=9999)
