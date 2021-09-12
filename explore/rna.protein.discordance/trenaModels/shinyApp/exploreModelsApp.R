library(shiny)
library(R6)
library(DataTableWidget)
library(HeatmapWidget)
library(MsgBoxWidget)
library(GOEnrichmentWidget)
library(shinyWidgets)   # for pickerInput
source("toHeatmapMatrix.R")
library(rsconnect)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
DemoApp = R6Class("DemoApp",

    #--------------------------------------------------------------------------------
    private = list(dtw = NULL,
                   tfGeneCountReadout = NULL,
                   heatmapWidget = NULL,
                   goWidget = NULL,
                   tbl.models= NULL,
                   tbl.sub = NULL,
                   tbl.currentSubset = NULL,
                      # calculated afresh on button click, from current filters:
                   heatmap.mtx = NULL,
                   eventMessage = NULL,
                   input = NULL,
                   output = NULL,
                   session = NULL

                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            printf("initializing demo")
            filename <- "tbl.trena-100-targets.RData"
            printf("--- loading %s", filename)
            tbl.tmp <- get(load(filename))
            printf("--- load complete")
            coi <- c("target", "gene", "class", "spearmanCoeff", "pearsonCoeff", "betaLasso",
                     "betaRidge", "rfScore", "xgboost")
            tbl.tmp <- tbl.tmp[, coi]
            colnames(tbl.tmp)[2] <- "TF"
            rownames(tbl.tmp) <- NULL
            tbl.tmp$pearsonCoeff <- round(tbl.tmp$pearsonCoeff, 2)
            tbl.tmp$betaLasso <- round(tbl.tmp$betaLasso, 2)
            tbl.tmp$betaRidge <- round(tbl.tmp$betaRidge, 2)
            tbl.tmp$xgboost <- round(tbl.tmp$xgboost, 2)
            private$tbl.models <- tbl.tmp
            private$tbl.sub <- private$tbl.models
            private$tbl.currentSubset <- private$tbl.models
            rownames(private$tbl.models) <- NULL
            private$dtw = dataTableWidget$new(id="dtw", private$tbl.models,
                                              pageLength=15,
                                              lengthMenu=c(4, 10, 15, 20, 25, 30, 50),
                                              #columnWidths=c(30, 30, 10, 30, 30),
                                              width="1600px", height="1000px")
            private$goWidget <- GOEnrichmentWidget$new(id="goWidget.1", title="GO Enrichment",
                                                      geneSymbols=sort(unique(private$tbl.models$target)))
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              titlePanel("100 trena models, mRNA only (so far)"),
              tabsetPanel(type="tabs",id="navigationTabSet", selected="trenaModelsTab",
                 # tabPanel(title="Introduction", value="introTab", includeHTML("intro.html")),
                 tabPanel(title="trena models", value="trenaModelsTab",
                      sidebarLayout(
                        sidebarPanel(
                            #private$tfGeneCountReadout$ui(),
                            h4(textOutput("tfGeneCountReadout")),
                            actionButton("displayHeatmapButton", "Selection to Heatmap"),
                            actionButton("sendGenesToGoWidgetButton", "Selected Genes to GO"),
                            br(), br(),
                            #actionButton("filterTableButton", "Filter by Current Settings"),
                            #br(), br(),
                            pickerInput(
                                inputId = "targetGenePicker",
                                label = "Select a Target Gene",
                                choices = sort(unique(private$tbl.models$target)),
                                selected = sort(unique(private$tbl.models$target)),
                                multiple = TRUE,
                                options = pickerOptions(
                                   actionsBox = TRUE,
                                   liveSearch=TRUE,
                                   `selected-text-format`="count > 4"
                                   )
                            ),
                            pickerInput(
                                inputId = "tfPicker",
                                label = "Select candidate Regulators",
                                choices = sort(unique(private$tbl.models$TF)),
                                selected = "STAU1", #sort(unique(private$tbl.models$TF)),
                                multiple = TRUE,
                                options = pickerOptions(
                                   actionsBox = TRUE,
                                   liveSearch=TRUE,
                                   `selected-text-format`="count > 4"
                                   )
                               ),

                            sliderInput("spearmanSlider", "Spearman Correlation", min = 0, max = 1.0,
                                        value=0.85),
                            radioButtons(
                                inputId="spearmanAbsoluteFilterRadioButtons",
                                label="abs(spearman)?",
                                choices = c("yes", "no"),
                                selected = "yes",
                                inline = TRUE,
                                width = NULL,
                                choiceNames = NULL,
                                choiceValues = NULL),
                            width=3), # sidebarPanel
                        mainPanel(
                          div(private$dtw$ui(), style="margin:20px;"),
                          width=9
                          ) # mainPanel
                      ) # sidebarLayout for openChromatin
                      ), # tabPanel
                 tabPanel(title="Heatmap", value="heatmapTab",
                    sidebarLayout(
                        sidebarPanel(
                            h4("choosers"),
                            width=3
                            ),
                        mainPanel(
                            div(id="heatmapContainer", div(id="heatmapWrapper")),
                            width=9
                            )
                        ) # sidebarLayout
                    ), # heatmap tabPanel
                 tabPanel(title="GO", value="GOtab",
                    private$goWidget$ui()
                    )
                 ) # tabsetPanel
            )},

        #------------------------------------------------------------
        server = function(input, output, session){

           private$input <- input
           private$output <- output
           private$session <- session

           private$dtw$server(input, output, session)
           private$goWidget$server(input, output, session)
           output$tfGeneCountReadout <- renderText(sprintf("%4d targetGene/s, %4d tf/s",
                                                           length(unique(private$tbl.models$target)),
                                                           length(unique(private$tbl.models$tf))))


           #---------------------------------------------------------------
           updateTable <- function(){
              private$tbl.sub <- private$tbl.models
              tfs <- input$tfPicker
              if(!is.null(tfs)){
                  printf("--- tfs is null")
                  private$tbl.sub <- subset(private$tbl.sub, TF %in% tfs)
                  }
              printf("--- tfs: %d", length(tfs))
              spearmanThreshold <- input$spearmanSlider
              spearmanAbsolute <- input$spearmanAbsoluteFilterRadioButtons
              if(spearmanAbsolute == "yes"){
                 private$tbl.sub <- subset(private$tbl.sub, abs(spearmanCoeff) >= spearmanThreshold)
              } else {
                 private$tbl.sub <- subset(private$tbl.sub, spearmanCoeff >= spearmanThreshold)
                 }
              targetGenes <- input$targetGenePicker
              if(!is.null(targetGenes)){
                private$tbl.sub <- subset(private$tbl.sub, target %in% targetGenes)
                printf("targetGenes is null")
                }

              printf("--- targetGenes: %d", length(targetGenes))
              output$tfGeneCountReadout <- renderText(sprintf("%4d targetGene/s, %4d tf/s",
                                                              length(unique(private$tbl.sub$target)),
                                                              length(unique(private$tbl.sub$TF))))

              rownames(private$tbl.sub) <- NULL
              private$dtw$setTable(private$tbl.sub)
              } # updateTable

           observe(updateTable())

           observeEvent(input$filterTableButton, ignoreInit=TRUE, {
              printf("--- filterTableButton clicked")
              updateTable()
              }) # observeEvent: filterTableButton



           #---------------------------------------------------------------------------------------
             # change the label of the two open chromatin reads, min or max, reflecting user's choice
           observeEvent(input$openChromatinPreference, ignoreInit=TRUE, {
              day0.closedChromatin.day2.open <- input$openChromatinPreference == "yes"
              slider.label.day0 <- "Day 0 open chromatin max reads:"
              slider.label.day2 <- "Day 2 open chromatin min reads:"
              if(!day0.closedChromatin.day2.open){
                 slider.label.day0 <- "Day 0 open chromatin min reads:"
                 slider.label.day2 <- "Day 2 open chromatin max reads:"
                 }
              #updateSliderInput(private$session, "day0.oc.slider", label = slider.label.day0)
              #updateSliderInput(private$session, "day2.oc.slider", label = slider.label.day2)
              })

           #---------------------------------------------------------------------------------------
           observeEvent(input$displayHeatmapButton, ignoreInit=TRUE, {
              tf.count <- length(unique(private$tbl.sub$TF))
              targetGene.count <- length(unique(private$tbl.sub$target))
              printf("--- displayHeatmapButton, targetGenes %d, tfs %d", targetGene.count, tf.count)
              tbl.tmp <- private$tbl.sub
              colnames(tbl.tmp)[1] <- "targetGene"   # match heatmap's expectations
              colnames(tbl.tmp)[2] <- "tf"           # match heatmap's expectations
              printf("--- about to clsasify for heatmap")
              print(head(private$tbl.sub))
              print(head(tbl.tmp))
              tbl.class <- private$tbl.sub[, "class", drop=FALSE]
              x <- HeatMapMatrixFromTMSTable$new(tbl.tmp)
              mtx.tmp <- x$calculate()
              printf(" xtab matrix, before, %d x %d", nrow(mtx.tmp), ncol(mtx.tmp))
              if(nrow(mtx.tmp) == 1)
                  mtx.tmp <- rbind(mtx.tmp, mtx.tmp)
              if(ncol(mtx.tmp) == 2){
                 printf("====== doubling columns")
                 mtx.tmp <- cbind(mtx.tmp, mtx.tmp)
                 }
              printf(" xtab matrix, after, %d x %d", nrow(mtx.tmp), ncol(mtx.tmp))
              private$heatmap.mtx <- mtx.tmp
              mtx.xtab <- private$heatmap.mtx
              private$heatmapWidget = HeatmapWidget$new(id="box1",
                                                        title="Coordinate Regulation",
                                                        private$heatmap.mtx,
                                                        rowClusters=5,
                                                        colClusters=5,
                                                        rowTitle="Gene",
                                                        columnTitle="TF",
                                                        colGroups=tbl.class,
                                                        width="100%",
                                                        height=700)
              removeUI(selector="#heatmapWrapper", immediate=TRUE)
              updateTabsetPanel(session, "navigationTabSet", selected="heatmapTab")
              insertUI(selector="#heatmapContainer", where="beforeEnd",
                       div(id="heatmapWrapper"), immediate=TRUE)
              insertUI(selector="#heatmapWrapper", where="beforeEnd", private$heatmapWidget$ui(), immediate=TRUE)
              private$heatmapWidget$server(input, output, session)
              })

           #---------------------------------------------------------------------------------------
           observeEvent(input$sendGenesToGoWidgetButton, ignoreInit=TRUE, {
               goi <- unique(private$tbl.sub$target) # [c(1:10, 1000:1010)])
               count <- length(goi)
               printf("genes for enrichment: %d", count)
               printf("tbl.sub: %d, %d", nrow(private$tbl.sub), ncol(private$tbl.sub))
               updateTabsetPanel(session, "navigationTabSet", selected="GOtab")
               private$goWidget$setGenes(goi)
              })


           #---------------------------------------------------------------------------------------
           #observeEvent(input$heatmap_click, ignoreInit=TRUE, {
           #    x <- input$heatmap_click
           #    click.info <- private$heatmapWidget$clickEventAsPosition(x)
           #    if(length(click.info) == 0) return()
           #    tbl.clickInfo <- as.data.frame(click.info)
           #    row <- tbl.clickInfo$row_index[1]
           #    col <- tbl.clickInfo$column_index[1]
           #      # update the reactiveVal:
           #    msg <- sprintf("%d %s binding sites for %s",
           #                   private$heatmap.mtx[row, col],
           #                   colnames(private$heatmap.mtx)[col],
           #                   rownames(private$heatmap.mtx)[row])
           #    output$heatmapClickReadout <- renderUI(HTML(msg))
           #    })

           #---------------------------------------------------------------------------------------
           #observeEvent(input$heatmap_brush, ignoreInit=TRUE, {
           #   #printf("=== tmsNavigator sees heatmap brush")
           #   #x <- input$heatmap_brush
           #   #print(x)
           #   })

            # })
           } # server
       ) # public
    ) # class
#--------------------------------------------------------------------------------
deploy <-function()
{
   repos <- options("repos")[[1]]
   stopifnot(sort(names(repos)) == c("BioCann", "BioCsoft", "CRAN"))
   stopifnot(repos$BioCann=="https://bioconductor.org/packages/3.13/data/annotation")
   stopifnot(repos$BioCsoft=="https://bioconductor.org/packages/3.13/bioc")
   stopifnot(repos$CRAN=="https://cran.microsoft.com")

   require(devtools)
   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
   install_github("PriceLab/BioShiny/MsgBoxWidget",      force=TRUE)
   install_github("PriceLab/BioShiny/HeatmapWidget",      force=TRUE)
   install_github("PriceLab/BioShiny/DataTableWidget",    force=TRUE)
   install_github("PriceLab/BioShiny/GOEnrichmentWidget", force=TRUE)
   install_github("dreamRs/shinyWidgets", force=TRUE)
   install_github("dreamRs/shinybusy", force=TRUE)
   #install_github("paul-shannon/igvShiny", force=TRUE)
   #install_github("PriceLab/BioShiny/igvWidget", force=TRUE)
   #install_github("PriceLab/BioShiny/GenomeTracksWidget", force=TRUE)


   require(rsconnect)
   deployApp(account="hoodlab",
             appName="trenaHematopoiesis",
             appTitle="Trena Hematopoiesis",
             appFiles=c("exploreModelsApp.R",
                        "toHeatmapMatrix.R",
                        "tbl.trena-100-targets.RData",
                        "intro.html"
                        ),
             appPrimaryDoc="exploreModelsApp.R",
             forceUpdate=TRUE
             )

} # deploy
#------------------------------------------------------------------------------------------------------------------------
app <- DemoApp$new()
if(grepl("hagfish", Sys.info()[["nodename"]]) & !interactive()){
   runApp(shinyApp(app$ui(), app$server), port=1119)
   } else {
   shinyApp(app$ui(), app$server)
   }
