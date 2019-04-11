library(RCyjs)
library(later)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("parseChromLocString"))
   source("~/github/trena/R/utils.R")

if(!exists("tbl.geneInfo"))
   tbl.geneInfo <- get(load((system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))))

if(!exists("mtx")){
   mtx.raw <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
   mtx.asinh <- asinh(mtx.raw)
   mtx.asinh.norm <- t(apply(mtx.asinh, 1, function(row) {
                                row.min <- min(row)
                                row.max <- max(row)
                                new.values <- (row-row.min)/(row.max - row.min)
                                new.values}))
   mtx.raw.norm <-  t(apply(mtx.raw, 1, function(row) {
                               row.min <- min(row)
                               row.max <- max(row)
                               new.values <- (row-row.min)/(row.max - row.min)
                               new.values}))
   }

#------------------------------------------------------------------------------------------------------------------------
required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                           "motifScore", "motifRelativeScore", "match",
                                           "distance.from.tss", "tf")
#------------------------------------------------------------------------------------------------------------------------
geneRegulatoryModelToGraph <- function (target.gene, tbl.gm, tbl.reg)
{

   required.geneModelColumnNames <- c("tf", "betaLasso", "lassoPValue", "pearsonCoeff", "rfScore", "betaRidge",
                                     "spearmanCoeff", "bindingSites",  "sample")

   stopifnot(all(required.geneModelColumnNames %in% colnames(tbl.gm)))
   stopifnot(all(required.regulatoryRegionsColumnNames %in% colnames(tbl.reg)))

   printf("genes: %d, %d occurences of %d motifs", length(tbl.gm$tf), length(tbl.reg$motifName),
          length(unique(tbl.reg$motifName)))

   g <- graphNEL(edgemode = "directed")

   nodeDataDefaults(g, attr = "type") <- "undefined"             # targetGene, tf, footprint
   nodeDataDefaults(g, attr = "label") <- "default node label"
   nodeDataDefaults(g, attr = "distance") <- 0
   nodeDataDefaults(g, attr = "pearson") <- 0
   nodeDataDefaults(g, attr = "motif") <- ""
   nodeDataDefaults(g, attr = "xPos") <- 0
   nodeDataDefaults(g, attr = "yPos") <- 0
   nodeDataDefaults(g, attr = "expression") <- 0
   nodeDataDefaults(g, attr = "strength") <- 0


   edgeDataDefaults(g, attr = "edgeType") <- "undefined"

   tfs <- sort(unique(c(tbl.gm$tf, tbl.reg$tf)))

   regRegions.names <- unlist(lapply(1:nrow(tbl.reg), function(i){
       distance.from.tss <- tbl.reg$distance.from.tss[i]
       region.size <- nchar(tbl.reg$match[i])
       motif.name <- tbl.reg$motifName[i]
       if(distance.from.tss < 0)
          sprintf("%s.fp.downstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
        else
          sprintf("%s.fp.upstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
        }))

   tbl.reg$regionName <- regRegions.names
   all.nodes <- unique(c(target.gene, tfs, regRegions.names))
   g <- addNode(all.nodes, g)

   browser()

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, regRegions.names, "type")  <- "regulatoryRegion"
   nodeData(g, all.nodes, "label")  <- all.nodes
   nodeData(g, regRegions.names, "label") <- tbl.reg$motifName
   nodeData(g, regRegions.names, "distance") <- tbl.reg$distance.from.tss
   nodeData(g, regRegions.names, "motif") <- tbl.reg$motifName

   nodeData(g, tfs, "pearson") <- tbl.gm$pearsonCoeff
   nodeData(g, tfs, "rfScore") <- tbl.gm$rfScore

   g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
   edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"

   g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
   edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"

   g

} # geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
buildMultiModelGraph <- function (models, targetGene)
{

    g <- graphNEL(edgemode = "directed")
    model.names <- names(models)

    node.attribute.specs <- list(type="undefined",
                                 label="default node label",
                                 distance=0,
                                 pearson=0,
                                 randomForest=0,
                                 pcaMax=0,
                                 concordance=0,
                                 betaLasso=0,
                                 motif="",
                                 xPos=0,
                                 yPos=0)
    edge.attribute.spec <- list(edgeType="undefined")
    attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes

      # create current version of these attributes, and then
      # per-model versions, which get mapped to current
      # in response to user's interactive choice on the cyjs user interface
      # the "current version" is, e.g., "distance".
      # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
      # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
      # attribute, upon which the cyjs style is defined

    #for(class.name in attribute.classes){
       #class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
       #if(nchar(class.name) > 0)
       #   class.name.prefix <- sprintf("%s.", class.name)
       #noa.names.without.prefix <- names(node.attribute.specs)
       noa.names <- names(node.attribute.specs)
       #noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
       #noa.count <- length(node.attribute.specs)
       for(noa.name in noa.names){
          printf("--- noa.name: %s", noa.name)
          #browser()
          nodeDataDefaults(g, attr=noa.name) <- node.attribute.specs[noa.name]
          }
       #for(i in 1:noa.count){
       #   nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
       #   }
       #} # for class

    edgeDataDefaults(g, attr = "edgeType") <- "undefined"


    all.tfs <- c()
    all.regulatoryRegions <- c()

    for(model in models){  # collect all the tf and regulatory region nodes
       tbl.model <- model$model
       all.tfs <- unique(c(all.tfs, tbl.model$tf))
       tbl.reg <- model$reg
       regRegions.with.distance.from.tss <- sprintf("%s:%d", tbl.reg$motifName, tbl.reg$distance.from.tss)
       all.regulatoryRegions <- unique(c(all.regulatoryRegions, regRegions.with.distance.from.tss))
       } # for model

    printf("total tfs: %d   total regs: %d", length(all.tfs), length(all.regulatoryRegions))

    all.nodes <- unique(c(targetGene, all.tfs, all.regulatoryRegions))
    g <- addNode(all.nodes, g)

    nodeData(g, targetGene, "type") <- "targetGene"
    nodeData(g, all.tfs, "type")         <- "TF"
    nodeData(g, all.regulatoryRegions, "type")  <- "regulatoryRegion"
    nodeData(g, all.nodes, "label")  <- all.nodes

      # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes

   xyz <- "buildMultiModelGraph, about to interate over models"

    for(model in models){
       tbl.model <- model$model
       sample <- tbl.model$sample[1]
       tbl.reg <- model$reg
       printf("==== %s: %d rows", unique(tbl.model$sample), nrow(tbl.model))
       tfs <- tbl.reg$tf
       regRegions <- sprintf("%s:%d", tbl.reg$motifName, tbl.reg$distance.from.tss)
       printf("--- working stage %s", sample)
       #browser()
       xyze <- "addEdge"
       suppressWarnings(g <- addEdge(tfs, regRegions, g))
       edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
       suppressWarnings(g <- addEdge(regRegions, targetGene, g))
       edgeData(g, regRegions, targetGene, "edgeType") <- "regulatorySiteFor"
       nodeData(g, regRegions, "label") <- regRegions
       nodeData(g, regRegions, "distance") <- tbl.reg$distance.from.tss
       nodeData(g, regRegions, "motif") <- tbl.reg$motifName
       } # for model

      # now copy in the first model's tf node data

    tbl.model <- models[[1]]$model
    nodeData(g, tbl.model$tf, attr="randomForest") <- tbl.model$rfScore
    nodeData(g, tbl.model$tf, attr="pearson") <- tbl.model$pearsonCoeff

     # now copy in each of the model's tf node data in turn
    #model.names <- names(models)
    #for(model.name in model.names){
    #   tbl.model <- models[[model.name]]$model
    #   noa.name <- sprintf("%s.%s", model.name, "randomForest")
    #   nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$rfScore
    #   noa.name <- sprintf("%s.%s", model.name, "pearson")
    #    nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$pearsonCoeff
    #  } # for model.name

    g

} # buildMultiModelGraph
#------------------------------------------------------------------------------------------------------------------------
addGeneModelLayout <-  function (g, xPos.span=1500)
{
    printf("--- addGeneModelLayout")
    all.distances <- sort(unique(unlist(nodeData(g, attr='distance'), use.names=FALSE)))
    print(all.distances)

    fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "regulatoryRegion")]
    tf.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "TF")]
    targetGene.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "targetGene")]

      # add in a zero in case all of the footprints are up or downstream of the 0 coordinate, the TSS
    span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="distance"))))
    span <- max(span.endpoints) - min(span.endpoints)
    footprintLayoutFactor <- 1
    printf("initial:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    footprintLayoutFactor <- xPos.span/span

    #if(span < 600)  #
    #   footprintLayoutFactor <- 600/span
    #if(span > 1000)
    #   footprintLayoutFactor <- span/1000

    printf("corrected:  span: %d  footprintLayoutFactor: %f", span, footprintLayoutFactor)

    xPos <- as.numeric(nodeData(g, fp.nodes, attr="distance")) * footprintLayoutFactor
    yPos <- 0
    nodeData(g, fp.nodes, "xPos") <- xPos
    nodeData(g, fp.nodes, "yPos") <- yPos

    adjusted.span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="xPos"))))
    printf("raw span of footprints: %d   footprintLayoutFactor: %f  new span: %8.0f",
           span, footprintLayoutFactor, abs(max(adjusted.span.endpoints) - min(adjusted.span.endpoints)))

    tfs <- names(which(nodeData(g, attr="type") == "TF"))

    for(tf in tfs){
       footprint.neighbors <- edges(g)[[tf]]
       if(length(footprint.neighbors) > 0){
          footprint.positions <- as.integer(nodeData(g, footprint.neighbors, attr="xPos"))
          new.xPos <- mean(footprint.positions)
          if(is.na(new.xPos)) browser()
          if(is.nan(new.xPos)) browser()
          #printf("%8s: %5d", tf, new.xPos)
          }
       else{
          new.xPos <- 0
          }
       nodeData(g, tf, "yPos") <- sample(300:1200, 1)
       nodeData(g, tf, "xPos") <- new.xPos
       } # for tf

    nodeData(g, targetGene.nodes, "xPos") <- 0
    nodeData(g, targetGene.nodes, "yPos") <- -200

    g

} # addGeneModelLayout
#------------------------------------------------------------------------------------------------------------------------
test_geneRegulatoryModelToGraph <- function()
{
   printf("--- test_geneRegulatoryModelToGraph")
   all.stages <- get(load("modelsAndRegionsAllSamples.RData"))
   tbl.model <- all.stages[[1]]$model
   tbl.reg   <- all.stages[[1]]$regulatoryRegions
   tbl.reg <- subset(tbl.reg, tf %in% tbl.model$tf)
   graph <- geneRegulatoryModelToGraph("GATA2", tbl.model, tbl.reg)
   graph.lo <- addGeneModelLayout(graph)
     # rcy <- RCyjs()
   deleteGraph(rcy)
   addGraph(rcy, graph.lo)
   loadStyleFile(rcy, "trenaStyle.js")
   fit(rcy, 200)
   browser()
   xyz <- 99

} # test_geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
showStage <- function(all.stages, stageNumber, targetGene)
{
   nodes.to.hide <- c()

   for(i in seq_len(length(all.stages))){
      tfs <- all.stages[[i]]$model$gene
      regs <- all.stages[[i]]$regulatoryRegions$motif
      nodes.to.hide <- unique(c(nodes.to.hide, tfs, regs))
      }

   nodes.to.show <- c(all.stages[[stageNumber]]$model$gene,
                      all.stages[[stageNumber]]$regulatoryRegions$motif,
                      targetGene)

   hideNodes(rcy, setdiff(nodes.to.hide, nodes.to.show))
   showNodes(rcy, nodes.to.show)

} # showStage
#------------------------------------------------------------------------------------------------------------------------
run.buildMultiModelGraph <- function(sampleNumbers)
{
   printf("--- test_geneRegulatoryModelToGraph")
   models <- get(load("modelsAndRegionsAllSamples.RData")) [sampleNumbers]
   #browser()

   g.big <- buildMultiModelGraph(models, "GATA2")
   if(!exists("rcy"))
      rcy <<- RCyjs()
   g.big.lo <- addGeneModelLayout(g.big)
     # rcy <- RCyjs()
   deleteGraph(rcy)
   addGraph(rcy, g.big.lo)
   loadStyleFile(rcy, "trenaStyle.js")
   fit(rcy, 200)

   invisible(list(model=models, graph=g.big.lo))

} # run_buildMultiModelGraph
#------------------------------------------------------------------------------------------------------------------------
showSample <- function(rcy, g, models, sampleNumber)
{
   sampleNames <- c("d04_rep1", "d04_rep2", "d08_rep1", "d10_rep1", "d10_rep2", "d11_rep1",
                    "d11_rep2", "d12_rep1", "d12_rep2", "d16_rep1", "d16_rep2")

   colnames(mtx) <- c("day0.r1", "day0.r2", "day2.r1", "day2.r2",
                      "day4.r1", "day4.r2",
                      "day6.r1", "day6.r2", "day7.r1", "day7.r2",
                      "day8.r1",
                      "day8.r2", "day8.r1", "day8.r2",
                      "day10.r1", "day10.r2",
                      "day10.r1", "day10.r2",
                      "day11.r1", "day11.r2",
                      "day11.r1", "day11.r2",
                      "day12.r1", "day12.r2",
                      "day14.r1", "day14.r2",
                      "day16.r1", "day16.r2")
    mtx.coi <- c(5, 6, 11, 15, 16, 19, 20, 23, 24, 27, 28)

   #browser()
   tbl.model <- models[[sampleNumber]]$model
   tbl.reg   <- models[[sampleNumber]]$regulatoryRegions
   regs <- sprintf("%s:%d", tbl.reg$motifName, tbl.reg$distance.from.tss)
   tfs <- unique(tbl.reg$tf)
   nodes.this.sample <- c("GATA2", regs, tfs)
   nodes.to.hide <- setdiff(nodes(g), nodes.this.sample)
   hideNodes(rcy, nodes.to.hide)
   showNodes(rcy, nodes.this.sample)

   #browser()
   nodes <- c(tbl.model$tf)
   values <- tbl.model$pearsonCoeff
   setNodeAttributes(rcy, "pearson", nodes, values)

   values <- tbl.model$rfScore
   #setNodeAttributes(rcy, "rfScore", nodes, values)

   nodes <- c(nodes, "GATA2")
   values <- as.numeric(mtx.asinh.norm[nodes, mtx.coi[sampleNumber]])
   values2 <- as.numeric(mtx.raw.norm[nodes, mtx.coi[sampleNumber]])
   #browser()
   setNodeAttributes(rcy, "expression", nodes, values)

   later(function() redraw(rcy), 0.25)
   xyz <- 99

} # showSample
#------------------------------------------------------------------------------------------------------------------------
talk <- function()
{
   x <- run.buildMultiModelGraph(1:11) # c(2,3))
   loadStyleFile(rcy, "trenaStyle.js")
   models <- x$model
   g <- x$graph

   showSample(rcy, g, models, 1)
   showSample(rcy, g, models, 2)
   showSample(rcy, g, models, 3)
   showSample(rcy, g, models, 4)
   showSample(rcy, g, models, 5)
   showSample(rcy, g, models, 6)
   showSample(rcy, g, models, 7)
   showSample(rcy, g, models, 8)
   showSample(rcy, g, models, 9)
   showSample(rcy, g, models, 10)
   showSample(rcy, g, models, 11)

} # talk
#------------------------------------------------------------------------------------------------------------------------
