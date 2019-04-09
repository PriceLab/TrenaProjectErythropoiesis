library(RCyjs)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("parseChromLocString"))
   source("~/github/trena/R/utils.R")
if(!exists("tbl.geneInfo"))
   tbl.geneInfo <- get(load((system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))))
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
   nodeDataDefaults(g, attr = "rfScore") <- 0
   nodeDataDefaults(g, attr = "pcaMax") <- 0
   nodeDataDefaults(g, attr = "concordance") <- 0
   nodeDataDefaults(g, attr = "betaLasso") <- 0
   nodeDataDefaults(g, attr = "motif") <- ""
   nodeDataDefaults(g, attr = "xPos") <- 0
   nodeDataDefaults(g, attr = "yPos") <- 0

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

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, regRegions.names, "type")  <- "regulatoryRegion"
   nodeData(g, all.nodes, "label")  <- all.nodes
   nodeData(g, regRegions.names, "label") <- tbl.reg$motifName
   nodeData(g, regRegions.names, "distance") <- tbl.reg$distance.from.tss
   nodeData(g, regRegions.names, "motif") <- tbl.reg$motifName

   nodeData(g, tfs, "pearson") <- tbl.gm$pearsonCoeff
   nodeData(g, tfs, "betaLasso") <- tbl.gm$betaLasso
   nodeData(g, tfs, "rfScore") <- tbl.gm$rfScore

   g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
   edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"

   g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
   edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"

   g

} # geneRegulatoryModelToGraph
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
fixModel <- function(tbl.model)
{
   colnames(tbl.model)[grep("^gene$", colnames(tbl.model))] <- "tf"
   colnames(tbl.model)[grep("^rfScore$", colnames(tbl.model))] <- "rfScore"
   return(tbl.model)

} # fixModel
#------------------------------------------------------------------------------------------------------------------------
fixReg <- function(tbl.reg, targetGene)
{
   colnames(tbl.reg)[grep("^start$", colnames(tbl.reg))] <- "motifStart"
   colnames(tbl.reg)[grep("^stop$", colnames(tbl.reg))] <- "motifEnd"
   chromLocStrings <- tbl.reg$sequence_name
   locs <- lapply(chromLocStrings, parseChromLocString)
   tss <- subset(tbl.geneInfo, geneSymbol==targetGene)$tss
   strand <- subset(tbl.geneInfo, geneSymbol==targetGene)$strand

   tbl.locs <- as.data.frame(do.call(rbind, locs))
   tbl.new <- cbind(tbl.locs, tbl.reg)
   tbl.new$chrom <- as.character(tbl.new$chrom)
   tbl.new$start <- as.numeric(tbl.new$start)
   tbl.new$end <- as.numeric(tbl.new$end)
   tbl.new$motifStart <- tbl.new$motifStart + tbl.new$start
   tbl.new$motifEnd   <- tbl.new$motifEnd   + tbl.new$start
   tbl.new$distance.from.tss <- tbl.new$motifStart - tss
   tbl.new$pValue <- -log10(tbl.new$pValue)

   tbl.new <- tbl.new[, c("motif", "chrom", "motifStart", "motifEnd", "strand",  "pValue", "score",
                          "matched_sequence", "distance.from.tss", "tf")]
   colnames(tbl.new) <- required.regulatoryRegionsColumnNames

   tbl.new

} # fixReg
#------------------------------------------------------------------------------------------------------------------------
test_fixReg <- function()
{
   printf("--- test_fixReg")
   targetGene <- "GATA2"
   all.stages <- get(load("modelsAndRegionsAllSamples.RData"))

   tbl.reg   <- all.stages[[1]]$regulatoryRegions
   tbl.fixed <- fixReg(tbl.reg, targetGene)
   checkEquals(colnames(tbl.fixed), required.regulatoryRegionsColumnNames)


} # test_fixReg
#------------------------------------------------------------------------------------------------------------------------
test_geneRegulatoryModelToGraph <- function()
{
   printf("--- test_geneRegulatoryModelToGraph")
   all.stages <- get(load("modelsAndRegionsAllSamples.RData"))
   tbl.model <- all.stages[[1]]$model
   tbl.reg   <- all.stages[[1]]$regulatoryRegions
   tbl.model <- fixModel(tbl.model)
   tbl.reg <- fixReg(tbl.reg, "GATA2")
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
