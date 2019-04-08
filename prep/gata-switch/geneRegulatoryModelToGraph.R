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
   browser()

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
   nodeDataDefaults(g, attr = "randomForest") <- 0
   nodeDataDefaults(g, attr = "pcaMax") <- 0
   nodeDataDefaults(g, attr = "concordance") <- 0
   nodeDataDefaults(g, attr = "betaLasso") <- 0
   nodeDataDefaults(g, attr = "motif") <- ""
   nodeDataDefaults(g, attr = "xPos") <- 0
   nodeDataDefaults(g, attr = "yPos") <- 0

   edgeDataDefaults(g, attr = "edgeType") <- "undefined"

   tfs <- tbl.gm$tf

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
   nodeData(g, regRegions.names, "distance") <- tbl.reg$distance
   nodeData(g, regRegions.names, "motif") <- tbl.reg$motifName

   nodeData(g, tfs, "pearson") <- tbl.gm$pearson
   nodeData(g, tfs, "betaLasso") <- tbl.gm$betaLasso
   nodeData(g, tfs, "randomForest") <- tbl.gm$randomForest
   nodeData(g, tfs, "pcaMax") <- tbl.gm$pcaMax
   nodeData(g, tfs, "concordance") <- tbl.gm$concordance

   #browser()
   g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
   edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"

   g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
   edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"

   g

} # geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
test_geneRegulatoryModelToGraph <- function()
{
   printf("--- test_geneRegulatoryModelToGraph")
   all.stages <- get(load("modelsAndRegionsAllSamples.RData"))
   tbl.model <- all.stages[[1]]$model
   tbl.reg   <- all.stages[[1]]$regulatoryRegions
   tbl.model <- fixModel(tbl.model)
   tbl.reg <- fixReg(tbl.reg)
   graph <- geneRegulatoryModelToGraph("GATA2", tbl.model, tbl.reg)

} # test_geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
fixModel <- function(tbl.model)
{
   colnames(tbl.model)[grep("^gene$", colnames(tbl.model))] <- "tf"
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
