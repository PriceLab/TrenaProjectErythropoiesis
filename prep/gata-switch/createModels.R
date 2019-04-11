library(TrenaProjectErythropoiesis)
library(igvR)
library(trenaSGM)
library(FimoClient)
library(RUnit)
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
mtx <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
mtx <- asinh(mtx)

# if(!exists("tpe")){
#    tpe <- TrenaProjectErythropoiesis()
#    mtx.name <- "brandLabDifferentiationTimeCourseRNA"
#    stopifnot(mtx.name %in% getExpressionMatrixNames(tpe))
#    mtx <- asinh(getExpressionMatrix(tpe, mtx.name))
#    }


if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   }
#------------------------------------------------------------------------------------------------------------------------
tbl.lcr <- data.frame(chrom="chr11", start=5269925, end=5304186, stringsAsFactors=FALSE)
tbl.hs2 <- data.frame(chrom="chr11", start=5280519, end=5281050, stringsAsFactors=FALSE)
#------------------------------------------------------------------------------------------------------------------------
displayGeneHancer <- function(gene)
{
   setTargetGene(tpe, gene)
   tbl.enhancers <- getEnhancers(tpe)
   track <- DataFrameQuantitativeTrack("GH", tbl.enhancers[, c("chrom", "start", "end", "combinedScore")],
                                       "brown", autoscale=FALSE, min=0, max=100)
   displayTrack(igv, track)
   with(tbl.enhancers, showGenomicRegion(igv, sprintf("%s:%d-%d", unique(chrom), min(start)-1000, max(end) + 1000)))

} # displayGeneHancer
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function(chromosome, start.loc, end.loc)
{
   directory <- "../import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
test_getATACseq <- function()
{
   printf("--- test_getATACseq")
   tbl.atac <- getATACseq("chr3", 128495142, 128498398)
   samples <- unique(tbl.atac$sample)
    # "d08_rep1" "d10_rep1" "d10_rep2" "d11_rep1" "d11_rep2" "d12_rep1" "d12_rep2" "d16_rep1" "d16_rep2"
   for(sample.x in samples){
      tbl.sample <- subset(tbl.atac, sample==sample.x)[, c("chrom", "start", "end")]
      write.table(tbl.sample, file=sprintf("tbl.%s.bed", sample.x), quote=FALSE, row.names=FALSE)
      } # for sample.x

} # test_getATACseq
#------------------------------------------------------------------------------------------------------------------------
displayATACseq <- function(chromosome, start.loc, end.loc)
{
   directory <- "../import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      dim(tbl.atac.region)
      #browser()
      #xyz <- "displayATACseq"
      track <- DataFrameQuantitativeTrack(track.name, tbl.atac.region[, c("chrom", "start", "end", "c10")],
                                          "blue", autoscale=FALSE, min=0, max=430)
      displayTrack(igv, track)
      } # files

} # displayATACseq
#------------------------------------------------------------------------------------------------------------------------
findTFs <- function(tbl.regions, threshold=1e-4)
{
    FIMO_HOST <- "khaleesi"
    FIMO_PORT <- 60000

   if(!exists("fc")){
      fc <<- FimoClient(FIMO_HOST, FIMO_PORT)
      }

    tbl.matches <- requestMatchForRegions(fc, tbl.regions, "hg38", threshold)
    motif.names <- tbl.matches$motif
    deleters <- setdiff(motif.names, names(MotifDb))
    for(deleter in deleters){
       #browser()
       deleter.indices <- grep(deleter, tbl.matches$motif)
       if(length(deleter.indices) > 0)
          tbl.matches <- tbl.matches[-deleter.indices,]
       }
    tfs <- mcols(MotifDb[tbl.matches$motif])$geneSymbol
    tbl.matches$tf <- tfs

    tbl.matches

} # findTFs
#------------------------------------------------------------------------------------------------------------------------
test_findTFs <- function()
{
   printf("--- test_findTFs")
   tbl.regions <- data.frame(chrom="chr3", start=128483072, end=128483461, stringsAsFactors=FALSE)
   tbl.matches <- findTFs(tbl.regions, threshold=1e-6)

} # test_findTFs
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(candidate.tfs, targetGene)
{
   genome <- "hg38"


   build.spec <- list(title="gata2.noDNA.allTFs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx,
                      tfPool=allKnownTFs(),
                      tfs=candidate.tfs,
                      tfPrefilterCorrelation=0.1,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   x

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
reformatMotifHits <- function(tbl.hits)
{
   tbl.out <- do.call(rbind, lapply(lapply(tbl.hits$sequence_name, parseChromLocString), as.data.frame))
   tbl.out$chrom <- as.character(tbl.out$chrom)
   tfGeneSymbols <- mcols(MotifDb[tbl.hits$motif])$geneSymbol
   tbl.out$tf <- tfGeneSymbols
   tbl.out$motif <- tbl.hits$motif

   tbl.out

} # reformatMotifHits
#------------------------------------------------------------------------------------------------------------------------
test_reformatMotifHits <- function(tbl.hits)
{
   printf("--- test_reformatMotifHits")
   load("tbl.motifHits.model.1.RData")
   tbl.fixed <- reformatMotifHits(tbl.motifHits)
   checkEquals(nrow(tbl.fixed), nrow(tbl.motifHits))
   checkEquals(colnames(tbl.fixed), c ("chrom", "start", "end", "tf", "motif"))

} # test_reformatMotifHits
#------------------------------------------------------------------------------------------------------------------------
test_buildModel <- function()
{
   printf("--- test_buildModel")

   tbl.regions.1 <- data.frame(chrom="chr3", start=128483072, end=128483461, stringsAsFactors=FALSE)
   tbl.matches.1 <- findTFs(tbl.regions.1, threshold=1e-4)
   x.1 <- buildModel(unique(tbl.matches.1$tf), "GATA2")
   x.1$regulatoryRegions <- reformatMotifHits(tbl.matches.1)
   tbl.model.1 <- x.1$model
   new.order.1 <- order(tbl.model.1$rfScore,decreasing=TRUE)
   tbl.model.1 <- tbl.model.1[new.order.1,]

   tbl.model <- tbl.model.1
   tbl.regions <- x.1$regulatoryRegions
   save(mtx, tbl.model, tbl.regions, file="gata2.model.RData")

   tbl.regions.2 <- data.frame(chrom="chr3", start=128487534, end=128488231, stringsAsFactors=FALSE)
   tbl.matches.2 <- findTFs(tbl.regions.2, threshold=1e-4)
   x.2 <- buildModel(unique(tbl.matches.2$tf), "GATA2")
   tbl.model.2 <- x.2$model
   new.order.2 <- order(tbl.model.2$rfScore,decreasing=TRUE)
   tbl.model.2 <- tbl.model2[new.order.2,]

   tbl.regions.3 <- data.frame(chrom="chr3", start=128481834, end=128493658, stringsAsFactors=FALSE)
   tbl.matches.3 <- findTFs(tbl.regions.3, threshold=1e-4)
   x.3 <- buildModel(unique(tbl.matches.3$tf), "GATA2")
   tbl.model.3 <- x.3$model
   new.order.3 <- order(tbl.model.3$rfScore,decreasing=TRUE)
   tbl.model.3 <- tbl.model3[new.order.3,]

   #gata
   #save(list(x1, x2, x3

} # test_buildModel
#------------------------------------------------------------------------------------------------------------------------
trimModel <- function(tbl.model, tbl.reg, tf.keepers=c(), votesNeeded=3)
{
   matched.keeper.rows <- unlist(lapply(tf.keepers, function(tf) grep(tf, tbl.model$gene)))

   # pvals.scores <- 1og10(tbl.model$lassoPValue) * -1
   pvals <- -log10(tbl.model$lassoPValue)
   good.lassoPval <- which(pvals > 5)

   good.betaLasso <- which(abs(tbl.model$betaLasso) > 0.1)
   good.betaRidge <- which(abs(tbl.model$betaRidge) > 0.1)

   spearman.cutoff <- fivenum(abs(tbl.model$spearmanCoeff))[4]
   good.spearmanCoeff <- which(abs(tbl.model$spearmanCoeff) >= spearman.cutoff)

   randomForest.cutoff <- fivenum(tbl.model$rfScore)[4]
   forest.tfs <- subset(tbl.model, rfScore >= randomForest.cutoff)$gene
   good.rfScore <- unlist(lapply(forest.tfs, function(tf) grep(tf, tbl.model$gene)))

   all.counts <- c(good.lassoPval, good.betaLasso, good.betaRidge, good.spearmanCoeff, good.rfScore)
   tbl.freq <- as.data.frame(table(all.counts), stringsAsFactors=FALSE)
   colnames(tbl.freq) <- c("rowNumber", "count")
   tbl.freq <- tbl.freq[order(tbl.freq$count, decreasing=TRUE),]
   tbl.freq$rowNumber <- as.integer(tbl.freq$rowNumber)
   good.tf.rows <- subset(tbl.freq, count >= votesNeeded)$rowNumber
   not.yet.included <- setdiff(matched.keeper.rows, good.tf.rows)
   if(length(not.yet.included) > 0)
      good.tf.rows <- c(good.tf.rows, not.yet.included)

   tbl.model <- tbl.model[good.tf.rows,]
   new.order <- order(abs(tbl.model$spearmanCoeff), decreasing=TRUE)
   tbl.model <- tbl.model[new.order,]
   browser()

   tbl.reg <- subset(tbl.reg, tf %in% tbl.model$gene)
   return(list(model=tbl.model, regulatoryRegions=tbl.reg))

   xyz <- 99

} # trimModel
#------------------------------------------------------------------------------------------------------------------------
test_trimModel <- function()
{
   printf("--- test_trimModel")
   x.raw <- get(load("modelsAndRegionsAllSamples-rawNeedTrimming.RData"))

   model.number <- 9
   tbl.model <- x.raw[[model.number]]$model
   tbl.reg <- x.raw[[model.number]]$regulatoryRegions
   x.trimmed <- trimModel(tbl.model, tbl.reg, c("GATA1", "GATA2"), votesNeeded=2)
   x.final <- fixModelAndRegions(x.trimmed, targetGene="GATA2")
   checkTrue(all(x.final$regulatoryRegions$tf %in% x.final$model$tf))
   browser()
   xyz <- "after trimming"

} # test_trimModel
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
fixModelAndRegions <- function(x, targetGene)
{
   result <- list(model=fixModel(x$model), regulatoryRegions=fixReg(x$regulatoryRegions, targetGene))
   result

} # fixModelAndRegions
#------------------------------------------------------------------------------------------------------------------------
buildModelsAtEachStage <- function()
{
   tbl.atac.all <- getATACseq("chr3", 128481798, 128498370)
   samples <- unique(tbl.atac.all$sample)
   models <- list()

   for(currentSample in samples){
      tbl.regions <- subset(tbl.atac.all, sample==currentSample)[, c("chrom", "start", "end")]
      tbl.matches <- findTFs(tbl.regions, threshold=1e-5)
      x <- buildModel(unique(tbl.matches$tf), "GATA2")
      x$regulatoryRegions <- tbl.matches
      browser()
      x2 <- trimModel(x$model, x$regulatoryRegions, c("GATA1", "GATA2"), votesNeeded=3)
      x3 <- fixModelAndRegions(x2, targetGene="GATA2")
      browser()
      x3$model$sample <- currentSample
         # just keep regulatory regions whose cognate TF are in the model
      #tbl.matches <- subset(tbl.matches, tf %in% x3$model$gene)
      #x3$regulatoryRegions <- tbl.matches
      models[[currentSample]] <- x3
      }

   filename <- "modelsAndRegionsAllSamples.RData"
   printf("saving %d models to %s", length(models), filename)
   save(models, file=filename)
   browser()
   xyz <- 99

} # buildModelsAtEachStage
#------------------------------------------------------------------------------------------------------------------------
