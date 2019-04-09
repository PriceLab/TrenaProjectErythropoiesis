library(TrenaProjectErythropoiesis)
library(igvR)
library(trenaSGM)
library(FimoClient)
library(RUnit)
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
   x <- getATACseq("chr3", 128495142, 128498398)

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
buildModelsAtEachStage <- function()
{
   tbl.atac.all <- getATACseq("chr3", 128481798, 128498370)
   samples <- unique(tbl.atac.all$sample)
   models <- list()

   for(currentSample in samples){
      tbl.regions <- subset(tbl.atac.all, sample==currentSample)[, c("chrom", "start", "end")]
      tbl.matches <- findTFs(tbl.regions, threshold=1e-5)
      x <- buildModel(unique(tbl.matches$tf), "GATA2")
      x$model$sample <- x$model
         # just keep regulatory regions whose cognate TF are in the model
      tbl.matches <- subset(tbl.matches, tf %in% x$model$gene)
      x$regulatoryRegions <- tbl.matches
      models[[currentSample]] <- x
      }

   filename <- "modelsAndRegionsAllSamples.RData"
   printf("saving %d models to %s", length(models), filename)
   save(models, file=filename)

} # buildModelsAtEachStage
#------------------------------------------------------------------------------------------------------------------------
#
#    candidate.tfs <- rownames(subset(tbl.cor, abs(cor) > 0.4))
#       #  [1] "BACH1"   "BATF"    "DMC1"    "ELF3"    "ELF5"    "ETV5"    "GTF2F1"  "IKZF2"   "IRF1"
#       # [10] "KLF1"    "KLF13"   "KLF16"   "MAFG"    "MAFK"    "MESP1"   "MTA3"    "NFE2"    "NFIC"
#       # [19] "PATZ1"   "PAX6"    "PLAGL1"  "PML"     "RARB"    "RARG"    "RXRA"    "SCRT2"   "SIN3A"
#       # [28] "SMARCC1" "SP2"     "STAT2"   "TAL1"    "TCF12"   "TCF4"    "VDR"     "WRNIP1"  "YY1"
#       # [37] "ZNF219"  "ZNF263"  "ZNF281"  "ZNF384"  "ZNF740"
#
#    build.spec <- list(title="hbb.hs2.motifs",
#                       type="noDNA.tfsSupplied",
#                       matrix=mtx.asinh,
#                       tfPool=allKnownTFs(),
#                       tfs=candidate.tfs,
#                       tfPrefilterCorrelation=0.2,
#                       annotationDbFile=dbfile(org.Hs.eg.db),
#                       orderModelByColumn="rfScore",
#                       solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
#                       quiet=TRUE)
#
#    builder <- NoDnaModelBuilder(genomeName, targetGene,  build.spec, quiet=TRUE)
#    x <- build(builder)
#
# } # buildModel
# #------------------------------------------------------------------------------------------------------------------------
# mdb.1 <- query(MotifDb, c("hsapiens", "jaspar2018", "gata2"))
# mdb.2 <- query(MotifDb, c("hsapiens", "gata1", "hocomoco"))
# mdb.3 <- query(MotifDb, c("hsapiens", "gata1", "tal1", "jaspar2018"))
# export(c(mdb.1, mdb.2, mdb.3), con="~/github/fimoService/pfms/gata1-tal1-gata2.meme", format="meme")
# explore.gata1.gata2.tal1 <- function()
# {
#    library(igvR)
#    igv <- igvR()
#    setGenome(igv, "hg38")
#    showGenomicRegion(igv, with(tbl.hs2, sprintf("%s:%d-%d", chrom, start, end)))
#    tbl.hs2.motifs <- requestMatchForRegions(fc, tbl.hs2, "hg38", 0.0013)
#    tbl.bedGraph <- tbl.hs2.motifs
#    tbl.bedGraph$start <- tbl.bedGraph$start + tbl.hs2$start
#    tbl.bedGraph$stop <- tbl.bedGraph$stop + tbl.hs2$start
#    tbl.bedGraph$chrom <- tbl.hs2$chrom
#    tbl.bedGraph <- tbl.bedGraph[, c("chrom", "start", "stop", "motif", "score")]
#
#
#    track <- DataFrameAnnotationTrack("GATA1", tbl.bedGraph[1,], color="brown", displayMode="EXPANDED", trackHeight=30)
#    displayTrack(igv, track)
#    track <- DataFrameAnnotationTrack("GATA2", tbl.bedGraph[4,], color="darkgreen", displayMode="EXPANDED", trackHeight=30)
#    displayTrack(igv, track)
#    track <- DataFrameAnnotationTrack("GATA1::TAL1", tbl.bedGraph[7,], color="blue", displayMode="EXPANDED", trackHeight=30)
#    displayTrack(igv, track)
#
#
#
# } # explore.gata1.gata2.tal1
#------------------------------------------------------------------------------------------------------------------------

