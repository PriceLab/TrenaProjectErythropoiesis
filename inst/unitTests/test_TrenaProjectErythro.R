library(TrenaProjectErythropoiesis)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tpe")) {
   message(sprintf("--- creating instance of TrenaProjectErythropoiesis"))
   tpe <- TrenaProjectErythropoiesis();
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectErythropoiesis", "TrenaProject") %in% is(tpe)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("Abca1")
   checkTrue(all(subset.expected %in% getSupportedGenes(tpe)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(getVariantDatasetNames(tpe), character(0))

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("extraembryonic_structure_wellington_16", "extraembryonic_structure_wellington_20",
                 "extraembryonic_structure_hint_16", "extraembryonic_structure_hint_20")
   checkTrue(is.na(getFootprintDatabaseNames(tpe)))
   checkTrue(is.na(getFootprintDatabaseHost(tpe)))

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c()
   checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))

   # mtx <- getExpressionMatrix(tpe, expected[1])
   # checkEquals(dim(mtx), c(6890, 14))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tpe, "HBB")
   checkEquals(getTargetGene(tpe), "HBB")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tpe)
   checkTrue(nrow(tbl.transcripts) >= 1)
   checkEquals(tbl.transcripts$chr, "chr11")
   checkEquals(tbl.transcripts$start, 5225464)
   checkEquals(tbl.transcripts$end,   5229395)
   checkEquals(tbl.transcripts$tss, 5227071)
   checkEquals(tbl.transcripts$strand, -1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tpe, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr11:5225464-5229395")

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tpe)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 20)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tpe, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr11:5227061-5379996")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tpe)
   checkTrue(nrow(tbl.dhs) > 100)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- getChipSeq(tpe, chrom="chr11", start=5225464, end=5229395, tfs=NA)
   checkTrue(nrow(tbl.chipSeq) > 130)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel_noDNA <- function()
{
   printf("--- test_buildSingleGeneModel_noDNA")
   require("trenaSGM")

   genome <- "hg38"
   targetGene <- "GATA2"

   mtx.name <- "brandLabDifferentiationTimeCourseRNA"
   checkTrue(mtx.name %in% getExpressionMatrixNames(tpe))
   mtx <- asinh(getExpressionMatrix(tpe, mtx.name))


   candidate.tfs <- allKnownTFs()
   candidate.tfs <- c("MAZ", "SP4", "SP2", "SP3", "SP1", "ZFX", "ZN148")

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
   tbl.model <- x$model
   dim(tbl.model)

   checkEquals(dim(x$regulatoryRegions), c(0,0))
   checkTrue(all(tbl.model$peasonCoeff > 0.7))
      # the order
   checkEquals(tbl.model$gene, c("PLEK", "IRF5", "IKZF1", "LYL1", "SPI1", "TFEC"))

} # test_buildSingleGeneModel_noDNA
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
