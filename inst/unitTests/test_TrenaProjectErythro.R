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
   test_expressionMatrices()
   test_setTargetGene()
   test_buildFimoDatabaseSingleGeneModel_TBX15()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectErythropoiesis", "TrenaProjectHG38") %in% is(tpe)))

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
test_bindingSiteDatabases <- function()
{
   message(sprintf("--- test_bindingSiteDatabases"))

   expected <- c("extraembryonic_structure_wellington_16", "extraembryonic_structure_wellington_20",
                 "extraembryonic_structure_hint_16", "extraembryonic_structure_hint_20")
   checkTrue(is.na(getFootprintDatabaseNames(tpe)))
   checkTrue(is.na(getFootprintDatabaseHost(tpe)))

} # test_bindingSiteDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("brandLabDifferentiationTimeCourse-filtered")
   checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
   mtx <- getExpressionMatrix(tpe, expected[1])
   checkEquals(dim(mtx), c(16173, 28))

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
test_buildFimoDatabaseSingleGeneModel_TBX15 <- function()
{
   printf("--- test_buildFimoDatabaseSingleGeneModel_TBX15")
   require("trenaSGM")

   genome <- "hg38"
   targetGene <- "GATA2"
   setTargetGene(tpe, targetGene)

   chromosome <- "chr3"
   tss <- 128493185
   upstream <- 5000
   downstream <- 5000
      # strand-aware start and end: GATA2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   mtx.name <- "brandLabDifferentiationTimeCourse-filtered"
   checkTrue(mtx.name %in% getExpressionMatrixNames(tpe))
   mtx <- getExpressionMatrix(tpe, mtx.name)

   db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs", "gata2.gh.fimoBindingSites.sqlite")
   checkTrue(file.exists(db.name))

   build.spec <- list(title="fimo.5000up.5000down",
                      type="fimo.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host="localhost",
                      db.port=NA_integer_,
                      databases=list(db.name),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="fimoDatabase",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      maxModelSize=10,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- FimoDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)

   checkTrue(all(c("model", "regulatoryRegions") %in% names(x)))
   tbl.model <- x$model
   checkEquals(nrow(tbl.model), 10)

} # test_buildFimoDatabaseSingleGeneModel_TBX15
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
