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
   test_genomicRegions()
   test_setTargetGene()

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

   subset.expected <- c("GATA2")
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
   expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")

   checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
   mtx <- getExpressionMatrix(tpe, expected[1])
   checkEquals(dim(mtx), c(16173, 28))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
test_genomicRegions <- function()
{
   expected <- c("ATAC_brand_d04_rep1", "ATAC_brand_d04_rep2", "ATAC_brand_d08_rep1", "ATAC_brand_d10_rep1",
                 "ATAC_brand_d10_rep2", "ATAC_brand_d11_rep1", "ATAC_brand_d11_rep2", "ATAC_brand_d12_rep1",
                 "ATAC_brand_d12_rep2", "ATAC_brand_d16_rep1", "ATAC_brand_d16_rep2", "prepFiles.R")

   checkTrue(all(expected %in% getGenomicRegionsDatasetNames(tpe)))

   tbl.atac <- getGenomicRegionsDataset(tpe, "ATAC_brand_d04_rep2")
   checkEquals(dim(tbl.atac), c(133208, 10))

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
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
