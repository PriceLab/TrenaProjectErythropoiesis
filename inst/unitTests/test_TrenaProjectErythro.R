library(TrenaProjectErythropoiesis)
library(RUnit)
library(trenaSGM)
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
   test_expressionMatrices()
   test_genomicRegions()
   test_setTargetGene()
   # test_singleGeneModel_NFE2()

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
                 "ATAC_brand_d12_rep2", "ATAC_brand_d16_rep1", "ATAC_brand_d16_rep2")

   dataset.names <- getGenomicRegionsDatasetNames(tpe)
   checkTrue(all(expected %in% dataset.names))
   checkTrue(!"prepFiles.R" %in% dataset.names)

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
   column.subset.expected <- c("chrom", "start", "end", "gene", "eqtl", "hic")
   checkTrue(all(column.subset.expected %in% colnames(tbl.enhancers)))
   checkTrue(nrow(tbl.enhancers) >= 20)

   message(sprintf("    geneGeneEnhancersRegion"))
   region <- getGeneEnhancersRegion(tpe, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "chr11:5227061-5384601")

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tpe)
   checkTrue(nrow(tbl.dhs) > 100)

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- getChipSeq(tpe, chrom="chr11", start=5225464, end=5229395, tfs=NA)
   checkTrue(nrow(tbl.chipSeq) > 130)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
demo_targetsOfGATA2Repression <- function()
{
   file <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression", "GTEX.lymphocyte.rna-seq-geneSymbols.21415x130.RData")
   mtx.lymph <- get(load(file))
   file <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression","GTEX.wholeBlood.rna-seq-geneSymbols.22330x407.RData")
   mtx.blood <- get(load(file))
   file <- system.file(package="TrenaProjectErythropoiesis", "extdata", "expression", "brandLabDifferentiationTimeCourse-27171x28.RData")
   mtx.marjorie <- get(load(file))

   cors.lymph <- lapply(rownames(mtx.lymph), function(gene) cor(mtx.lymph["GATA2",], mtx.lymph[gene,]))
   names(cors.lymph) <- rownames(mtx.lymph)
   cors.blood <- lapply(rownames(mtx.blood), function(gene) cor(mtx.blood["GATA2",], mtx.blood[gene,]))
   names(cors.blood) <- rownames(mtx.blood)
   cors.marjorie <- lapply(rownames(mtx.marjorie), function(gene) cor(mtx.marjorie["GATA2",], mtx.marjorie[gene,]))
   names(cors.marjorie) <- rownames(mtx.marjorie)

   all.genes <- sort(unique(c(names(cors.lymph), names(cors.blood), names(cors.marjorie))))
   count <- length(all.genes)
   tbl <- data.frame(marj=rep(0, count), lymph=rep(0, count), blood=rep(0, count), stringsAsFactors=FALSE)
   rownames(tbl) <- all.genes
   tbl[names(cors.marjorie), "marj"] <- as.numeric(cors.marjorie)
   tbl[names(cors.lymph), "lymph"] <- as.numeric(cors.lymph)
   tbl[names(cors.blood), "blood"] <- as.numeric(cors.blood)

   subset(tbl, marj < -0.5 & blood < -0.35)
     #              marj        lymph      blood
     # ACSM4   -0.7167251 -0.039075656 -0.3653664
     # ADAMTS2 -0.5013062  0.288090374 -0.5290857
     # ALPL    -0.6282120  0.328399161 -0.3974193
     # AP3B2   -0.5258294  0.153368658 -0.3829907
     # ARG1    -0.7845579  0.193346004 -0.3680406
     # COX6B2  -0.5461555 -0.009669172 -0.3524632
     # NFE2    -0.5019009  0.301234970 -0.3636447
     # NMNAT2  -0.8958629  0.201230740 -0.4288662
     # SLC22A4 -0.7884134  0.318294382 -0.3947438

   subset(tbl, marj > 0.85 & blood > 0.7)


} # demo_targetsOfGATA2Repression
#------------------------------------------------------------------------------------------------------------------------
demo_NFE2_models <- function()
{
   library(TrenaProjectLymphocyte)
   library(org.Hs.eg.db)
   tp <- TrenaProjectLymphocyte();

   genome <- "hg38"
   targetGene <- "NFE2"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)

   chromosome <- tbl.info$chrom
   tss <- tbl.info$tss
      # strand-aware start and end: atf1 is on the + strand
   start <- tss - 50000
   end   <- tss + 50000

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)
   file <- system.file(package="TrenaProjectLymphocyte", "extdata", "expression","GTEX.wholeBlood.rna-seq-geneSymbols.22330x407.RData")
   mtx.blood <- get(load(file))
   file <- system.file(package="TrenaProjectErythropoiesis", "extdata", "expression", "brandLabDifferentiationTimeCourse-27171x28.RData")
   mtx.marjorie <- get(load(file))


      #----------------------------------------------------------------------------------------------------
      # first, build a model with "placenta2", an early version of the placenta footprint database
      #----------------------------------------------------------------------------------------------------

   recipe <- list(title="NFE2",
                  type="footprint.database",
                  regions=tbl.regions,
                  geneSymbol=targetGene,
                  tss=tss,
                  matrix=mtx.blood,
                  db.host="khaleesi.systemsbiology.net",
                  db.port=5432,
                  databases=list("lymphoblast_hint_16", "lymphoblast_hint_20"),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="builtinFimo",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.1,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "xgboost"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  recipe, quiet=FALSE)
   x <- build(fpBuilder)

   recipe$matrix <- mtx.marjorie
   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  recipe, quiet=FALSE)
   x2 <- build(fpBuilder)


   file <- "~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74246_RNAseq_All_Counts.txt"
   file.exists(file)
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE,nrow=-1)
   rownames(tbl) <- tbl[, 1]
   tbl <- tbl[, -1]
   mtx <- asinh(as.matrix(tbl))
   fivenum(mtx)

   suppressWarnings(
      logTimingInfo("placenta2 db, +/- 5bk generic promoter on ATF1", system.time(x2 <- build(fpBuilder)))
      )


} # demo_NFE2_models
#------------------------------------------------------------------------------------------------------------------------
demo_extractEarlyTimePointHighlyExpressedGenes <- function()
{
   message(sprintf("--- demo_extractEarlyTimePointHighlyExpressedGenes"))
   library(jsonlite)
   library(httr)

   expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")
   checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
   mtx <- getExpressionMatrix(tpe, expected[1])

   cutoff <- 8
   goi <- names(which(sapply(rownames(mtx), function(r) all(mtx[r, 1:4] > cutoff))))
   goi.string <- toJSON(goi)
   uri <- sprintf("http://localhost:8000/goEnrich")
   body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=goi)))

   r <- POST(uri, body=body.jsonString)

      #sprintf('{"geneSymbols": "%s"}', goi.string))
   tbl <- fromJSON(content(r)[[1]])
   dim(tbl)
   tbl[c(1,6,8,21),-1][,c(6,5,1, 7)]
   tbl[grep("erythro", tbl$Term, ignore.case=TRUE),]

   uri <- sprintf("http://localhost:8000/keggEnrich")
   body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=goi)))

   r <- POST(uri, body=body.jsonString)

      #sprintf('{"geneSymbols": "%s"}', goi.string))
   tbl.kegg <- fromJSON(content(r)[[1]])



   cutoff.1 <- 10
   cutoff.2 <- 7

   goi <- names(which(sapply(rownames(mtx), function(geneName) all(mtx[geneName, 1:2] > cutoff.1))))
   print(length(goi))
   goi <- names(which(sapply(goi, function(geneName) all(mtx[geneName, 3:4] < cutoff.2))))
   print(length(goi))

   goi.string <- toJSON(goi)
   uri <- sprintf("http://localhost:8000/goEnrich")
   body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=goi)))

   r <- POST(uri, body=body.jsonString)

      #sprintf('{"geneSymbols": "%s"}', goi.string))
   tbl <- fromJSON(content(r)[[1]])
   dim(tbl)
   tbl[c(1,6,8,21),-1][,c(6,5,1, 7)]





} # demo_extractEarlyTimePointHighlyExpressedGenes
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
