if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 1)
   targetGene <- args[1]
   }

library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")

tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
tbl.fimo <- get(load(sprintf("tbl.fimo.%s.RData", targetGene)))
dim(tbl.fimo)                 # 114021 9
trenaProject <- TrenaProjectErythropoiesis()
mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-rna-pkfm-27170x14.RData"))
mtx.srm <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-srm-copyNumber-100x13.RData"))

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_promoterOnlyFimoTable()
    test_inclusiveFimoTable()

} # runTests
#----------------------------------------------------------------------------------------------------
test_promoterOnlyFimoTable <- function()
{

   message(sprintf("--- test_promoterOnlyFimoTable"))

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo.promoter.only, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo.promoter.only))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.7 & gh > 1)$tf
   length(tfs)

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna)

   checkEquals(unique(tbl.trena.tf$class), "tf")
   checkTrue(all(c("RXRA", "ZBTB7A", "SP2", "NFE2", "KLF1", "KLF16") %in% tbl.trena.tf$gene[1:10]))

   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)

   checkTrue(all(c("tf", "rbp") %in% tbl.trena.both$class))
   checkTrue(all(c("STAU1", "RBM15", "MBNL2", "RXRA", "ALKBH5", "ZBTB7A") %in% tbl.trena.both$gene[1:10]))
   checkEquals("STAU1", tbl.trena.both$gene[1])

} # test_inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
build.models.with.inclusiveFimoTable <- function()
{
   message(sprintf("--- test_inclusiveFimoTable: %s", targetGene))

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-3 &
                          abs(cor.all) > 0.5 &
                          (gh > 10 | oc))$tf
   printf("running trena with candidate tfs: %d", length(unique(tfs)))

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna)
   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)
   tbl.trena.tf$target <- targetGene
   tbl.trena.both$target <- targetGene
   save(tbl.trena.tf, tbl.trena.both, file=sprintf("trena.%s.RData", targetGene))

} # build.models.with.inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
build.models.with.inclusiveFimoTable()
