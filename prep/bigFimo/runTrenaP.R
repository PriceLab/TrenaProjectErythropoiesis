args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
targetGene <- args[1]
rnaGene <- sub("p$", "", targetGene)
library(RUnit)

source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
f <- sprintf("../%s/tbl.fimo.%s.RData", rnaGene, rnaGene)
f
file.exists(f)
tbl.fimo <- get(load(f))

dim(tbl.fimo)                 # 114021 9
trenaProject <- TrenaProjectErythropoiesis()
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx.rna <- get(load(file.path(data.dir, "mtx-rna-pkfm-27170x14.RData")))
mtx.srm <- get(load(file.path(data.dir, "mtx-srm-copyNumber-100x13.RData")))
mtx.both <- get(load(file.path(data.dir, "mtx-rna-srm-2720x13.RData")))

early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")

mtx.both <- mtx.both[, late]
dim(mtx.both)

#----------------------------------------------------------------------------------------------------
build.models.with.inclusiveFimoTable <- function()
{
   message(sprintf("--- test_inclusiveFimoTable"))

   tms <- TMS$new(trenaProject, rnaGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.both, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   tfs <- unique(subset(tbl.tms, fimo_pvalue < 1e-3 & abs(cor.all) > 0.3 &
                                 (gh > 10 | chip | oc))$tf)
   printf("candidate tf count: %d", length(tfs))

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.both, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)
   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))

      #------------------------------------------------------------
      # build the protein model, late time points only
      #------------------------------------------------------------
   tbl.trena.p <- tms$build.trena.model(tfs, rbps, mtx.both, alternateTarget=targetGene,
                                        order.by="rfScore")
   tbl.trena.p <- tbl.trena.p[order(abs(tbl.trena.p$pearsonCoeff), decreasing=TRUE),]
   tbl.trena.p$rank <- seq_len(nrow(tbl.trena.p))
   rownames(tbl.trena.p) <- NULL
   tbl.trena.p$target <- targetGene
   for(col in 2:7)
      tbl.trena.p[, col] <- round(tbl.trena.p[, col], digits=3)

   print(head(tbl.trena.p, n=10))
   save(tbl.trena.p, file=sprintf("trena.p.lateTimepoints.%s.RData", targetGene))

      #------------------------------------------------------------
      # build the gene (rna) model, late time points only
      #------------------------------------------------------------
   tbl.trena.g <- tms$build.trena.model(tfs, rbps, mtx.both, alternateTarget=rnaGene,
                                        order.by="rfScore")
   tbl.trena.g <- tbl.trena.g[order(abs(tbl.trena.g$pearsonCoeff), decreasing=TRUE),]
   rownames(tbl.trena.g) <- NULL
   tbl.trena.g$target <- rnaGene
   tbl.trena.g$rank <- seq_len(nrow(tbl.trena.g))
   for(col in 2:7)
      tbl.trena.g[, col] <- round(tbl.trena.g[, col], digits=3)

   print(head(tbl.trena.g, n=12))
   save(tbl.trena.g, file=sprintf("trena.g.lateTimepoints.%s.RData", rnaGene))


} # build.models.with.inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
if(!interactive())
    build.models.with.inclusiveFimoTable()
