if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 1)
   targetProtein <- args[1]
   }



#data.set <- "rna.srm"
data.set <- "fpkm.tmt"
#data.set <- "asinh-rna.tmt"
#------------------------------------------------------------
#
#------------------------------------------------------------
if(data.set == "rna.srm"){
   targetProtein <- "KLF1p"
   rnaGene <- sub("p$", "", targetProtein)
   data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
   mtx.both <- get(load(file.path(data.dir, "mtx-rna-srm-27270x13.RData")))
   late <- c("D8", "D8_5", "D10", "D10_5", "D11", "D11_5", "D12", "D14")
   mtx.both <- mtx.both[, late]
   }

if(data.set == "fpkm.tmt"){
   targetProtein <- "KLF1.Q13351"
   rnaGene <- sub("\\..*$", "", targetProtein)
   data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression"
   mtx.both <- get(load(file.path(data.dir, "mtx.fpkm.tmt.31760x6.RData")))
   late  <- c("d8", "d10", "d12", "d14")
   mtx.both <- mtx.both[, late]
   }

if(data.set == "asinh-rna.tmt"){
   targetProtein <- "KLF1.Q13351"
   rnaGene <- sub("\\..*$", "", targetProtein)
   data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression"
   mtx.both <- get(load(file.path(data.dir, "mtx.rna.tmt.31761x6.RData")))
   late  <- c("d8", "d10", "d12", "d14")
   mtx.both <- mtx.both[, late]
   }

#----------------------------------------------------------------------------------------------------
library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#----------------------------------------------------------------------------------------------------

tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
f <- sprintf("../../../prep/bigFimo/%s/tbl.fimo.%s.RData", rnaGene, rnaGene)
file.exists(f)
tbl.fimo <- get(load(f))

dim(tbl.fimo)                 # 114021 9
trenaProject <- TrenaProjectErythropoiesis()

    #---------------------------------------------------------------------
    # ignore RBP binding sites, use only expression levels, mean and sd
    #---------------------------------------------------------------------
data.dir <- "~/github/rnaBindingProteins/explore/vonNostrand"
tbl.rbp <-
    read.table(file.path(data.dir, "RBPsFromVonNostrand_2020.txt"), sep="\t",
               skip=2, header=FALSE, nrows=-1, fill=TRUE)[, 1:2]
colnames(tbl.rbp) <- c("geneSymbol", "ensg")
dim(tbl.rbp)
head(tbl.rbp)

rbp.candidates <- intersect(tbl.rbp$geneSymbol, rownames(mtx.both))

mean(mtx.both)  # 213
rbp.means <- as.integer(apply(mtx.both[rbp.candidates,], 1, mean))
rbp.sd    <- as.integer(apply(mtx.both[rbp.candidates,], 1, sd))

tbl.rbp.stats <- data.frame(mean=rbp.means, sd=rbp.sd, gene=rbp.candidates,
                            row.names=rbp.candidates, stringsAsFactors=FALSE)

mean.threshold <- 50
sd.threshold   <- 5
tbl.rbp.filtered <- subset(tbl.rbp.stats, mean >= 200 & sd >= sd.threshold)
dim(tbl.rbp.filtered)

rbp.candidates.filtered <- rownames(tbl.rbp.filtered)
length(rbp.candidates.filtered)

#----------------------------------------------------------------------------------------------------
regulators.from.filtered.matrix <- function(mtx, mean.cutoff, sd.cutoff)
{

   means <- as.integer(apply(mtx, 1, mean))
   sds   <- as.integer(apply(mtx, 1, sd))
   genes <- rownames(mtx)
   tbl.stats <- data.frame(mean=means, sd=sds, gene=genes,
                          row.names=genes, stringsAsFactors=FALSE)

   tbl.filtered <- subset(tbl.stats, mean >= mean.cutoff & sd >= sd.cutoff)
   return(tbl.filtered$gene)

} # filter.matrix
#----------------------------------------------------------------------------------------------------
build.models.with.inclusiveFimoTable <- function()
{
   message(sprintf("--- test_inclusiveFimoTable"))

   mean.cutoff <- 50
   sd.cutoff <- 10

   tms <- TMS$new(trenaProject, rnaGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$addRBP()
   #tms$addRBP(tbl.rbp.filtered)
   tms$add.rbp.mrna.correlations(mtx.both, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)
   rbps <- unique(tbl.rbp$gene)
   rbps <- regulators.from.filtered.matrix(mtx.both[rbps,], mean.cutoff, sd.cutoff)
   printf("candidate rbps: %d", length(rbps))

   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.both, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   tfs <- unique(subset(tbl.tms, fimo_pvalue < 1e-3 &  & (gh > 10 | chip | oc))$tf)

   #tf.means <- as.integer(apply(mtx.both[tfs,], 1, mean))
   #tf.sd    <- as.integer(apply(mtx.both[tfs,], 1, sd))
   #tbl.tf.stats <- data.frame(mean=tf.means, sd=tf.sd, gene=tfs, stringsAsFactors=FALSE)
   #dim(tbl.tf.stats)
   #tbl.tf.stats.filtered <- subset(tbl.tf.stats, mean >= mean.threshold & sd >= sd.threshold)
   #dim(tbl.tf.stats.filtered)  # 19
   #tfs <- unique(tbl.tf.stats.filtered$gene)
   length(tfs)
   tfs <- regulators.from.filtered.matrix(mtx.both[tfs,], mean.cutoff, sd.cutoff)
   length(tfs)

   printf("candidate tf count: %d", length(tfs))


      #------------------------------------------------------------
      # build the protein model, late time points only
      #------------------------------------------------------------
   tbl.trena.p <- tms$build.trena.model(c(tfs, rbps),
                                        alternateTarget=targetProtein,
                                        mtx.both, order.by="rfScore")
   tbl.trena.p <- tbl.trena.p[order(abs(tbl.trena.p$pearsonCoeff), decreasing=TRUE),]
   tbl.trena.p$rank <- seq_len(nrow(tbl.trena.p))
   rownames(tbl.trena.p) <- NULL
   tbl.trena.p$target <- targetProtein
   for(col in 2:7)
      tbl.trena.p[, col] <- round(tbl.trena.p[, col], digits=3)

   printf("------------ tbl.trena.p, data.set %s, target: %s",
          data.set, targetProtein)
   print(head(tbl.trena.p, n=12))
   table(tbl.trena.p$class)
   browser()

   #save(tbl.trena.p, file=sprintf("trena.p.lateTimepoints.%s.RData", targetProtein))

      #------------------------------------------------------------
      # build the gene (rna) model, late time points only
      #------------------------------------------------------------
   tbl.trena.g <- tms$build.trena.model(c(tfs, rbps), mtx.both, alternateTarget=rnaGene,
                                        order.by="rfScore")
   tbl.trena.g <- tbl.trena.g[order(abs(tbl.trena.g$pearsonCoeff), decreasing=TRUE),]
   rownames(tbl.trena.g) <- NULL
   tbl.trena.g$rank <- seq_len(nrow(tbl.trena.g))
   tbl.trena.g$target <- rnaGene
   for(col in 2:7)
      tbl.trena.g[, col] <- round(tbl.trena.g[, col], digits=3)

   printf("------------ tbl.trena.g, data.set %s, target: %s",
          data.set, rnaGene)
   print(head(tbl.trena.g, n=12))
   browser()
   #save(tbl.trena.g, file=sprintf("trena.g.lateTimepoints.%s.RData", rnaGene))


} # build.models.with.inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
if(!interactive())
    build.models.with.inclusiveFimoTable()
