library(RUnit)
#----------------------------------------------------------------------------------------------------
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx <- get(load(file.path(data.dir, "mtx-rna-srm-27270x13.RData")))
dim(mtx)
proteins <- grep("p$", rownames(mtx), v=TRUE)
checkEquals(length(proteins), 100)
early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")
#----------------------------------------------------------------------------------------------------
learn.discordant.strength.categories <- function()
{
   tbls <- list()
   for(protein in proteins){
       gene <- sub("p$", "", protein)

       gene.late <- mtx[gene, late]
       prot.late <- mtx[protein, late]

       gene.early <- mtx[gene, early]
       prot.early<- mtx[protein, early]

       gene.all <- mtx[gene,]
       prot.all <- mtx[protein,]

       cor.early <- cor(gene.early, prot.early, method="spearman", use="pairwise.complete")
       cor.late <- cor(gene.late, prot.late, method="spearman", use="pairwise.complete")
       cor.all <- cor(gene.all, prot.all, method="spearman", use="pairwise.complete")

       tbl <- data.frame(gene=gene, protein=protein,
                         cor.early = cor.early,
                         cor.late = cor.late,
                         cor.all = cor.all,
                         stringsAsFactors=FALSE)
       tbls[[protein]] <- tbl
       } # for protein

   tbl <- do.call(rbind, tbls)
   rownames(tbl) <- NULL
   dim(tbl)
   tbl$corDelta <- tbl$cor.late - tbl$cor.early
   fivenum(tbl$corDelta) # -1.6857143 -0.7285714 -0.0500000  0.5785714  1.6928571
     # we are only interested in gene/protein pairs which are anti-correlated late
     # include everything that is negative, even if small negative values will
     # be uninteresting
   tbl.negCor <- subset(tbl, corDelta < 0)

   tbl$cor.early <- round(tbl$cor.early, digits=2)
   tbl$cor.late <- round(tbl$cor.late, digits=2)
   tbl$cor.all <- round(tbl$cor.all, digits=2)
   tbl$corDelta <- round(tbl$corDelta, digits=2)

      # divide the discordant proteins into 4 groups
   sprintf("%5.2f", (fivenum(subset(tbl, corDelta < 0)$corDelta)))
      # "-1.69" "-1.11" "-0.71" "-0.51" "-0.00"

   divergent.late.strong4 <- subset(tbl, corDelta <= -1.11)$gene
   divergent.late.strong3 <- subset(tbl, corDelta > -1.11 & corDelta <= -0.71 )$gene
   divergent.late.strong2 <- subset(tbl, corDelta > -0.71 & corDelta <= -0.51 )$gene
   divergent.late.weak <- subset(tbl, corDelta > -0.51 & corDelta <= -0.3)$gene

   length(divergent.late.strong4)
   length(divergent.late.strong3)
   length(divergent.late.strong2)
   length(divergent.late.weak)

   #save(tbl.negCor,
   #     divergent.late.strong4,
   #     divergent.late.strong3,
   #     divergent.late.strong2,
   #     divergent.late.weak,
   #     file="gene.protein.divergence.categories.RData")

} # learn.discordant.strength.categories
#----------------------------------------------------------------------------------------------------

