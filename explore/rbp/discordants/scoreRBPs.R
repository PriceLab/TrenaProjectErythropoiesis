library(RUnit)
library(RPostgreSQL)
library(BSgenome.Hsapiens.UCSC.hg38)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
#----------------------------------------------------------------------------------------------------
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx <- get(load(file.path(data.dir, "mtx-rna-srm-27270x13.RData")))
dim(mtx)
early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")
#----------------------------------------------------------------------------------------------------
learn.discordant.strength.categories <- function()
{
   tbls <- list()
   proteins <- grep("p$", rownames(mtx), v=TRUE)
   checkEquals(length(proteins), 100)

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
       # round off the values for greater scrutability, no loss of relevant information
   tbl$cor.early <- round(tbl$cor.early, digits=2)
   tbl$cor.late <- round(tbl$cor.late, digits=2)
   tbl$cor.all <- round(tbl$cor.all, digits=2)
   tbl$corDelta <- round(tbl$corDelta, digits=2)

   fivenum(tbl$corDelta) # -1.69 -0.73 -0.05  0.58  1.69

      # for the present, we are only interested in gene/protein pairs which are
      # anti-correlated late in the time course.
      # include everything that is negative, even if small negative values will
      # be uninteresting

      # divide the discordant proteins into 4 groups
   sprintf("%5.2f", (fivenum(subset(tbl, corDelta < 0)$corDelta)))
      # "-1.69" "-1.11" "-0.71" "-0.51" "-0.00"

   divergent.late.strong1 <- subset(tbl, corDelta <= -1.11)$gene
   divergent.late.strong2 <- subset(tbl, corDelta > -1.11 & corDelta <= -0.71 )$gene
   divergent.late.strong3 <- subset(tbl, corDelta > -0.71 & corDelta <= -0.51 )$gene
   divergent.late.weak <- subset(tbl, corDelta > -0.51 & corDelta < 0)$gene



   length(divergent.late.strong1)  # 13
      # BACH1 CTCF  TCF3  E2F4  E2F8  GATA1 HLTF  KLF1  LDB1  SMC3  STAT2 TAL1  ZBTB7A
   length(divergent.late.strong2)  # 13
      # BCL11A FOXO3 GFI1B KMT2C NR2C2 PSIP1 RAD21 RXRB  STAT1 STAT3 SUZ12 TFDP1  TRIM33
   length(divergent.late.strong3)  # 13
      # CEBPB NR3C1 HCFC1 HDAC3 IKZF1 KLF13 KLF3 DR1 POU2F1 SCML2 SSRP1 USF1 WDHD1
   length(divergent.late.weak)     # 12
      # CREBBP CHD4 DOT1L FLI1 HDAC2 MTA1 NRF1 OGT SIRT6 SMARCB1 SPI1 KDM6A


   list(strong1=divergent.late.strong1,
        strong2=divergent.late.strong2,
        strong3=divergent.late.strong3,
        weak=divergent.late.weak,
        tbl.divergence=tbl)

} # learn.discordant.strength.categories
#----------------------------------------------------------------------------------------------------
assess.rbps <- function()
{
    x <- learn.discordant.strength.categories()

    divergent.late.strong1 <- sprintf("%sp", x$strong1)
    divergent.late.strong2 <- sprintf("%sp", x$strong2)
    divergent.late.strong3 <- sprintf("%sp", x$strong3)
    divergent.late.weak    <- sprintf("%sp", x$weak)
    divergent.late.all  <- c(divergent.late.strong1,
                             divergent.late.strong2,
                             divergent.late.strong3,
                             divergent.late.weak)
    tbl.models <- get(load("trena-92-ProteinModels-32603x10.RData"))
    dim(tbl.models)
    head(tbl.models)
    subset(tbl.models, gene=="DDX3X" & rank <=3 & target %in% divergent.late.strong1)

    as.data.frame(sort(
        table(subset(tbl.models, class=="rbp" &
                                 rank <= 5 &
                                 target %in% c(divergent.late.strong1, divergent.late.strong2))$gene),
        decreasing=TRUE))
   #   rank <= 5                   rank <= 3       rank <= 3, divergene.late.all   rank <= 5, all
   #       Var1 Freq                Var1 Freq                  Var1 Freq               Var1 Freq
   # 1    DDX3X    9          1    DDX3X    8            1    DDX3X   12         1    DDX3X   16
   # 2     FXR2    2          2    BUD13    1            2    BUD13    2         2    RBM47    3
   # 3    RBM47    2          3    CPSF1    1            3    DGCR8    2         3    BUD13    2
   # 4    BUD13    1          4    DGCR8    1            4    RBM47    2         4    DGCR8    2
   # 5    CPSF1    1          5   DICER1    1            5    CPSF1    1         5     FXR2    2
   # 6    DGCR8    1          6  IGF2BP3    1            6   DICER1    1         6     AARS    1
   # 7   DICER1    1          7    NCBP3    1            7  IGF2BP3    1         7    CPSF1    1
   # 8   DIS3L2    1          8    RBM47    1            8    NCBP3    1         8    DDX42    1
   # 9    EIF3D    1          9    SF3A3    1            9     PUM2    1         9   DICER1    1
   # 10   GPKOW    1          10   SRRM4    1            10   SF3A3    1         10  DIS3L2    1
   # 11 IGF2BP3    1          11  ZNF622    1            11   SRRM4    1         11   EIF3D    1
   # 12   NCBP3    1                                     12  ZNF622    1         12   GPKOW    1
   # 13    NPM1    1                                                             13  HNRNPC    1
   # 14   SF3A3    1                                                             14 IGF2BP3    1
   # 15   SRRM4    1                                                             15   NCBP3    1
   # 16  YTHDC1    1                                                             16    NPM1    1
   # 17  ZNF622    1                                                             17    PUM2    1
   #                                                                             18   SF3A3    1
   #                                                                             19   SRRM4    1
   #                                                                             20    YBX3    1
   #                                                                             21  YTHDC1    1
   #                                                                             22  ZNF622    1

    as.data.frame(sort(
        table(subset(tbl.models, class=="tf" &
                                 rank <= 5 &
                                 target %in% c(divergent.late.strong1, divergent.late.strong2))$gene),
        decreasing=TRUE))

} # assess.rbps
#----------------------------------------------------------------------------------------------------
# This led us to denote these RBPs as “multifunctional RBPs” - context-specific RBPs whose
# functional outcome depends on the set of mRNAs it targets. A key example appears to be the
# multifunctional RBP DDX3X (Mo et al., 2021; Soto-Rifo et al., 2012), whose abundance correlates
# significantly with the mRNA levels of 339 target genes (padj = 2.83 × 10−5; Glass’ ∆ = 6.9) and
# the translational efficiency of 730 target genes (padj = 5.25 × 10−5; Glass’ ∆= 11.89), of which
# only 43 targets overlap between both sets (Figure 3A and 3C). The consequences of DDX3X binding
# for mRNA abundance (positive correlation) or TE (negative correlation) are opposite, though this
# is not the case for all multifunctional RBPs (Figure S3B).

ddx3x.bindingSitesAndStructure <- function()
{
    subset(tbl.models, gene=="DDX3X" & rank <=3 & target %in% divergent.late.strong1)
      # LDB1 has some top scores: spearman, rf, ridge

    targetGenes <- c("CTCF", "E2F4", "KLF1", "LDB1", "SMC3", "TAL1", "ZBTB7A")
    targetGene <- targetGenes[1]

    rbp <- "DDX3X"
    # igv <- start.igv(targetGene, "hg38")
    showGenomicRegion(igv, targetGene)
    setBrowserWindowTitle(igv, targetGene)
    zoomOut(igv)
    query <- sprintf("select * from rbp where gene='%s' and target='%s'", rbp, targetGene)
    tbl.out <- dbGetQuery(geneRegDB, query)
    dim(tbl.out)
    tbl.out <- tbl.out[order(tbl.out$start, decreasing=FALSE),]
    tbl.out.wide <- subset(tbl.out, width >= 20 & celltype=="K562")
    dim(tbl.out.wide)
    rownames(tbl.out.wide) <- NULL

    tbl.track <- tbl.out.wide[, c("chrom", "start", "endpos", "score")]
    rownames(tbl.track) <- NULL
    track <- DataFrameQuantitativeTrack(sprintf("%s.wide", rbp), tbl.track, autoscale=TRUE, color="brown")
    displayTrack(igv, track)

    tbl.out.narrow <- subset(tbl.out, width < 20 & celltype=="K562")
    dim(tbl.out.narrow)
    track <- DataFrameQuantitativeTrack(sprintf("%s.narrow", rbp), tbl.out.narrow, autoscale=TRUE, color="brown")
    displayTrack(igv, track)


} # ddx3x.bindingSitesAndStructure
#----------------------------------------------------------------------------------------------------
# http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi
fold.rna <- function()
{
   goi <- getGenomicRegion(igv)
   as.character(with(goi, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end)))


   row <- 1
   as.character(with(tbl.out.wide[row,,drop=FALSE], getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, endpos)))

} # fold.rna
#----------------------------------------------------------------------------------------------------
determine.background.binding.rates <- function()
{
    query <- sprintf("select * from rbp where gene='%s' and celltype='K562'", rbp)
    system.time(tbl.out <- dbGetQuery(geneRegDB, query))

    dim(tbl.out)

    fivenum(tbl.out$score)  # 3.123918   6.024251  11.712333  25.737027 400.000000
    tbl.freq <- as.data.frame(sort(table(subset(tbl.out, width >= 20)$target), decreasing=TRUE))
    dim(tbl.freq)  # 9125 2

    tbl.freq <- as.data.frame(sort(table(subset(tbl.out)$target), decreasing=TRUE))
    dim(tbl.freq)  # 9125 2


    tbl.freq <- as.data.frame(sort(table(subset(tbl.out, width >= 20 & score > 5)$target), decreasing=TRUE))
    dim(tbl.freq)  # 9125 2

    query <- "select target from rbp where celltype='K562'"
    system.time(tbl.out <- dbGetQuery(geneRegDB, query))
    dim(tbl.out)
    length(unique(tbl.out$target))  #16168

} # determine.background.binding.rates
#----------------------------------------------------------------------------------------------------
