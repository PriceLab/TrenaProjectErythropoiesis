library(RUnit)
file <- "tbl.trena-47-targets.RData"
tbl.models <- get(load(file))
dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx.rna <- get(load(file.path(dir, "mtx-rna-pkfm-27170x14.RData")))[,-14]
mtx.srm <- get(load(file.path(dir, "mtx-srm-copyNumber-100x13.RData")))

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_assignRank()
   test_assignRankAllTables()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
assignRank <- function(tbl, scoreName)
{
   stopifnot(length(unique(tbl$target)) == 1)
   stopifnot(scoreName %in% colnames(tbl))

      # for instance, for OGT, the largest rf value is in row 95
      # head(tbl$rf, n=3):  95  85  91
      # tbl$rfScore[head(new.order, n=3)]  146.5121 117.2359 103.3986

   new.order <- tbl[, scoreName]
   tbl.reordered <- tbl[new.order,]
   tbl.reordered$rank <- seq_len(nrow(tbl.reordered))

   tbl.reordered

} # assignRank
#------------------------------------------------------------------------------------------------------------------------
test_assignRank <- function()
{
    message(sprintf("--- test_assignRandk"))

   tbl.sub <- subset(tbl.models, target=="OGT")
   dim(tbl.sub)  # 117 16

   tbl.ranked <- assignRank(tbl.sub, "rf")
   checkEquals(nrow(tbl.ranked), nrow(tbl.sub))
   checkEquals(tbl.ranked$gene[1:3], c("MAZ", "NR2C2", "NCBP3"))
   checkEqualsNumeric(tbl.ranked$rfScore[1:3], c(146.5121, 117.2359, 103.3986), tol=0.001)

   tbl.ranked <- assignRank(tbl.sub, "pear.abs")
   checkEquals(nrow(tbl.ranked), nrow(tbl.sub))
   checkEquals(tbl.ranked$gene[1:3], c("CSTF2T", "ELAVL1", "UCHL5"))
   checkEqualsNumeric(tbl.ranked$pearsonCoeff[1:3], c(-0.9153902, -0.9117892, -0.9094098), tol=0.001)

   tbl.ranked <- assignRank(tbl.sub, "spear.pos")
   checkEquals(nrow(tbl.ranked), nrow(tbl.sub))
   checkEquals(tbl.ranked$gene[1:3], c("PBX2", "STAU1", "AKAP8L"))
   checkEqualsNumeric(tbl.ranked$spearmanCoeff[1:3], c(0.8769231, 0.8505495, 0.7978022), tol=0.001)

} # test_assignRank
#------------------------------------------------------------------------------------------------------------------------
assignRankAllTables <- function(tbl, scoreName)
{
   targetGenes <- unique(tbl$target)
   tbls <- list()
   for(targetGene in targetGenes){
      tbl.sub <- subset(tbl, target==targetGene)
      tbl.ranked <- assignRank(tbl.sub, scoreName)
      tbls[[targetGene]] <- tbl.ranked
      }
   tbl.all.ranked <- do.call(rbind, tbls)
   rownames(tbl.all.ranked) <- NULL
   invisible(tbl.all.ranked)

} # assignRankAllTables
#------------------------------------------------------------------------------------------------------------------------
test_assignRankAllTables <- function()
{
   message(sprintf("--- test_assignRankAllTables"))
   tbl.ranked <- assignRankAllTables(tbl.models, "spear.pos")
   checkEquals(nrow(tbl.models), nrow(tbl.ranked))
   checkEquals(subset(tbl.ranked, gene=="STAU1" & rank <=10)$target, c("BACH1", "OGT", "RCOR1"))

   tbl.ranked <- assignRankAllTables(tbl.models, "rf")
   x <- as.list(head(sort(table(subset(tbl.ranked, rank <= 5)$gene), decreasing=TRUE), n=3))
   x <- x[sort(names(x))]
   checkEquals(x, list(FOSL1=5, HINFP=4, KLF9=5))

} # test_assignRankAllTables
#------------------------------------------------------------------------------------------------------------------------
combine <- function()
{

    tbls <- list()
    i <- 0

    for(file in files){
      i <- i + 1
      printf("file: %s", file)
      load(file)
      tbl <- tbl.trena.both
      new.order <- order(abs(tbl$spearmanCoeff), decreasing=TRUE)
      tbl <- tbl[new.order,]
      if("rank" %in% colnames(tbl)){
          rank.col <- grep("rank", colnames(tbl))
          tbl <- tbl[, -rank.col]
      }
      tbl$spear.pos <- order(tbl$spearmanCoeff, decreasing=TRUE)
      tbl$spear.abs  <- order(abs(tbl$spearmanCoeff), decreasing=TRUE)
      tbl$pear.pos  <- order(tbl$pearsonCoeff, decreasing=TRUE)
      tbl$pear.abs  <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
      tbl$rf <- order(tbl$rfScore, decreasing=TRUE)
      tbl$lasso.pos <- order(tbl$betaLasso, decreasing=TRUE)
      tbl$lasso.abs <- order(abs(tbl$betaLasso), decreasing=TRUE)
      # if(tbl$target[1]=="OGT") browser()
      tbls[[i]] <- tbl
      }

    tbl <- do.call(rbind, tbls)
    dim(tbl)

    tbl.ogt <- subset(tbl, target=="OGT")
    sort(table(subset(tbl, rank <= 10 & class=="rbp")$gene), decreasing=TRUE)

} # combine
#------------------------------------------------------------------------------------------------------------------------
top.regulators <- function()
{
   max <- 10

   tbl.ranked <- assignRankAllTables(tbl.models, "rf")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=3)

   tbl.ranked <- assignRankAllTables(tbl.models, "spear.pos")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "spear.abs")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "pear.pos")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "pear.abs")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "lasso.abs")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "lasso.pos")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)

   tbl.ranked <- assignRankAllTables(tbl.models, "rf")
   head(as.data.frame(sort(table(subset(tbl.ranked, rank <= max)$gene), decreasing=TRUE)), n=10)



} # top.6
#----------------------------------------------------------------------------------------------------
# marjorie:  The ProEB to EB transition is roughly from Day 10 to Day 14.
#
# More precisely (based on cytof) it is from Day 12 to Day 14, but a
# lot of proteins start to change at Day 11. So I think it is fair
# to look from Day10 to 14.
find.discordant.genes <- function()
{
  source("~/s/data/public/human/symToGeneID.R"); test_assignGeneIDs()

  proteins <- rownames(mtx.srm)
  length(proteins)
  #length(intersect(
  setdiff(proteins, rownames(mtx.rna))
  no.gene <- setdiff(proteins, rownames(mtx.rna))
  length(no.gene)
  x <- assignGeneIDs(no.gene)

  map <- c("BCL11A_XL_L" = "BCL11A",
           "CBP" = "CREBBP",
           "coREST" = "RCOR1",
           "E12/E47" =  "TCF3",
           "GR" = "NR3C1",
           "HXB4" = "HOXB4",
           "MLL3" = "KMT2C",
           "MLL4 (KMT2D)" = "KMT2D",
           "NC2B" = "DR1",
           "PO2F1" = "POU2F1",
           "RPB1" = "POLR2A",
           "SET1B" = "SETD1B",
           "SETB1" = "SETDB1",
           "SIR6" = "SIRT6",
           "SMCA4" = "SMARCA4",
           "SMRC1" = "SMARCC1",
           "SNF5" = "SMARCB1",
           "STA5A" = "STAT5A",
           "T2FA" = "GTF2F1",
           "TF2B" = "GTF2F2",
           "TF3C2" = "GTF3C2",
           "UBF1" = "UBTF",
           "UTX"  = "KDM6A",
           "ZC11A" =  "ZC3H11A")

   all(names(map) %in% rownames(mtx.srm))
   all(as.character(map) %in% rownames(mtx.rna))
   indices <- match(names(map), rownames(mtx.srm))

   rownames(mtx.srm)[indices] <- as.character(map)
   dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"

   save(mtx.srm, file=file.path(dir, "mtx-srm-copyNumber-100x13.RData"))

} # find.discordant.genes
#----------------------------------------------------------------------------------------------------
gpkow <- function()
{
    subset(tbl, gene=="GPKOW" & spear.abs < 10)$target  # "RAD21"  "RCOR1"  "SAP130" "TAL1"
    gpkow <- mtx.rna["GPKOW",]
    gpkow <- gpkow/max(gpkow)
    plot(gpkow, type="b", col="black", ylim=c(0,1))

    rad21 <- mtx.rna["RAD21",]
    rad21 <- rad21/max(rad21)
    lines(rad21, type="b", col="blue")

    rcor1 <- mtx.rna["RCOR1",]
    rcor1 <- rcor1/max(rcor1)
    lines(rcor1, type="b", col="green")

    sap130 <- mtx.rna["SAP130",]
    sap130 <- sap130/max(sap130)
    lines(sap130, type="b", col="brown")

    tal1 <- mtx.rna["TAL1",]
    tal1 <- tal1/max(tal1)
    lines(tal1, type="b", col="red")

    tal1.p <- mtx.srm["TAL1",]
    tal1.p <- tal1.p/max(tal1.p)
    lines(tal1.p, type="b", col="blue")

    legend(1, 0.98, c("GPKOW mRNA", "TAL1 mRNA", "TAL1 protein"),
           c("black", "red",       "blue")
           )

} # gpkow
#----------------------------------------------------------------------------------------------------
