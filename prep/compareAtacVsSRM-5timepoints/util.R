# atac-seq, srm, jeff & marjorie: only 5 shared time points: days 4, 8*, 10, 11, 12, FLI1 and KLF1
#--------------------------------------------------------------------------------------------------------------
# use these to form the columns of interest (coi) in the assembled big matrix
days <- c(4, 4, 8, 10, 10, 11, 11, 12, 12)
reps <- c(1, 2, 1, 1,   2,  1,  2,  1,  2)
#--------------------------------------------------------------------------------------------------------------
library(RUnit)
library(FimoClient)
library(BSgenome.Hsapiens.UCSC.hg38)

FIMO_HOST <- "khaleesi"
FIMO_PORT <- 5558
if(!exists("fc"))
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
#--------------------------------------------------------------------------------------------------------------
createPfmsFileNeededForFIMO <- function()
{
   # from manual inspection:
   #   for KLF1, load up fimoServer with Hsapiens-SwissRegulon-KLF1.SwissRegulon
   #   for FLI1, Hsapiens-jaspar2018-FLI1-MA0475.1

   library(MotifDb)
   pfm.names <- c("Hsapiens-SwissRegulon-KLF1.SwissRegulon", "Hsapiens-jaspar2018-FLI1-MA0475.1")
   moi <- MotifDb[pfm.names]
   export(moi, con="~/github/fimoService/pfms/fli1.klf1.human.meme", format="meme")
   # consensusString(moi[[1]]) "?AGGGTG?GGC"  -> "ATCGATCGAAGGGTGAGGCATCGATCG"
   # consensusString(moi[[2]]) "?CAGGAAGTGG"  -> "ACTACAGGAAGTGGCAGCAGCAGCAGCAG"


} # createPfmsFileNeededForFIMO
#--------------------------------------------------------------------------------------------------------------
test_fimoService <- function()
{
   # start or restart fimoServer:
   #   cd ~/github/fimoService/server 
   #   python -i runFimoServer.py 5558 "/users/pshannon/meme/bin/fimo" "../pfms/fli1.klf1.human.meme"
   #
  tbl <- requestMatch(fc, list(klf1="ATCGATCGAAGGGTGAGGCATCGATCG", fli1="ACTACAGGAAGTGGCAGCAGCAGCAGCAG"), 1e-4)
  checkEquals(dim(tbl), c(2, 9))

} # test_fimoService
#--------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_readAtacSeq()
  test_findFimoHitsInAtacSeqRegions()
  test_lineCountAtacSeq()
  test_reshapeAtacSeqData()
  test_extractHitData()
   
} # runTests
#--------------------------------------------------------------------------------------------------------------
readAtacSeq <- function(dayNumber, replicateNumber)
{
   dir <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   atac.files <- grep("narrowPeak$", list.files(dir), v=TRUE)
   dayString <- sprintf("d%02d", dayNumber)
   repString <- sprintf("rep%1d", replicateNumber)
   searchString <- sprintf("%s_%s", dayString, repString)

   stopifnot(any(grepl(searchString, atac.files)))

   foi <- grep(searchString, atac.files, value=TRUE)
   full.path <- file.path(dir, foi)
   stopifnot(file.exists(full.path))

   tbl <- read.table(full.path, sep="\t", as.is=TRUE)
     # from http://ginolhac.github.io/chip-seq/peak/
     # chromosome
     # start
     # end
     # peak name
     # integer score for display
     # strand
     # fold-change
     # -log10 pvalue
     # -log10 qvalue
     # relative summit position to peak start
   macs2.colnames <-  c("chrom", "start", "end", "name", "displayScore", "strand", "fc", "pvalScore", "qvalScore", "peakOffset")
   colnames(tbl) <- macs2.colnames
   invisible(tbl)

} # readAtacSeq
#--------------------------------------------------------------------------------------------------------------
test_readAtacSeq <- function()
{
   printf("--- test_readAtacSeq")

   tbl.4.1 <- readAtacSeq(4, 1)

   expected <-  c("chrom", "start", "end", "name", "displayScore", "strand", "fc", "pvalScore", "qvalScore", "peakOffset")
   checkEquals(colnames(tbl.4.1), expected)
   checkEquals(dim(tbl.4.1), c(305762, 10))

   tbl.12.2 <- readAtacSeq(12, 2)
   checkEquals(colnames(tbl.12.2), expected)
   checkEquals(dim(tbl.12.2), c(107772, 10))

} # test_readAtacSeq
#--------------------------------------------------------------------------------------------------------------
# when sanity checking the open chromatin regions in each file, the count of those regions should
# roughly correlate with the number of lines in the file.   that number is found here.
lineCountAtacSeq <- function(dayNumber, replicateNumber)
{
   dir <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   atac.files <- grep("narrowPeak$", list.files(dir), v=TRUE)
   dayString <- sprintf("d%02d", dayNumber)
   repString <- sprintf("rep%1d", replicateNumber)
   searchString <- sprintf("%s_%s", dayString, repString)

   stopifnot(any(grepl(searchString, atac.files)))

   foi <- grep(searchString, atac.files, value=TRUE)
   full.path <- file.path(dir, foi)
   stopifnot(file.exists(full.path))
   cmd <- sprintf("wc -l %s", full.path)
   x <- system(cmd, intern=TRUE)
   as.numeric(strsplit(x, " ")[[1]][1])

} # lineCountAtacSeq
#--------------------------------------------------------------------------------------------------------------
test_lineCountAtacSeq <- function()
{
   printf("--- test_lineCountAtacSeq")

   checkEquals(lineCountAtacSeq(4,1), 305762)
   checkEquals(lineCountAtacSeq(4,2), 133208)
           
} # test_lineCountAtacSeq
#--------------------------------------------------------------------------------------------------------------
# for all TFs currently loaded into the FimoServer
findFimoHitsInAtacSeqRegions <- function(tbl.atac, threshold)
{
   requestMatchForRegions(fc, tbl.atac[, 1:3], "hg38", 1e-4)

} # findFimoHitsInAtacSeqRegions
#--------------------------------------------------------------------------------------------------------------
test_findFimoHitsInAtacSeqRegions <- function()
{
   printf("--- test_findFimoHitsInAtacSeqRegions")

   tbl.4.1 <- readAtacSeq(4, 1)
   dim(tbl.4.1)
   tbl.test <- tbl.4.1[1:100,]
   tbl.hits <- findFimoHitsInAtacSeqRegions(tbl.test, 1e-4)
   checkTrue(nrow(tbl.hits) > 75)

} # test_findFimoHitsInAtacSeqRegions
#--------------------------------------------------------------------------------------------------------------
reshapeAtacSeqData <- function()
{
   f <- "../import/srm-from-jeff/srm-copyNumberMatrices-rep1-rep2.RData"
   printf(load(f))
   stopifnot(all(rownames(mtx.cn1) == rownames(mtx.cn2)))
   stopifnot(all(colnames(mtx.cn1) == colnames(mtx.cn2)))

   coi <- sprintf("day.%02d.%1d", days, reps)
   mtx <- matrix(0, nrow=nrow(mtx.cn1), ncol=length(coi),
                 dimnames=list(rownames(mtx.cn1), coi))
    #"day.04.1" "day.04.2" "day.08.1" "day.10.1" "day.10.1" "day.11.1" "day.11.2" "day.12.1" "day.12.2"
   mtx[, "day.04.1"] <- mtx.cn1[, "Day4"]
   mtx[, "day.04.2"] <- mtx.cn2[, "Day4"]
     #mtx[, "day.08.1"] <- (mtx.cn1[, "Day8"] + mtx.cn2[, "Day8"])/2
   mtx[, "day.08.1"] <- mtx.cn1[, "Day8"]
   mtx[, "day.10.1"] <- mtx.cn1[, "Day10"]
   mtx[, "day.10.2"] <- mtx.cn2[, "Day10"]
   mtx[, "day.11.1"] <- mtx.cn1[, "Day11"]
   mtx[, "day.11.2"] <- mtx.cn2[, "Day11"]
   mtx[, "day.12.1"] <- mtx.cn1[, "Day12"]
   mtx[, "day.12.2"] <- mtx.cn2[, "Day12"]
   
   invisible(mtx)   

} # reshapeAtacSeqData
#--------------------------------------------------------------------------------------------------------------
test_reshapeAtacSeqData <- function()
{
   mtx <- reshapeAtacSeqData()
   checkEquals(dim(mtx), c(107, 9))

} # test_reshapeAtacSeqData
#--------------------------------------------------------------------------------------------------------------
extractHitData <- function(regions, hits, day, rep, geneSymbol, mode)
{
   name <- sprintf("day.%d.%1d", day, rep)
   tbl <- hits[[name]]
   tbl.sub <- tbl[grep(geneSymbol, tbl$motif),]
   hit.count.total <- nrow(tbl.sub)
   topQuartileScore <- fivenum(tbl.sub$score)[4]
   hit.count.topScore <- nrow(subset(tbl.sub, score >= topQuartileScore))
   
   tbl.region <- regions[[name]]
   region.count <- nrow(tbl.region)

   result <- switch(mode,
                "total" = hit.count.total,
                "top"   = hit.count.topScore,
                "total-density" = {round(100 * hit.count.total/region.count, digits=2)},
                "top-density" = {round(100* hit.count.topScore/region.count, digits=2)}
                )
   result
   
} # extractHitData
#--------------------------------------------------------------------------------------------------------------
test_extractHitData <- function()
{
   printf("--- test_extractHitData")
   gene <- "KLF1"

   tbl.sub <- hits[[1]][grep("KLF1", hits[[1]]$motif),]
   expected <- nrow(tbl.sub)
   total.klf1.hits.day.4.1 <- extractHitData(regions, hits, days[1], reps[1], gene, "total")
   checkEquals(expected, total.klf1.hits.day.4.1)

   topQuartile <- fivenum(tbl.sub$score)[4]
   tbl.sub.sub <- subset(tbl.sub, score >= topQuartile)
   expected <- nrow(tbl.sub.sub)
   top.klf1.hits.day.4.1 <- extractHitData(regions, hits, days[1], reps[1], gene, "top")
   checkEquals(expected, top.klf1.hits.day.4.1)

      # density is the percentage of motif matches/open regions
   total.density.klf1.hits.day.4.1 <- extractHitData(regions, hits, days[1], reps[1], gene, "total-density")
   top.density.klf1.hits.day.4.1 <- extractHitData(regions, hits, days[1], reps[1], gene, "top-density")

      # we expect the ratio of top/total hits to be the same as topDensity/totalDensity
      # since the same denominator is used in each case: the count of open chromatin regions
      # identified by ChIP-seq
   ratio.top.hits.total.hits <- top.klf1.hits.day.4.1/total.klf1.hits.day.4.1
   ratio.top.dentisty.total.density <- top.density.klf1.hits.day.4.1/total.density.klf1.hits.day.4.1

   checkEqualsNumeric(ratio.top.hits.total.hits, ratio.top.density.total.density, tol=0.001)

} # test_extractHitData
#--------------------------------------------------------------------------------------------------------------
