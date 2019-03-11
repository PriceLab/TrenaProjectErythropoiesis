# atac-seq, srm, jeff & marjorie: only 5 shared time points: days 4, 8*, 10, 11, 12, FLI1 and KLF1

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
   tbl.hits <- findFimoHitsInAtacSeq(tbl.test, 1e-4)
   checkTrue(nrow(tbl.hits) > 75)

} # test_findFimoHitsInAtacSeqRegions
#--------------------------------------------------------------------------------------------------------------
# FIMO_HOST <- "khaleesi"
# FIMO_PORT <- 5558
# fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)
# requestMatch(fc, list(seq0="ATCGATTTTGGGGGAAAAAATTTT"), 10e-2)

# getFimoMatches <- function(day, replicate,  threshold=10e-4)
# {
#    tbl.day.rep <- read.table("../import/atacPeaks/ATAC_Cord_d04_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
# colnames(tbl.d04.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
# dim(tbl.d04.1)
# tbl.hits.d04.1 <- getMatches(tbl.d04.1, 10e-5)
# printf("d04.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d04.1), nrow(tbl.d04.1), 100 * nrow(tbl.hits.d04.1)/ nrow(tbl.d04.1))
# 
#    chroms <- tbl$chrom
#    starts <- tbl$start
#    ends   <- tbl$end
# 
#    printf("--- requesting %d sequences", nrow(tbl))
#    seqs <- as.list(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, chroms, starts, ends)))
#    names(seqs) <- sprintf("seq%06d", 1:length(seqs))
#    length(seqs)
#    printf("--- requesting matches in %d regions", nrow(tbl))
#    tbl.out <- requestMatch(fc, seqs, threshold)
#    invisible(tbl.out)
# 
# } # getFimoMatches
# #--------------------------------------------------------------------------------------------------------------
# test_getFimoMatches <- function()
# {
#    printf("--- test_getFimoMatches")
#    getFimoMatches(day=4, rep=1)
# 
# } # test_getFimoMatches
# #--------------------------------------------------------------------------------------------------------------
# 
# # load("~/github/TrenaProjectErythropoiesis/prep/import/srm-from-jeff/srm-copyNumberMatrices-rep1-rep2.RData") # mtx.cn1, mtx.cn2
# # dim(mtx.cn1)
# # mtx.cn1[1:10, 1:10]
# 
# 
# 
