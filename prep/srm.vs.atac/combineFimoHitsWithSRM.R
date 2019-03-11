library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
x <- org.Hs.egCHRLENGTHS[1:24]
tbl.chromRegions <- data.frame(chrom=names(x), start=1, end=as.numeric(x), stringsAsFactors=FALSE)
tbl.chromRegions$chrom <- paste("chr", tbl.chromRegions$chrom, sep="")

load("~/github/TrenaProjectErythropoiesis/prep/import/srm-from-jeff/srm-copyNumberMatrices-rep1-rep2.RData")

# use two tfs for first look:
#  FLI1: possibly a pioneer, cor(rep1, rep2): 0.965
#   KLF: late stage, may depend on other proteins to open up the chromatin. cor(rep1, rep2): 0.832

stopifnot(dim(tbl.cor) == c(107, 2))
poi <- subset(tbl.cor, cor > 0.9)$protein
# choose FLI1: a suppressor of erythroid differentiation
motif.fli1 <- query(MotifDb, c("sapiens", "FLI1", "jaspar2018", "MA0475.1"))
export(motif.fli1, con="~/github/fimoService/pfms/fli1.meme", format="meme")
plot(mtx.rep1["FLI1",], mtx.rep2["FLI1",])

library(FimoClient)
library(BSgenome.Hsapiens.UCSC.hg38)
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)

load("~/github/TrenaProjectErythropoiesis/prep/import/srm-from-jeff/matrices.rep1.and.rep2.RData")
atac.files <- grep("narrowPeak$", list.files("../import/atacPeaks"), v=TRUE)

#-------------------------------------------------------------------------------------------------------------
getMatches <- function(tbl, threshold=10e-4)
{
   chroms <- tbl$chrom
   starts <- tbl$start
   ends   <- tbl$end

   printf("--- requesting %d sequences", nrow(tbl))
   seqs <- as.list(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, chroms, starts, ends)))
   names(seqs) <- sprintf("seq%06d", 1:length(seqs))
   length(seqs)
   printf("--- requesting matches in %d regions", nrow(tbl))
   tbl.out <- requestMatch(fc, seqs, threshold)
   invisible(tbl.out)

} # getMatches
#--------------------------------------------------------------------------------------------------------------
tbl.hits.genome <- getMatches(tbl.chromRegions[1,], 10e-5)

# d04
tbl.d04.1 <- read.table("../import/atacPeaks/ATAC_Cord_d04_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d04.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d04.1)
tbl.hits.d04.1 <- getMatches(tbl.d04.1, 10e-5)
printf("d04.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d04.1), nrow(tbl.d04.1), 100 * nrow(tbl.hits.d04.1)/ nrow(tbl.d04.1))


tbl.d04.2 <- read.table("../import/atacPeaks/ATAC_Cord_d04_rep2_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d04.2) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d04.2)
tbl.hits.d04.2 <- getMatches(tbl.d04.2, 10e-5)
printf("d04.2 matches: %d/%d: %5.2f", nrow(tbl.hits.d04.2), nrow(tbl.d04.2), 100 * nrow(tbl.hits.d04.2)/ nrow(tbl.d04.2))


tbl.d08.1 <- read.table("../import/atacPeaks/ATAC_Cord_d08_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d08.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d08.1)
tbl.hits.d08.1 <- getMatches(tbl.d08.1, 10e-5)
printf("d08.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d08.1), nrow(tbl.d08.1), 100 * nrow(tbl.hits.d08.1)/ nrow(tbl.d08.1))


tbl.d10.1 <- read.table("../import/atacPeaks/ATAC_Cord_d10_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d10.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d10.1)
tbl.hits.d10.1 <- getMatches(tbl.d10.1, 10e-5)
printf("d10.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d10.1), nrow(tbl.d10.1), 100 * nrow(tbl.hits.d10.1)/ nrow(tbl.d10.1))

tbl.d10.2 <- read.table("../import/atacPeaks/ATAC_Cord_d10_rep2_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d10.2) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d10.2)
tbl.hits.d10.2 <- getMatches(tbl.d10.2, 10e-5)
printf("d10.2 matches: %d/%d: %5.2f", nrow(tbl.hits.d10.2), nrow(tbl.d10.2), 100 * nrow(tbl.hits.d10.2)/ nrow(tbl.d10.2))

tbl.d11.1 <- read.table("../import/atacPeaks/ATAC_Cord_d11_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d11.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d11.1)
tbl.hits.d11.1 <- getMatches(tbl.d11.1, 10e-5)
printf("d11.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d11.1), nrow(tbl.d11.1), 100 * nrow(tbl.hits.d11.1)/ nrow(tbl.d11.1))

tbl.d11.2 <- read.table("../import/atacPeaks/ATAC_Cord_d11_rep2_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d11.2) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d11.2)
tbl.hits.d11.2 <- getMatches(tbl.d11.2, 10e-5)
printf("d11.2 matches: %d/%d: %5.2f", nrow(tbl.hits.d11.2), nrow(tbl.d11.2), 100 * nrow(tbl.hits.d11.2)/ nrow(tbl.d11.2))

tbl.d16.1 <- read.table("../import/atacPeaks/ATAC_Cord_d16_rep1_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d16.1) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d16.1)
tbl.hits.d16.1 <- getMatches(tbl.d16.1, 10e-5)
printf("d16.1 matches: %d/%d: %5.2f", nrow(tbl.hits.d16.1), nrow(tbl.d16.1), 100 * nrow(tbl.hits.d16.1)/ nrow(tbl.d16.1))

tbl.d16.2 <- read.table("../import/atacPeaks/ATAC_Cord_d16_rep2_hg38_macs2_default_peaks.narrowPeak", sep="\t", as.is=TRUE)
colnames(tbl.d16.2) <- c("chrom", "start", "end", "name", "score", "strand", "score2", "score3", "score4", "score5")
dim(tbl.d16.2)
tbl.hits.d16.2 <- getMatches(tbl.d16.2, 10e-5)
printf("d16.2 matches: %d/%d: %5.2f", nrow(tbl.hits.d16.2), nrow(tbl.d16.2), 100 * nrow(tbl.hits.d16.2)/ nrow(tbl.d16.2))

print(load("../import/srm-from-jeff/srm-copyNumberMatrices-rep1-rep2.RData"))

timepoints <- c("Day4",  "Day8", "Day10",  "Day11", "Day16")
count <- length(timepoints)


tbl.fli1.hits <- data.frame(atac.rep1=rep(0, count),
                            fli1.match1=rep(0, count),
                            atac.rep2=rep(0, count),
                            fli1.match2=rep(0, count),
                            fli1.srm1=rep(0, count),
                            fli1.srm2=rep(0, count),
                            stringsAsFactors=0)
rownames(tbl.fli1.hits) <- timepoints

tbl.fli1.hits["Day4", "atac.rep1"] <- nrow(tbl.d04.1)
tbl.fli1.hits["Day4", "fli1.match1"] <- nrow(tbl.hits.d04.1)
tbl.fli1.hits["Day4", "atac.rep2"] <- nrow(tbl.d04.2)
tbl.fli1.hits["Day4", "fli1.match2"] <- nrow(tbl.hits.d04.2)

tbl.fli1.hits["Day8", "atac.rep1"] <- nrow(tbl.d08.1)
tbl.fli1.hits["Day8", "fli1.match1"] <- nrow(tbl.hits.d08.1)

tbl.fli1.hits["Day10", "atac.rep1"] <- nrow(tbl.d10.1)
tbl.fli1.hits["Day10", "fli1.match1"] <- nrow(tbl.hits.d10.1)
tbl.fli1.hits["Day10", "atac.rep2"] <- nrow(tbl.d10.2)
tbl.fli1.hits["Day10", "fli1.match2"] <- nrow(tbl.hits.d10.2)

tbl.fli1.hits["Day11", "atac.rep1"] <- nrow(tbl.d11.1)
tbl.fli1.hits["Day11", "fli1.match1"] <- nrow(tbl.hits.d11.1)
tbl.fli1.hits["Day11", "atac.rep2"] <- nrow(tbl.d11.2)
tbl.fli1.hits["Day11", "fli1.match2"] <- nrow(tbl.hits.d11.2)

tbl.fli1.hits["Day16", "atac.rep1"] <- nrow(tbl.d16.1)
tbl.fli1.hits["Day16", "fli1.match1"] <- nrow(tbl.hits.d16.1)
tbl.fli1.hits["Day16", "atac.rep2"] <- nrow(tbl.d16.2)
tbl.fli1.hits["Day16", "fli1.match2"] <- nrow(tbl.hits.d16.2)

tbl.fli1.hits["Day4", "fli1.srm1"] <- mtx.cn1["FLI1", "Day4"]
tbl.fli1.hits["Day4", "fli1.srm2"] <- mtx.cn2["FLI1", "Day4"]

tbl.fli1.hits["Day8", "fli1.srm1"] <- mtx.cn1["FLI1", "Day8"]
tbl.fli1.hits["Day8", "fli1.srm2"] <- mtx.cn2["FLI1", "Day8"]

tbl.fli1.hits["Day10", "fli1.srm1"] <- mtx.cn1["FLI1", "Day10"]
tbl.fli1.hits["Day10", "fli1.srm2"] <- mtx.cn2["FLI1", "Day10"]

tbl.fli1.hits["Day11", "fli1.srm1"] <- mtx.cn1["FLI1", "Day11"]
tbl.fli1.hits["Day11", "fli1.srm2"] <- mtx.cn2["FLI1", "Day11"]

tbl.fli1.hits["Day16", "fli1.srm1"] <- mtx.cn1["FLI1", "Day14"]   # note: atac-seq but no srm at day16. substituting day14
tbl.fli1.hits["Day16", "fli1.srm2"] <- mtx.cn2["FLI1", "Day14"]


# tbl.fli1.hits
#       atac.rep1 fli1.match1 atac.rep2 fli1.match2 fli1.srm1 fli1.srm2
# Day4     305762       32305    133208       18000      3929  1282.478
# Day8     192979       31820         0           0       690     0.000
# Day10     84179       16875    116255       21485       334     0.000
# Day11     83692       15898     98661       18091       328     0.000
# Day16     90002       20620     58240       16477         0     0.000




