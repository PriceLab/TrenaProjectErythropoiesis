library(igvR)
library(MotifDb)
library(TrenaProjectErythropoiesis)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38

source("~/github/fimoService/batchMode/fimoBatchTools.R")
#------------------------------------------------------------------------------------------------------------------------
conservationTrack <- function()
{
   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)
   tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.blocks)), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

} # conservationTrack
#------------------------------------------------------------------------------------------------------------------------


if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, "GATA2")
}

tfs <- c("TAL1", "LYL1", "LMO2", "GATA2", "RUNX1", "MEIS1", "SPI1", "FLI1", "GFI1B", "ERG")
motifs <- query(MotifDb, c("hsapiens"), tfs[10])
length(motifs)
names(motifs)
displayMotifs(motifs, 3)

gata2 <- c("Hsapiens-jaspar2018-GATA2-MA0036.1",
           "Hsapiens-jaspar2018-GATA2-MA0036.2",
           "Hsapiens-jaspar2018-GATA2-MA0036.3")

runx1 <- "Hsapiens-jaspar2018-RUNX1-MA0002.1"

meis1 <- c("Hsapiens-HOCOMOCOv10-MEIS1_HUMAN.H10MO.C",
           "Hsapiens-jaspar2018-MEIS1-MA0498.2")

spi1 <- c("Hsapiens-HOCOMOCOv10-SPI1_HUMAN.H10MO.A",
          "Hsapiens-jaspar2018-SPI1-MA0080.1")

fli1 <- c("Hsapiens-jaspar2018-FLI1-MA0475.2")

gfi1b <- c("Hsapiens-HOCOMOCOv10-GFI1B_HUMAN.H10MO.C")

erg <-  "Hsapiens-jaspar2018-ERG-MA0474.2"

motif.names <- c(gata2, runx1, meis1, spi1, fli1, gfi1b, erg)
meme.file <- "tmp.meme"
export(MotifDb[motif.names], con=meme.file, format="meme")

tpe <- TrenaProjectErythropoiesis()
targetGene <- "GATA2"
setTargetGene(tpe, targetGene)
tbl.geneInfo <- getTranscriptsTable(tpe)
tbl.regions <- with(tbl.geneInfo, data.frame(chrom=chrom, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE))
with(tbl.regions, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end)))


conservationTrack()
tbl.fimo <- fimoBatch(tbl.regions, matchThreshold=1e-4, genomeName="hg38", pwmFile=meme.file)
tbl.fimo$pvalScore <- -log10(tbl.fimo$p.value)
dim(tbl.fimo)

for(TF in tbl.fimo$tf){
   tbl.sub <- subset(tbl.fimo, tf==TF)[, c("chrom", "start", "end", "pvalScore")]
   track <- DataFrameQuantitativeTrack(TF, tbl.sub, autoscale=TRUE, color="random")
   displayTrack(igv, track)
   }


