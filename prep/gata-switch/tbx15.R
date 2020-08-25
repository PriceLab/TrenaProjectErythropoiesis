library(TrenaProjectErythropoiesis)
library(igvR)
library(trenaSGM)
library(FimoClient)
library(RUnit)
library(AnnotationHub)
source(file.path(Sys.getenv("HOME"), "github", "fimoService", "batchMode", "fimoBatchTools.R"))
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
epigenetic.markers <- c("H3K4me3", "H3K27me3", "H3K36me3", "H3k4me1", "H3K9me3", "H3K27ac")

#------------------------------------------------------------------------------------------------------------------------
library (RColorBrewer)
totalColorCount <- 12
#colors <- brewer.pal(12, "Paired")
colors <- brewer.pal(8, "Dark2")
currentColorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
if(!exists("parseChromLocString"))
   source("~/github/trena/R/utils.R")
if(!exists("tbl.geneInfo"))
   tbl.geneInfo <- get(load((system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))))
#------------------------------------------------------------------------------------------------------------------------
required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                           "motifScore", "motifRelativeScore", "match",
                                           "distance.from.tss", "tf")
#------------------------------------------------------------------------------------------------------------------------
mtx <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/mtx-rna.RData"))
mtx <- asinh(mtx)
tpe <- TrenaProjectErythropoiesis()
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   }
#------------------------------------------------------------------------------------------------------------------------
displayGeneHancer <- function(gene="GATA2")
{
   setTargetGene(tpe, gene)
   tbl.enhancers <- getEnhancers(tpe)
   track <- DataFrameQuantitativeTrack("GH", tbl.enhancers[, c("chrom", "start", "end", "combinedscore")],
                                       "brown", autoscale=FALSE, min=0, max=100)
   displayTrack(igv, track)
   with(tbl.enhancers, showGenomicRegion(igv, sprintf("%s:%d-%d", unique(chrom), min(start)-1000, max(end) + 1000)))

} # displayGeneHancer
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function(chromosome, start.loc, end.loc)
{
   directory <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
test_getATACseq <- function()
{
   printf("--- test_getATACseq")
   tbl.atac <- getATACseq("chr3", 128495142, 128498398)
   samples <- unique(tbl.atac$sample)
    # "d08_rep1" "d10_rep1" "d10_rep2" "d11_rep1" "d11_rep2" "d12_rep1" "d12_rep2" "d16_rep1" "d16_rep2"
   for(sample.x in samples){
      tbl.sample <- subset(tbl.atac, sample==sample.x)[, c("chrom", "start", "end")]
      write.table(tbl.sample, file=sprintf("tbl.%s.bed", sample.x), quote=FALSE, row.names=FALSE)
      } # for sample.x

} # test_getATACseq
#------------------------------------------------------------------------------------------------------------------------
displayATACseq <- function(tbl.all)
{
   samples <- unique(tbl.all$sample)
   current.day.string <- ""

   for(current.sample in samples){
      this.day.string <- strsplit(current.sample, "_")[[1]][1]
      if(this.day.string != current.day.string){
         currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
         color <- colors[currentColorNumber]
         printf("new day  %s   new color: %s", this.day.string, color)
         current.day.string <- this.day.string
         }
      tbl.atac.sub <- subset(tbl.all, sample == current.sample)
      track.name <- current.sample
      track <- DataFrameQuantitativeTrack(track.name, tbl.atac.sub[, c("chrom", "start", "end", "c10")],
                                          color, autoscale=FALSE, min=0, max=430)
      displayTrack(igv, track)
      } # for samples


   tbl.regions.condensed <- as.data.frame(union(GRanges(tbl.all[, c("chrom", "start", "end")]),
                                                GRanges(tbl.all[, c("chrom", "start", "end")])))[, c("seqnames", "start", "end")]
   colnames(tbl.regions.condensed) <- c("chrom", "start", "end")
   tbl.regions.condensed$chrom <- as.character(tbl.regions.condensed$chrom)
   lapply(tbl.regions.condensed, class)

   state$tbl.regions.condensed <- tbl.regions.condensed
   track <- DataFrameAnnotationTrack("atac combined", tbl.regions.condensed, color="black")
   displayTrack(igv, track)


} # displayATACseq
#------------------------------------------------------------------------------------------------------------------------
test_displayATACseq <- function()
{
  tbl.all <- getATACseq("chr3", 128459896, 128522959)
  displayATACseq(tbl.all)
  state$tbl.all <- tbl.all

} # test_displayATACseq
#------------------------------------------------------------------------------------------------------------------------
old.findTFs <- function(tbl.regions, threshold=1e-4)
{
    FIMO_HOST <- "khaleesi"
    FIMO_PORT <- 60000

   if(!exists("fc")){
      fc <<- FimoClient(FIMO_HOST, FIMO_PORT)
      }

    tbl.matches <- requestMatchForRegions(fc, tbl.regions, "hg38", threshold)
    motif.names <- tbl.matches$motif
    deleters <- setdiff(motif.names, names(MotifDb))
    for(deleter in deleters){
       #browser()
       deleter.indices <- grep(deleter, tbl.matches$motif)
       if(length(deleter.indices) > 0)
          tbl.matches <- tbl.matches[-deleter.indices,]
       }
    tfs <- mcols(MotifDb[tbl.matches$motif])$geneSymbol
    tbl.matches$tf <- tfs

    tbl.matches

} # old.findTFs
#------------------------------------------------------------------------------------------------------------------------
findTFs <- function(tbl.regions, fimo.threshold=1e-4, name, display, memeFile)
{

   tbl.motifMatches <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold,
                                 genomeName="hg38", pwmFile=memeFile)

   if(nrow(tbl.motifMatches) > 0){
     if(display){
        track <- DataFrameQuantitativeTrack(name,
                                            tbl.motifMatches[, c("chrom", "start", "end", "score")],
                                            autoscale=TRUE, color="blue")
        displayTrack(igv, track)
        } # if display
     } # if nrow

   invisible(tbl.motifMatches)

} # findTFs
#------------------------------------------------------------------------------------------------------------------------
test_findTFs <- function()
{
   printf("--- test_findTFs")
   tbl.regions <- data.frame(chrom="chr3", start=128482906, end=128483594, stringsAsFactors=TRUE)

   tbl.1 <- findTFs(tbl.regions, fimo.threshold=1e-3)
   dim(tbl.1)
   tbl.2 <- findTFs(tbl.regions, fimo.threshold=1e-6)

   track <- DataFrameQuantitativeTrack("fimo 1e-3", tbl.1[, c(1,2,3,6)], autoscale=TRUE)
   displayTrack(igv, track)


   tbl.bindingSites <- fixReg(tbl.bindingSites, "GATA2")
   state$tbl.bindingSites <- tbl.bindingSites

} # test_findTFs
#------------------------------------------------------------------------------------------------------------------------
displayBindingSites <- function()
{
   motifs <- grep("TBX15", tbl.bs$motifName, v=TRUE)
   table(motifs)  # yikes!
                  # Hsapiens-HOCOMOCOv10-TBX15_HUMAN.H10MO.D 99
                  # Hsapiens-jaspar2018-TBX15-MA0803.1    1
                  # Hsapiens-SwissRegulon-TBX15.SwissRegulon  1

   motifs <- unique(motifs)
   for(this.motif in motifs){
      tbl.motif <- subset(tbl.bs, motifName==this.motif & motifScore > 6)
      printf("%40s: %d", this.motif, nrow(tbl.motif))
      currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
      color <- colors[currentColorNumber]
      track <- DataFrameQuantitativeTrack(this.motif,
                                          tbl.motif[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                          autoscale=FALSE, min=0, max=30, color=color)
      displayTrack(igv, track)
      }

} # displayBindingSites
#------------------------------------------------------------------------------------------------------------------------
displayAllHocomocoBindingSitesInGenerousRegion <- function()
{
   tbl.hocomoc.tbx15 <- findTFs(data.frame(chrom="chr3", start=128394282, end=128587177, stringsAsFactors=FALSE))
   state$tbl.hocomoc.tbx15 <- tbl.hocomoc.tbx15
   save(state, file="tbx15.state.RData")

   tbl.tbx15 <- tbl.hocomoc.tbx15[grep("TBX15", tbl.hocomoc.tbx15$motif),]
   dim(tbl.tbx15)
   motifs <- tbl.tbx15$motif
   table(motifs)  # 1389, 88, 88
   tbl.tbx15 <- fixReg(tbl.tbx15, "GATA2")
   state$tbl.tbx15 <- tbl.tbx15

   motifs <- unique(motifs)
   length(motifs)

   for(this.motif in motifs){
      tbl.motif <- subset(tbl.tbx15, motifName==this.motif & motifScore > 6)
      dim(tbl.motif)
      currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
      color <- colors[currentColorNumber]
      trackName <- sprintf("big.%s", this.motif)
      track <- DataFrameQuantitativeTrack(trackName,
                                          tbl.motif[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                          autoscale=FALSE, min=0, max=30, color=color)
      displayTrack(igv, track)
      }


} # displayAllHocomocoBindingSitesInGenerousRegion
#------------------------------------------------------------------------------------------------------------------------
findLowerFidelityBindingSites <- function()
{
   pfm <- query(MotifDb, c("sapiens", "jaspar2018", "TBX15"))
   mm <- MotifMatcher("hg38", as.list(pfm), quiet=TRUE)
   tbl.matches <- findMatchesByChromosomalRegion(mm, state$tbl.regions.condensed, 85)
   dim(tbl.matches)
   trackName <- names(pfm)
   track <- DataFrameQuantitativeTrack(trackName,
                                       tbl.matches[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                       autoscale=FALSE, min=0.5, max=1.0, color="magenta")
   displayTrack(igv, track)


} # findLowerFidelityBindingSites
#------------------------------------------------------------------------------------------------------------------------
lookForChIPseq <- function()
{
    library(RPostgreSQL)
      # chipseq, the remap database
    db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="khaleesi")
    dbGetQuery(db, "select count(*) from chipseq")  # 80,143,591
    colnames(dbGetQuery(db, "select * from chipseq limit 3"))
      # "chrom"     "start"     "endpos"    "tf"        "name"      "strand"    "peakStart" "peakEnd"
   roi <- "chr3:128,471,624-128,502,021"
   tbl <- dbGetQuery(db, "select * from chipseq where chrom='chr3' and start >= 128471624 and endpos <= 128502021")
   dim(tbl)
   tbl.freq <- as.data.frame(table(tbl$tf))
   colnames(tbl.freq) <- c("tf", "count")
   tbl.freq <- tbl.freq[order(tbl.freq$count, decreasing=TRUE),]

    # look for chip-atlas
    dbGetQuery(db, "select * from pg_database")
    grep("^chip", dbGetQuery(db, "select datname from pg_database")$datname, value=TRUE)

    db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="chipatlas", host="khaleesi")
    dbListTables(db)  # experiments peaaks
    dbGetQuery(db, "select * from experiments limit 3")
    dbGetQuery(db, "select * from experiments limit 3")[, c("antigen", "cellType")]

} # lookForChIPseq
#------------------------------------------------------------------------------------------------------------------------
displayMethylationMarks <- function()
{
   print(load("h3k4me3.gr.RData"))  # created ~/github/trena/documents/gata2/tbx15.R
   track <- GRangesQuantitativeTrack("h3k4me3.gr", gr.gata2, autoscale=TRUE, color="green")
   displayTrack(igv, track)


} # displayMethylationMarks
#------------------------------------------------------------------------------------------------------------------------
addMethylationTrack <- function()
{
  loc <- getGenomicRegion(igv)
  saved.data.file <- "gr.hg38.generous.gata2.RData"

  if(!file.exists(saved.data.file)){
     aHub <- AnnotationHub()
     key <- "AH36576"
     x <- aHub[[key]]
     xi <- import(x)
     system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
     system("gunzip hg19ToHg38.over.chain.gz")
     chain <- import.chain("hg19ToHg38.over.chain")
     x.hg38 <- liftOver(xi, chain)
     gr.hg38 <- unlist(x.hg38)
     seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
     gr.hg38.gata2.region <- gr.hg38[chrom(gr.hg38)==loc$chrom &
                                     start(gr.hg38) >= loc$start &
                                     end(gr.hg38)   <= loc$end]
     length(gr.hg38.gata2.region)  # 13279
     save(gr.hg38.gata2.region, file=saved.data.file)
     } else {
        load(saved.data.file)
        }

  track <- GRangesQuantitativeTrack("h3k27ac", gr.hg38.gata2.region, autoscale=TRUE, color="green")
  displayTrack(igv, track)

} # addMethylationTrack
#------------------------------------------------------------------------------------------------------------------------
displayEpigeneticTrack <- function(markerName, loc)
{
   printf("---- marker: %s", markerName)


  gr.hg38 <- get(load(sprintf("gr.%s.RData", markerName)))

  gr.sub <- gr.hg38[chrom(gr.hg38)==loc$chrom &
                    start(gr.hg38) >= loc$start &
                    end(gr.hg38)   <= loc$end]

  track <- GRangesQuantitativeTrack(markerName, gr.sub, autoscale=TRUE, color="green")
  displayTrack(igv, track)

} # displayEpigeneticTrack
#------------------------------------------------------------------------------------------------------------------------
displayAtacSeqBamFiles <- function()
{
   bamFiles <- sort(list.files("gata2-bamFiles")) # "ATAC_Cord_d10_rep1_hg38_chr3_gata2_sorted.bam")
   names <- sub("_hg38.*$", "", sub("ATAC_Cord_", "", bamFiles))
   for(i in seq_len(length(bamFiles))){
      bamFile <- file.path("gata2-bamFiles", bamFiles[i])
      trackName <- sprintf("bam-%s", names[i])
      x <- readGAlignments(bamFile, use.names=TRUE) #, param=param)
      track <- GenomicAlignmentTrack(trackName, x, trackHeight=50, visibilityWindow=1000000, color="random")
      displayTrack(igv, track)
      }

} # displayAtacSeqBamFiles
#------------------------------------------------------------------------------------------------------------------------
#   print(system.time(tbl.match <- fimoBatch(tbl.np, matchThreshold=fimo.threshold, genomeName="hg19",
#                                            pwmFile="ctcf-human.meme")))
reproduce <- function()
{
  displayGeneHancer(gene="GATA2")
  loc <- getGenomicRegion(igv)
  tbl.all <- getATACseq(loc$chrom, loc$start, loc$end)
  displayATACseq(tbl.all)
  # state$tbl.all <- tbl.all

  #displayBindingSites()

  tbl.regions <- data.frame(chrom=loc$chrom, start=loc$start, end=loc$end, stringsAsFactors=FALSE)

  addMethylationTrack()

  loc <- getGenomicRegion(igv)

  for(marker in epigenetic.markers){
     displayEpigeneticTrack(marker, loc)
     }

  pfms.hocomoco <- MotifDb::query(MotifDb, c("sapiens", "tbx15", "hocomoco"))
  pfms.jaspar2018 <- MotifDb::query(MotifDb, c("sapiens", "tbx15", "jaspar2018"))

  export(pfms.hocomoco, "tbx15-hocomoco-human.meme", 'meme')
  export(pfms.jaspar2018, "tbx15-jaspar2018-human.meme", 'meme')

  tbl.hocomoco <- findTFs(tbl.regions, 1e-4, "TBX15:hocomco", TRUE, "tbx15-hocomoco-human.meme")
  tbl.jaspar   <- findTFs(tbl.regions, 1e-3, "TBX15:jaspar:1e-3",  TRUE, "tbx15-jaspar2018-human.meme")

} # reproduce
#------------------------------------------------------------------------------------------------------------------------
