source("~/github/regulatoryGenomePaper/demos/common.R")

if(!exists("tbl.atac")){
   load("~/github/TrenaProjectErythropoiesis/misc/multiScore/brandAtacSeqCuration/tbl.atac.RData")
   printf("brand tbl.atac loaded, %d x %d", nrow(tbl.atac), ncol(tbl.atac))
   }

if(!exists("igv")){
  igv <- igvR()
  setGenome(igv, "hg38")
  }
if(!exists("tpe")){
   tpe <- TrenaProjectErythropoiesis()
   }

if(!exists("mtx")){
   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   }

#------------------------------------------------------------------------------------------------------------------------
# 1) find tss, chrom, strand, min and max of genehancer regions (aka, "geneRegion")
# 2) load atac-seq regions in geneRegion
# 3) find all TFBS in atac-seq regions
# 4) score all TFBS for conservation
# 5) determine distance from TSS for all TFBS
build.multiScoreTable <- function(targetGene, display=FALSE)
{
  if(display)
     showGenomicRegion(igv, targetGene)

  setTargetGene(tpe, targetGene)
  targetGene.info <- as.list(getTranscriptsTable(tpe)[c("tss", "strand", "chrom", "start", "end")])
  tbl.gh <- getEnhancersForTissues("all", TRUE)
  gh.loc.start <- min(tbl.gh$start) - 1000
  gh.loc.end   <- max(tbl.gh$end) + 1000
  gh.chrom <- tbl.gh$chrom[1]
  showGenomicRegion(igv, sprintf("%s:%d-%d", gh.chrom, gh.loc.start, gh.loc.end))

  tbl.fp <- subset(tbl.atac, chrom==gh.chrom & start > gh.loc.start & end < gh.loc.end)
  dim(tbl.fp)

  track <- DataFrameQuantitativeTrack("fp", tbl.fp[, c("chrom", "start", "end", "meanScore", "sdScore")],
                                      autoscale=FALSE, min=0, max=max(tbl.fp$meanScore), color="random")
  displayTrack(igv, track)

  meme.file <- "jaspar2018.meme"
  motifs <- query(MotifDb, c("sapiens", "jaspar2018"))
  length(motifs)
  export(motifs, con=meme.file, format="meme")
  #tbl.geneRegion <- with(targetGene.info, data.frame(chrom=chrom, start=start-10000,
  #                                                   end=end+ 10000, stringsAsFactors=FALSE))
  #with(tbl.geneRegion, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end)))
  tbl.fimo <- fimoBatch(tbl.fp, matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file)
  dim(tbl.fimo)

  with(tbl.fimo, plot(score, -log10(p.value)))
  tbl.fimo$fimo <- round(-log10(tbl.fimo$p.value), digits=2)
  track <- DataFrameAnnotationTrack("fimo", tbl.fimo[,c("chrom", "start", "end", "tf", "fimo")], col="red",
                                    displayMode="EXPANDED", trackHeight=200)
  displayTrack(igv, track)
  track <- DataFrameAnnotationTrack("tbx15", subset(tbl.fimo, tf=="TBX15")[,c("chrom", "start", "end", "tf", "fimo")],
                                    col="blue", trackHeight=30, displayMode="EXPANDED")
  displayTrack(igv, track)

  tbl.oi <- subset(tbl.fimo, tf=="TBX15")
  dim(tbl.oi)
  tbl.cons <- as.data.frame(gscores(phast.7, GRanges(tbl.oi[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
  dim(tbl.cons)
  colnames(tbl.cons)[1] <- "chrom"
  tbl.cons$chrom <- as.character(tbl.cons$chrom)

  tbl.track <- tbl.cons[, c(1,2,3,6)]
  track <- DataFrameQuantitativeTrack("phast7", tbl.track, autoscale=FALSE, min=0, max=1.0, color="random")
  displayTrack(igv, track)

  tbl.fimo$phast7 <- round(tbl.cons$default, digits=2)

  tbl.cons <- as.data.frame(gscores(phast.100, GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
  tbl.fimo$phast100 <- round(tbl.cons$default, digits=2)
  colnames(tbl.cons)[1] <- "chrom"
  tbl.cons$chrom <- as.character(tbl.cons$chrom)
  tbl.track <- tbl.cons[, c(1,2,3,6)]

  track <- DataFrameQuantitativeTrack("phast100", tbl.track, autoscale=FALSE, min=0, max=1.0, color="random")
  displayTrack(igv, track)

  tbl.cons <- as.data.frame(gscores(phast.30, GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
  tbl.fimo$phast30 <- round(tbl.cons$default, digits=2)
  colnames(tbl.cons)[1] <- "chrom"
  tbl.cons$chrom <- as.character(tbl.cons$chrom)
  tbl.track <- tbl.cons[, c(1,2,3,6)]
  track <- DataFrameQuantitativeTrack("phast30", tbl.track, autoscale=FALSE, min=0, max=1.0, color="random")
  displayTrack(igv, track)


  cor <- unlist(lapply(tbl.fimo$tf, function(tf) if(tf %in% rownames(mtx)) return(cor(mtx["GATA2",], mtx[tf,])) else return(NA)))
  cor <- round(cor, digits=2)
  tbl.fimo$cor <- cor
  tbl.fimo$cor.abs <- abs(cor)
  tbl.fimo$target <- "GATA2"
  tbl.fimo$width <- 1 + tbl.fimo$end - tbl.fimo$start
  coi <- c("chrom", "start", "end", "width", "tf", "target", "cor", "cor.abs", "fimo", "phast7", "phast100")

  tbl <- tbl.fimo[, coi]
  dim(tbl)

  tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl), GRanges(tbl.gh)))
  colnames(tbl.ov) <- c("fimo", "gh")
  tbl$gh <- rep(0, nrow(tbl))
  tbl$gh[tbl.ov$fimo] <- round(tbl.gh$combinedscore[tbl.ov$gh], digits=2)

  tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl), GRanges(tbl.fp)))
  colnames(tbl.ov) <- c("fimo", "fp")
  tbl$fp <- rep(0, nrow(tbl))
  tbl$fp[tbl.ov$fimo] <- tbl.fp$score[tbl.ov$fp]

  tbl$to.tss <- (targetGene.info$tss - tbl$start)  * (targetGene.info$strand) * -1
  rownames(tbl) <- NULL

} # buildMultiScoreTable
#------------------------------------------------------------------------------------------------------------------------
gata2.tbx15.tfbs.for.marjorie <- function()
{

  # these are scored with bioc pwmmatch.  fimo rates them at 1.25e-4, 2.08e-4, 1.67e-4, 2.08e-4
  # which, by mariam's work, are good enough values
  # 1: chr3:128,483,003-128,483,009    Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128483002 128483009      +   6.642737          0.8767398 AGGGGTGA
  # 2: chr3:128,483,452-128,483,457    Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128483451 128483458      +   6.658030          0.8787582 CGGTGTGA
  # 3: chr3:128,486,956-128,486,962    Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128486956 128486962      -   ?
  # 4: chr3:128,497,528-128,497,534    Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128497527 128497534      +   6.658030          0.8787582 CGGTGTGA


   tbl.track <- data.frame(chrom=c("chr3", "chr3", "chr3", "chr3"),
                           start=c(128483003, 128483452, 128486956, 128497528),
                           end=c(128483009, 128483457, 128486962, 128497534),
                           stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("tbx15-old", tbl.track, color="random")
   displayTrack(igv, track)



} # gata2.tbx15.tfbs.for.marjorie
#------------------------------------------------------------------------------------------------------------------------



