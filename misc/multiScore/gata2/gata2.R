source("~/github/regulatoryGenomePaper/demos/common.R")
if(!exists("igv")){
  igv <- igvR()
  setGenome(igv, "hg38")
  showGenomicRegion(igv, "GATA2")
  }
if(!exists("tpe")){
   tpe <- TrenaProjectErythropoiesis()
   setTargetGene(tpe, "GATA2")
   }

targetGene.info <- as.list(getTranscriptsTable(tpe)[c("tss", "strand")])


if(!exists("mtx")){
   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   }

tbl.gh <- getEnhancersForTissues("all", TRUE)
with(tbl.gh, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom[1], min(start) - 1000, max(end) + 1000)))
tbl.fp <- get(load("tbl.gata2.d12.fp.RData"))
track <- DataFrameQuantitativeTrack("fp", tbl.fp[, c("chrom", "start", "end", "score")], autoscale=FALSE,
                                    min=0, max=max(tbl.fp$score), color="random")
displayTrack(igv, track)


meme.file <- "jaspar2018.meme"
motifs <- query(MotifDb, c("sapiens", "jaspar2018"))
length(motifs)
export(motifs, con=meme.file, format="meme")
tbl.fimo <- fimoBatch(tbl.fp[, c("chrom", "start", "end")], matchThreshold=1e-6, genomeName="hg38", pwmFile=meme.file)
dim(tbl.fimo)

with(tbl.fimo, plot(score, -log10(p.value)))
tbl.fimo$fimo <- round(-log10(tbl.fimo$p.value), digits=2)
track <- DataFrameAnnotationTrack("fimo", tbl.fimo[,c("chrom", "start", "end", "tf", "fimo")], col="red",
                                  displayMode="EXPANDED", trackHeight=200)
displayTrack(igv, track)
tbl.cons <- as.data.frame(gscores(phast.7, GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
tbl.fimo$phast7 <- round(tbl.cons$default, digits=2)

tbl.cons <- as.data.frame(gscores(phast.100, GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
tbl.fimo$phast100 <- round(tbl.cons$default, digits=2)

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




