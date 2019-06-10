library(TrenaProjectErythropoiesis)
library(FimoClient)
library(MotifDb)

if(!exists("tp"))
   tp <- TrenaProjectErythropoiesis()
if(!exists("tbl.regions")){
   setTargetGene(tp, "GATA2")
   tbl.enhancers <- getEnhancers(tp)
   tbl.gata2 <- getTranscriptsTable(tp)
   tss <- tbl.gata2$tss
   chrom <- tbl.gata2$chrom
   loc.min <- min(tbl.enhancers$start)
   loc.max <- max(tbl.enhancers$end)
   printf("query fimo across %d bases", loc.max - loc.min)
   tbl.regions <- data.frame(chrom="chr3", start=loc.min, end=loc.max, stringsAsFactors=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
findBindingSites <- function(tbl.regions, threshold=1e-2)
{
    FIMO_HOST <- "localhost"
    FIMO_PORT <- 60000

   if(!exists("fc")){
      fc <<- FimoClient(FIMO_HOST, FIMO_PORT)
      }

    tbl.matches <- requestMatchForRegions(fc, tbl.regions, "hg38", threshold)
    motif.names <- tbl.matches$motif
    deleters <- setdiff(motif.names, names(MotifDb))
    xyz <- "delete?"
    for(deleter in deleters){
       deleter.indices <- grep(deleter, tbl.matches$motif)
       if(length(deleter.indices) > 0)
          tbl.matches <- tbl.matches[-deleter.indices,]
       }
    tfs <- mcols(MotifDb[tbl.matches$motif])$geneSymbol
    tbl.matches$tf <- tfs

    tbl.matches

} # findBindingSites
#------------------------------------------------------------------------------------------------------------------------
tbl.fimo <- findBindingSites(tbl.regions, threshold=1e-2)
save(tbl.fimo, file="tbl.fimo.RData")


