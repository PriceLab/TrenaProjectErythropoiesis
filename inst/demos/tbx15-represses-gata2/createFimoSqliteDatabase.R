library(TrenaProjectErythropoiesis)
library(FimoClient)
library(MotifDb)
library(RSQLite)
#------------------------------------------------------------------------------------------------------------------------
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
run <- function()
{
   tbl.regions.gata2.enhancer <- data.frame(chrom="chr3", start=128481976, end=128500133, stringsAsFactors=FALSE)
   tbl.fimo <- findBindingSites(tbl.regions.gata2.enhancer, threshold=2e-4)
   save(tbl.fimo, file="tbl.fimo.gata2.enhancer.19k.RData")

   tbl.regions.gata2.enhancer.214k <- data.frame(chrom="chr3", start=128407127, end=128621364, stringsAsFactors=FALSE)
   tbl.fimo <- findBindingSites(tbl.regions.gata2.enhancer.214k, threshold=2e-4)
     # dim(tbl.fimo)  # [1] 302522     10
     # length(grep("TBX15", tbl.fimo$tf)) 2663   10
   save(tbl.fimo, file="tbl.regions.gata2.enhancer.214k.RData")
   db <- dbConnect(SQLite(), "fimoResults-2e-4-chr3-128074944-128620958-214k.sqlite")
   dbWriteTable(db, name="fimoBindingSites", value=tbl.fimo, overwrite=TRUE)
   dbDisconnect(db)

     # try it out
   db <- dbConnect(SQLite(), "fimoResults-2e-4-chr3-128074944-128620958-214k.sqlite")
   tbl.atac.538bp <- dbGetQuery(db, "select * from fimoBindingSites where chrom='chr3' and start >= 128497516 and end <= 128498053")
   dim(tbl.atac.538bp)  # 1748 x 10


} # run
#------------------------------------------------------------------------------------------------------------------------


