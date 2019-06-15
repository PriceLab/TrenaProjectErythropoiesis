#------------------------------------------------------------------------------------------------------------------------
getATACseq.byRegion <- function(chromosome, start.loc, end.loc)
{
   directory <- "../../../prep/import/atacPeaks"
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

} # getATACseq.byRegion
#------------------------------------------------------------------------------------------------------------------------
to.RData.files <- function()
{
   directory <- "../../../prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      save(tbl.atac, file=sprintf("ATAC_brand_%s.RData", track.name))
      } # for file



} # to.RData.files
#------------------------------------------------------------------------------------------------------------------------
