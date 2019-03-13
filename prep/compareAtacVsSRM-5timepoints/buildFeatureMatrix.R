# atac-seq, srm, jeff & marjorie: only 5 shared time points: days 4, 8*, 10, 11, 12, FLI1 and KLF1
# call FLI1 and KLF1 hits in all ATAC-seq in the 5 shared time points, 4 of which have replicates
#--------------------------------------------------------------------------------------------------------------
if(!exists("readAtacSeq"))
    source("util.R")
#--------------------------------------------------------------------------------------------------------------
# atac-seq, srm, jeff & marjorie: only 5 shared time points: days 4, 8*, 10, 11, 12, FLI1 and KLF1
# call FLI1 and KLF1 hits in all ATAC-seq in the 5 shared time points, 4 of which have replicates
#--------------------------------------------------------------------------------------------------------------
if(!exists("readAtacSeq"))
    source("util.R")
#--------------------------------------------------------------------------------------------------------------
load("regionsAndHits.RData")

#--------------------------------------------------------------------------------------------------------------
assessTranscriptionFactorHits <- function(gene)
{
   timepoints <- seq_len(length(days))     # days and reps define in util.R
   x <- lapply(timepoints, function(i){
       total.hits <- extractHitData(regions, hits, days[i], reps[i], gene, "total")
       top.hits <- extractHitData(regions, hits, days[i], reps[i], gene, "top")
       total.density <- extractHitData(regions, hits, days[i], reps[i], gene, "total-density")
       top.density <- extractHitData(regions, hits, days[i], reps[i], gene, "top-density")
       tbl <- data.frame(total.hits=total.hits, top.hits=top.hits,
                         total.density=total.density, top.density=top.density)
       colnames(tbl) <- sprintf("%s.%s", gene, colnames(tbl))
       tbl
       })
   t(do.call(rbind, x))

} # assessTranscriptionFactorHits
#--------------------------------------------------------------------------------------------------------------
mtx <- reshapeAtacSeqData()[c("KLF1", "FLI1"),]
mtx <- rbind(mtx, assessTranscriptionFactorHits("KLF1"))
mtx <- rbind(mtx, assessTranscriptionFactorHits("FLI1"))
save(mtx, file="KLF1.FLI1.atac.srm.RData")



# 
# regions <- lapply(span, function(i) readAtacSeq(days[i], reps[i]))
# names(regions) <- sprintf("day.%d.%d", days, reps)
# regionCounts <- unlist(lapply(regions, nrow), use.names=FALSE)
#    # quick QC: do these open chromatin region counts correspond to file size?
# lineCounts <- unlist(lapply(span, function(i) lineCountAtacSeq(days[i], reps[i])))
# checkTrue(all(regionCounts == lineCounts))
# hits <- lapply(regions, function(tbl.region) findFimoHitsInAtacSeqRegions(tbl.region, 1e-4))
# names(hits) <- sprintf("day.%d.%d", days, reps)
# save(regions, hits, file="regionsAndHits.RData")
# 
# 
# 
# span <- seq_len(length(days))     # days and reps define in util.R
# 
#    total.klf1.hits.day.4.1 <- extractHitData(regions, hits, days[1], reps[1], gene, "total")
# 
