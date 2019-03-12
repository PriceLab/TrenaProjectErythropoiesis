# atac-seq, srm, jeff & marjorie: only 5 shared time points: days 4, 8*, 10, 11, 12, FLI1 and KLF1
# call FLI1 and KLF1 hits in all ATAC-seq in the 5 shared time points, 4 of which have replicates
#--------------------------------------------------------------------------------------------------------------
if(!exists("readAtacSeq"))
    source("util.R")
#--------------------------------------------------------------------------------------------------------------
days <- c(4, 4, 8, 10, 10, 11, 11, 12, 12)
reps <- c(1, 2, 1, 1,   1,  1,  2,  1,  2)
span <- seq_len(length(days))
regions <- lapply(span, function(i) readAtacSeq(days[i], reps[i]))
names(regions) <- sprintf("day.%d.%d", days, reps)
regionCounts <- unlist(lapply(regions, nrow), use.names=FALSE)
   # quick QC: do these open chromatin region counts correspond to file size?
lineCounts <- unlist(lapply(span, function(i) lineCountAtacSeq(days[i], reps[i])))
checkTrue(all(regionCounts == lineCounts))
hits <- lapply(regions, function(tbl.region) findFimoHitsInAtacSeqRegions(tbl.region, 1e-4))
names(hits) <- sprintf("day.%d.%d", days, reps)
save(regions, hits, file="regionsAndHits.RData")

