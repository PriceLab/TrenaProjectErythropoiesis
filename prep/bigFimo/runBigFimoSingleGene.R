args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1)
targetGene <- args[1]
if(file.exists(targetGene)) {
    message(sprintf("%s directory already exists, exiting", targetGene))
    quit("no")
    }            


source("~/github/bigFimo/R/BigFimo.R")

files <- list.files(path=targetGene, pattern="*.RData")

if(length(files) > 0)
   unlink(file.path(targetGene, files))

processCount <- 30
fimoThreshold <- 1e-3
gh.elite.only <- FALSE
maxGap.between.oc.and.gh <- 5000
chrom <- NA
start <- NA
end   <- NA

bf <-  BigFimo$new(targetGene,
                   project="BrandLabErythropoiesis",
                   processCount,
                   fimoThreshold,
                   gh.elite.only,
                   maxGap.between.oc.and.gh,
                   chrom=chrom, start=start, end=end)
tbl.gh <- bf$get.tbl.gh()

bf$calculateRegionsForFimo()
tbl.gh.oc <- bf$get.tbl.gh.oc()

filenames.roi <- bf$createFimoTables()
bf$runMany()


Sys.sleep(10)
completed <- FALSE
actual.processes.needed <- length(bf$getFimoRegionsFileList())

while(!completed){
    file.count <- length(list.files(path=targetGene, pattern="^fimo.*"))
    completed <- (file.count == actual.processes.needed)
    if(!completed){
        printf("waiting for completion: %d/%d", file.count, actual.processes.needed)
        Sys.sleep(3)
    }
} # while

printf("complete %d/%d", actual.processes.needed, actual.processes.needed)

result.files <- list.files(path=targetGene, pattern="^fimo.*")
tbls <- list()
for(file in result.files){
    tbl <- get(load(file.path(targetGene, file)))
    tbls[[file]] <- tbl
    }
tbl.fimo <- do.call(rbind, tbls)
tbl.fimo <- tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),]
rownames(tbl.fimo) <- NULL

