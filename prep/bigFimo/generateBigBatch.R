load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-srm-copyNumber-100x13.RData")
genes.todo <- setdiff(rownames(mtx.srm), list.files(include.dirs=TRUE))
minutes <- 40
interval <- 30
for(gene in genes.todo){
    minutes <- minutes + interval
    printf("echo /usr/bin/Rscript runBigFimoSingleGene.R %s | at now + %d minute", gene, minutes)
    }
