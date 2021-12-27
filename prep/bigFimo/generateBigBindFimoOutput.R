load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-srm-copyNumber-100x13.RData")
genes.todo <- intersect(rownames(mtx.srm), list.files(include.dirs=TRUE))
for(gene in genes.todo){
    if(!file.exists(sprintf("%s/tbl.fimo.%s.RData", gene, gene)))
        printf("cd %s; Rscript ~/github/bigFimo/R/rbind.R %s; cd ..", gene, gene)
    }
