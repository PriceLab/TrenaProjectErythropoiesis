load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-srm-copyNumber-100x13.RData")
genes.todo <- intersect(rownames(mtx.srm), list.files(include.dirs=TRUE))
for(gene in genes.todo){
    fimo.complete <- file.exists(sprintf("%s/tbl.fimo.%s.RData", gene, gene))
    tms.complete  <- file.exists(sprintf("%s/trena.%s.RData", gene, gene))
    if(fimo.complete & !tms.complete)
        printf("cd %s; Rscript ../runTMS.R %s; ls -l trena.%s.RData; cd ..", gene, gene, gene)
    }
