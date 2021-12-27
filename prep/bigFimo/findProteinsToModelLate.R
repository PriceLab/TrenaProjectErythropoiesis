all.gene.dirs <- system('find . -name "trena.*.RData" | egrep -v imepoints', intern=TRUE)
length(all.gene.dirs)     # 100
head(all.gene.dirs)

tokens <- strsplit(all.gene.dirs, "/", fixed=TRUE)
all.genes <- sort(unlist(lapply(tokens, "[", 2)))
head(all.genes)

already.complete.dirs <- system('find . -name "trena.p.*.RData"', intern=TRUE)
length(already.complete.dirs)  # 23

tokens <- strsplit(already.complete.dirs, "/", fixed=TRUE)
complete.genes <- sort(unlist(lapply(tokens, "[", 2)))

complete.genes <- sub("p$", "", complete.genes)
head(complete.genes)
genes.to.do <- setdiff(all.genes, complete.genes)
length(genes.to.do) # 77 
