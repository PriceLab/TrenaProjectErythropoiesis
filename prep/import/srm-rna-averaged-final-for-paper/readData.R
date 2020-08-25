mtx.rna <- as.matrix(read.table("rna.tsv", sep="\t", as.is=TRUE, header=TRUE, row.names=1, nrow=-1))
mtx.srm <- as.matrix(read.table("srm.tsv", sep="\t", as.is=TRUE, header=TRUE, row.names=1, nrow=-1))
dim(mtx.rna)
dim(mtx.srm)
preferred.colnames <-  c("d_0", "d_2", "d_4", "d_6", "d_7.5", "d_8", "d_8.5", "d_10", "d_10.5", "d_11", "d_11.5", "d_12", "d_14")
colnames(mtx.rna) <- preferred.colnames
colnames(mtx.srm) <- preferred.colnames

setdiff(rownames(mtx.rna), rownames(mtx.srm))
setdiff(rownames(mtx.srm), rownames(mtx.rna))

common.names <- sort(intersect(rownames(mtx.rna), rownames(mtx.srm)))
length(common.names)

mtx.rna <- mtx.rna[common.names,]
mtx.srm <- mtx.srm[common.names,]
save(mtx.rna, mtx.srm, file="srm.rna.averaged.clean.RData")
fivenum(mtx.rna)
fivenum(mtx.srm)
boxplot(asinh(as.numeric(mtx.rna)), asinh(as.numeric(mtx.srm)))



