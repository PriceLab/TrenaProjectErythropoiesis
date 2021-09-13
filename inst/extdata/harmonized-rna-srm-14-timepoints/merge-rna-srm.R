f <- "mtx-rna-pkfm-27170x14.RData"
mtx.rna <- get(load(f))
dim(mtx.rna)
f <- "mtx-srm-copyNumber-100x13.RData"
mtx.srm <- get(load(f))
dim(mtx.srm)

coi <- intersect(colnames(mtx.rna), colnames(mtx.srm))
length(coi)
head(rownames(mtx.srm))
rownames(mtx.srm) <- sprintf("%sp", rownames(mtx.srm))
head(rownames(mtx.srm))
mtx.both <- rbind(mtx.rna[, coi], mtx.srm[, coi])
grep("BACH", rownames(mtx.both), v=TRUE)

mtx.rna.srm <- mtx.both
dim(mtx.rna.srm)
save(mtx.rna.srm, file="mtx-rna-srm-2720x13.RData")

