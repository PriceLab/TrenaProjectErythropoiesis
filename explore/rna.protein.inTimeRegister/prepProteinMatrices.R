srm.file <- "~/github/TrenaProjectErythropoiesis/prep/import/srm-rna-averaged-final-for-paper/srm.tsv"
tbl.srm <- read.table(srm.file, sep="\t", header=TRUE, nrow=-1)
rownames(tbl.srm) <- tbl.srm$Protein.Name

mtx.srm <- as.matrix(tbl.srm[, -1])
save(mtx.srm, file="mtx.srm.RData")

itraq.file <- "~/github/TrenaProjectErythropoiesis/prep/iTRAQ/iTRAQ-erythro.RData"
tbl.itraq <- get(load(itraq.file))
dim(tbl.itraq)  # 3831 12
tbl.itraq[1:10, 1:10]


no.name <- grep("No Name", tbl.itraq$Gene.Name)
tbl.itraq <- tbl.itraq[-no.name,]
uniprot.dups <- which(duplicated(tbl.itraq$Uniprot.Accession.Number..AC.))
length(uniprot.dups)
tbl.itraq <- tbl.itraq[-uniprot.dups,]

dups <- which(duplicated(tbl.itraq$Gene.Name))
length(dups) # 19
dup.genes <- tbl.itraq$Gene.Name[dups]
tbl.dups <- subset(tbl.itraq, Gene.Name %in% dup.genes)

tbl.keep <- subset(tbl.itraq, !Gene.Name %in% dup.genes)
nrow(tbl.itraq)
nrow(tbl.keep)
nrow(tbl.dups)

tbl.restore <- tbl.dups[grep("^[PQO]", tbl.dups[, 1]),]
length(which(duplicated(tbl.keep$Gene.Name))) # 0
dups2 <- which(duplicated(tbl.restore$Gene.Name)) # 1: TMPO, just drop the second one
tbl.restore <- tbl.restore[-dups2,]
length(which(duplicated(tbl.restore$Gene.Name))) # 0

length(which(duplicated(rbind(tbl.keep, tbl.restore)$Gene.Name))) # 0
tbl.itraq <- rbind(tbl.keep, tbl.restore)
tbl.itraq <- tbl.itraq[order(tbl.itraq$Gene.Name),]
rownames(tbl.itraq) <- tbl.itraq$Gene.Name
mtx.itraq <- as.matrix(tbl.itraq[, -(1:3)])
save(mtx.itraq, file="mtx.itraq.RData")


#--------------------------------------------------------------------------------
# now ensure that the columns of these two matrices are expanded to match
# the brand lab expression matrices
#--------------------------------------------------------------------------------
f <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"
f <- "~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-27171x28.RData"
mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/expression/brandLabDifferentiationTimeCourse-27171x28-namesCorrected.RData"))
dim(mtx.rna)
wdth(80)
colnames(mtx.rna)
#  [1] "day0.r1"  "day0.r2"  "day2.r1"  "day2.r2"  "day4.r1"  "day4.r2"
#  [7] "day6.r1"  "day6.r2"  "day7.r1"  "day7.r2"  "day8.r1"  "day8.r2"
# [13] "day8.r1"  "day8.r2"  "day10.r1" "day10.r2" "day10.r1" "day10.r2"
# [19] "day11.r1" "day11.r2" "day11.r1" "day11.r2" "day12.r1" "day12.r2"
# [25] "day14.r1" "day14.r2" "day16.r1" "day16.r2"

#                          rna             itraq       srm
preferred.colnames <- c("day.00.r1",   #    day.0       day.0
                        "day.00.r2",   #    day.0       day.0
                        "day.02.r1",   #    day.2       day.2
                        "day.02.r2",   #    day.2       day.2
                        "day.04.r1",   #    day.4       day.4
                        "day.04.r2",   #    day.4       day.4
                        "day.06.r1",   #    day.6       day.6
                        "day.06.r2",   #    day.6       day.6
                        "day.07_5.r1", #     -          day.75
                        "day.07_5.r2", #     -          day.75
                        "day.08.r1",   #    day.8       day.8
                        "day.08.r2",   #    day.8       day.8
                        "day.08_5.r1", #     -          day.85
                        "day.08_5.r2", #     -          day.85
                        "day.10.r1",   #    day.10      day.10
                        "day.10.r2",   #    day.10      day.10
                        "day.10_5.r1", #      -         day.105
                        "day.10_5.r2", #      -         day.105
                        "day.11.r1",   #      -         day.11
                        "day.11.r2",   #      -         day.11
                        "day.11.5.r1", #      -         day.115
                        "day.11.5.r2", #      -         day.115
                        "day.12.r1",   #    day.12      day.12
                        "day.12.r2",   #    day.12      day.12
                        "day.14.r1",   #    day.14      day.14
                        "day.14.r2")   #    day.14      day.14
                        #"day.16.r1",  #     -            -
                        #"day.16.r2")  #     -            -

c.drop <- c(9,10,13,14,17:22,27,28)
mtx.rna.conformant <- mtx.rna[, -c.drop]
dim(mtx.rna.conformant)   # 27171 16
save(mtx.rna.conformant,
     file="brandLabDifferentiationTimeCourse-27171x28-namesCorrected-conformant.RData")

mtx.itraq <- get(load("~/github/CHD1/tms/mtx.itraq.RData"))
dim(mtx.itraq)
colnames(mtx.itraq)
mtx.itraq <- mtx.itraq[, -1] # get rid of percent.coverage
dim(mtx.itraq)
colnames(mtx.itraq)

coi <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
mtx.itraq.conformant <- mtx.itraq[, coi]
colnames(mtx.itraq.conformant) <- colnames(mtx.rna.conformant)
save(mtx.itraq.conformant, file="mtx.itraq-conformant.RData")

mtx.srm <- get(load("~/github/CHD1/tms/mtx.srm.RData"))
dim(mtx.srm)  # 105 13
colnames(mtx.srm)
coi <- c(1,1,2,2,3,3,4,4,6,6,8,8,12,12,13,13)
mtx.srm.conformant <- mtx.srm[, coi]
colnames(mtx.srm.conformant) <- colnames(mtx.rna.conformant)
save(mtx.srm.conformant, file="mtx.srm-conformant.RData")

dim(mtx.rna.conformant)     # 27171 16
dim(mtx.srm.conformant)     #   105 16
dim(mtx.itraq.conformant)   #  3798 16

colnames(mtx.rna.conformant)
colnames(mtx.srm.conformant)
colnames(mtx.itraq.conformant)



