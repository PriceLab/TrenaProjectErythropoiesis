library(plyr)
library(randomForest)
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx.rna <- get(load(file.path(data.dir, "mtx-rna-pkfm-27170x14.RData")))
mtx.srm <- get(load(file.path(data.dir, "mtx-srm-copyNumber-100x13.RData")))
mtx.both <- get(load(file.path(data.dir, "mtx-rna-srm-2720x13.RData")))

early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")


files <- system("find . -name 'trena.*.RData'", intern=TRUE)
length(files)
excluders <- grep("imepoints", files)
files <- files[-excluders]
length(files)

tbls <- list()
for(file in files){
    load(file)
    name <- strsplit(file, "\\/")[[1]][2]
    tbl.trena.both$rank <- seq_len(nrow(tbl.trena.both))
    tbls[[name]] <- tbl.trena.both
    }

tbl <- do.call(rbind.fill, tbls)
rownames(tbl) <- NULL
dim(tbl)  # 17410 10

head(subset(tbl, class=="rbp" & rank < 3), n=20)


