library(plyr)
files <- sort(system("find . -name 'trena*RData'", intern=TRUE))
tbls <- list()
i <- 0

files.test <- grep("BACH1", files, v=TRUE)

for(file in files){
  i <- i + 1
  printf("file: %s", file)
  load(file)
  tbl <- tbl.trena.both
  new.order <- order(abs(tbl$spearmanCoeff), decreasing=TRUE)
  tbl <- tbl[new.order,]
  if("rank" %in% colnames(tbl)){
      rank.col <- grep("rank", colnames(tbl))
      tbl <- tbl[, -rank.col]
      }
  tbl$spear.pos <- order(tbl$spearmanCoeff, decreasing=TRUE)
  tbl$spear.abs  <- order(abs(tbl$spearmanCoeff), decreasing=TRUE)
  tbl$pear.pos  <- order(tbl$pearsonCoeff, decreasing=TRUE)
  tbl$pear.abs  <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
  tbl$rf <- order(tbl$rfScore, decreasing=TRUE)
  tbl$lasso.pos <- order(tbl$betaLasso, decreasing=TRUE)
  tbl$lasso.abs <- order(abs(tbl$betaLasso), decreasing=TRUE)
  tbl$sourceFile <- file
  # if(tbl$target[1]=="OGT") browser()
  tbls[[i]] <- tbl
  }

tbl <- do.call(rbind.fill, tbls)
dim(tbl)
save(tbl, file=sprintf("tbl.trena-%d-targets.RData", length(files)))
