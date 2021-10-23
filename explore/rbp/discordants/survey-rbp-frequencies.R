library(RPostgreSQL)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                       dbname="genereg2021", host="khaleesi")
query <- sprintf("select * from rbp limit 3")
dbGetQuery(geneRegDB, query)

rbps <- c("DDX3X",         # 1684009
          "FXR2",          #   69694
          "RBM47",         #   45759
          "BUD13",         #   96181
          "DGCR8",         #   86482
          "DICER1",        #   18449
          "GPKOW",         #   13999
          "HNRNPA1",       # 1275055
          "HNRNPC",        # 7995744
          "IGF2BP1",       #  164484
          "IGF2BP3",       #  155409
          "LIN28B",        #  564575
          "NONO",          #   48197
          "NPM1",          #     897
          "PUM2",          #   44530
          "SF3A3",         #   33000
          "TROVE2")        #    6945

tbl.counts <- data.frame(
     rbp=c("DDX3X", "FXR2", "RBM47", "BUD13", "DGCR8", "DICER1", "GPKOW", "HNRNPA1",
           "HNRNPC", "IGF2BP1", "IGF2BP3", "LIN28B", "NONO", "NPM1", "PUM2",
           "SF3A3", "TROVE2"),
    model.count=c(11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    stringsAsFactors=FALSE)


freqs <- list()
for(rbp in rbps){
   query <- sprintf("select count(*) from rbp where gene='%s'", rbp)
   freqs[[rbp]] <- dbGetQuery(geneRegDB, query)
   printf("%20s: %8d", rbp, freqs[[rbp]][1,1])
   }
genome.count <- unlist(lapply(freqs, function(freq) return(freq[1,1])), use.names=FALSE)

tbl.counts$genome <- genome.count
tbl.counts$x <- with(tbl.counts, model.count/genome*1000000)
tbl.counts$x <- round(tbl.counts$x, digits=2)
colnames(tbl.counts)[4] <- "per million total binding sites"

