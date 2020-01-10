tbl <- read.table("network.tsv", sep="\t", as.is=TRUE)
colnames(tbl) <- c("a", "b")
a <- unique(tbl$a)
b.raw <- strsplit(tbl$b, " ")
b <- unlist(b.raw)
length(b)
deleters.and <- grep("AND", b)
deleters.or <-   grep("OR", b)
length(deleters.and);
length(deleters.or)
deleters <- sort(c(deleters.and, deleters.or))
if(length(deleters) > 0)
   b <- b[-deleters]

