library(GenomicRanges)
tbl.041 <- get(load("ATAC_brand_d04_rep1.RData"))[, c("chrom", "start", "end")]
tbl.042 <- get(load("ATAC_brand_d04_rep2.RData"))[, c("chrom", "start", "end")]
tbl.08  <- get(load("ATAC_brand_d08_rep1.RData"))[, c("chrom", "start", "end")]
tbl.101 <- get(load("ATAC_brand_d10_rep1.RData"))[, c("chrom", "start", "end")]
tbl.102 <- get(load("ATAC_brand_d10_rep2.RData"))[, c("chrom", "start", "end")]
tbl.111 <- get(load("ATAC_brand_d11_rep1.RData"))[, c("chrom", "start", "end")]
tbl.112 <- get(load("ATAC_brand_d11_rep2.RData"))[, c("chrom", "start", "end")]
tbl.121 <- get(load("ATAC_brand_d12_rep1.RData"))[, c("chrom", "start", "end")]
tbl.122 <- get(load("ATAC_brand_d12_rep2.RData"))[, c("chrom", "start", "end")]
tbl.161 <- get(load("ATAC_brand_d16_rep1.RData"))[, c("chrom", "start", "end")]
tbl.162 <- get(load("ATAC_brand_d16_rep2.RData"))[, c("chrom", "start", "end")]


gr.041 <- GRanges(tbl.041)
gr.042 <- GRanges(tbl.042)
gr.08 <- GRanges(tbl.08)
gr.101 <- GRanges(tbl.101)
gr.102 <- GRanges(tbl.102)
gr.111 <- GRanges(tbl.111)
gr.112 <- GRanges(tbl.112)
gr.121 <- GRanges(tbl.121)
gr.122 <- GRanges(tbl.122)
gr.161 <- GRanges(tbl.161)
gr.162 <- GRanges(tbl.162)

gr.list <- c(gr.041, gr.042, gr.08 , gr.101, gr.102, gr.111, gr.112, gr.121, gr.122, gr.161, gr.162)

x <- unlist(as(gr.list, "GRangesList"))
length(x)  # 1,391,318

head(x)
tbl.atacMerged <- as.data.frame(x)
head(tbl.atacMerged)
save(tbl.atacMerged, file="tbl.atacMerged.RData")
