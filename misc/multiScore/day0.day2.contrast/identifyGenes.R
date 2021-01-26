library(jsonlite)
library(httr)
library(TrenaProjectErythropoiesis)
library(RUnit)

tpe <- TrenaProjectErythropoiesis();

expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")
checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
mtx <- getExpressionMatrix(tpe, expected[1])

cutoff.1 <- 0.25
cutoff.2 <- 2

goi.0 <- names(which(sapply(rownames(mtx), function(geneName) all(mtx[geneName, 1:2] < cutoff.1))))
print(length(goi.0))  # 902
goi.2 <- names(which(sapply(rownames(mtx), function(geneName) all(mtx[geneName, 3:4] > cutoff.2))))
print(length(goi.2))  # 10299

goi <- intersect(goi.0, goi.2)
length(goi)  # 11

goi.string <- toJSON(goi)
uri <- sprintf("http://localhost:8000/goEnrich")
body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=goi)))

r <- POST(uri, body=body.jsonString)

   #sprintf('{"geneSymbols": "%s"}', goi.string))
tbl <- fromJSON(content(r)[[1]])
dim(tbl)

require(RColorBrewer)
colors <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
lineType <- "b"
plot(mtx[goi[1],], col=sample(colors, size=1), type=lineType, ylim=c(0,9),
     main=sprintf("expression levels of %d genes unexpressed at day 0", length(goi)))
for(i in 2:length(goi)){
    lines(mtx[goi[i],], col=sample(colors, size=1), type=lineType)
    }

filename <- sprintf("day0.day2-%d-topGenes.RData", length(goi))
f <- file.path("~/github/TrenaProjectErythropoiesis/inst/extdata/geneSets/", filename)
save(goi, file=f)



