library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
tpe <- TrenaProjectErythropoiesis()
print(load("~/github/TrenaProjectErythropoiesis/inst/extdata/geneSets/day0.day2-11-topGenes.RData"))
goi

tbls <- list()

for(gene in goi){
   tms <- TrenaMultiScore(tpe, gene)
   getGeneHancerRegion(tms)
   findOpenChromatin(tms)
   tbl.oc <- getOpenChromatin(tms)
   if(nrow(tbl.oc) == 0){
       printf("bailing out of %s", gene)
       tbls[[gene]] <- getMultiScoreTable(tms)
   }else{
       findFimoTFBS(tms, fimo.threshold=1e-3)
       scoreMotifHitsForConservation(tms)
       scoreMotifHitsForGeneHancer(tms)
       addDistanceToTSS(tms)
       mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
       addGeneExpressionCorrelations(tms, mtx)
       addGenicAnnotations(tms)
       addChIP(tms)
       tbl <- getMultiScoreTable(tms)
       dim(tbl)
       tbl$cor[which(is.na(tbl$cor))] <- 0
       tbl$motifScore <- round(-log10(tbl$p.value), 2)
       dim(tbl)
       tbls[[gene]] <- tbl
       } # else: some open chromatin found
   } # for gene

