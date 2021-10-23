library(RUnit)
library(RPostgreSQL)
library(BSgenome.Hsapiens.UCSC.hg38)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx <- get(load(file.path(data.dir, "mtx-rna-srm-27270x13.RData")))
dim(mtx)
early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")
#----------------------------------------------------------------------------------------------------
calculate.discordance <- function()
{
   tbls <- list()
   proteins <- grep("p$", rownames(mtx), v=TRUE)
   checkEquals(length(proteins), 100)

   for(protein in proteins){
       gene <- sub("p$", "", protein)

       gene.late <- mtx[gene, late]
       prot.late <- mtx[protein, late]

       gene.early <- mtx[gene, early]
       prot.early<- mtx[protein, early]

       gene.all <- mtx[gene,]
       prot.all <- mtx[protein,]

       cor.early <- cor(gene.early, prot.early, method="spearman", use="pairwise.complete")
       cor.late <- cor(gene.late, prot.late, method="spearman", use="pairwise.complete")
       cor.all <- cor(gene.all, prot.all, method="spearman", use="pairwise.complete")

       tbl <- data.frame(gene=gene, protein=protein,
                         cor.early = cor.early,
                         cor.late = cor.late,
                         cor.all = cor.all,
                         stringsAsFactors=FALSE)
       tbls[[protein]] <- tbl
       } # for protein

   tbl <- do.call(rbind, tbls)
   rownames(tbl) <- NULL
   dim(tbl)
   tbl$corDelta <- tbl$cor.late - tbl$cor.early
       # round off the values for greater scrutability, no loss of relevant information
   tbl$cor.early <- round(tbl$cor.early, digits=2)
   tbl$cor.late <- round(tbl$cor.late, digits=2)
   tbl$cor.all <- round(tbl$cor.all, digits=2)
   tbl$corDelta <- round(tbl$corDelta, digits=2)
   new.order <- order(tbl$corDelta, decreasing=FALSE)
   new.order <- order(tbl$corDelta, decreasing=FALSE)
   tbl <- tbl[new.order,]
   rownames(tbl) <- NULL
   tbl.cor <- tbl
   save(tbl.cor, file="tbl.srm.rna.cor.early.late.all.delta.100genes.RData")
   dim(tbl.cor)  # 100 6
   tbl.cor

} # calculate.discordance
#----------------------------------------------------------------------------------------------------
get.ddx3x.binding.sites <- function(tbl.discordance)
{
   save.name <- "tbl.ddx3x.hits.100.genes.RData"
   if(file.exists(save.name)){
      tbl.hits <- get(load(save.name))
      return(tbl.hits)
      }

   targetGenes <- tbl.discordance$gene
   rbp <- "DDX3X"
   tbls <- list()

   i <- 0
   for(gene in targetGenes){
      query <- sprintf("select * from rbp where gene='%s' and target='%s' and celltype='K562'",
                       rbp, gene)
      i <- i + 1
      printf("   %2d) %s", i, gene)
      tbls[[gene]] <- dbGetQuery(geneRegDB, query)
      }

   tbl.hits <- do.call(rbind, tbls)
   rownames(tbl.hits) <- NULL
   save(tbl.hits, file=save.name)

   return(tbl.hits)

} # get.ddx3x.binding.sites
#----------------------------------------------------------------------------------------------------
build.discordanceAndHitsTable <- function()
{
   tbl.discordance <- calculate.discordance()
   dim(tbl.discordance)
   head(tbl.discordance)
   tbl.hits <- get.ddx3x.binding.sites(tbl.discordances)
   tbl.hitCount <- data.frame(gene=tbl.discordance$gene, count=0, stringsAsFactors=FALSE)
   list.hitcount <- rep(0, nrow(tbl.discordance))
   names(list.hitcount) <- tbl.discordance$gene
   list.hitcount
   list.hits <- as.list(sort(table(tbl.hits$target), decreasing=TRUE))
   list.hitcount[names(list.hits)] <- as.integer(list.hits)
   tbl.discordance$ddx3x.hits <- as.integer(list.hitcount)

      # add a column 0/1, indicating the presence or absence of any ddx3x hit

   tbl.discordance$hit.01 <- 0
   tbl.discordance$hit.01[tbl.discordance$ddx3x.hits > 0] <- 1
   table(tbl.discordance$hit.01)   #  0  1
                                   # 26 74

} # build.discordanceAndHitsTable
#----------------------------------------------------------------------------------------------------
compare.ddx3x.rank.to.discordance <- function()
{
   tbl.trena <- get(load("trena-92-ProteinModels-32603x10.RData"))
   tbl.ddx3x <- subset(tbl.trena, gene=="DDX3X")[, c("target", "rank")]
   dim(tbl.ddx3x)
   ddx3x.rank <- rep(NA, nrow(tbl.discordance))
   names(ddx3x.rank) <-  tbl.discordance$protein
   ddx3x.rank[tbl.ddx3x$target] <- tbl.ddx3x$rank
   tbl.discordance

} # compare.ddx3x.rank.to.discordance
#----------------------------------------------------------------------------------------------------
linear.models <- function()
{
   tbl.x <- subset(tbl.trena, target=="KLF1p")[1:10,]
   summary(lm(KLF1p ~ 0 + TCF7 + DDX3X +  BUD13 +  FXR2, data=as.data.frame(t(mtx))))

} # linear.models
#----------------------------------------------------------------------------------------------------
# using this function (21 oct 2021) to introduce the ddx3 slides
# ~/github/TrenaProjectErythropoiesis/explore/slides/ddx3x.pptx
klf1.as.example <- function()
{
   tbl.discordance <- calculate.discordance()
   head(tbl.discordance)

   tbl.trena <- get(load("trena-92-ProteinModels-32603x10.RData"))
   tbl.klf1 <- subset(tbl.trena, target=="KLF1p")
   dim(tbl.klf1)   # 353 10
   new.order <- order(abs(tbl.klf1$spearmanCoeff), decreasing=TRUE)
   tbl.klf1 <- tbl.klf1[new.order,]
   head(tbl.klf1)
   tbl.bach1p <- subset(tbl.trena, target=="BACH1p")
   dim(tbl.bach1p)   #330 10
   new.order <- order(abs(tbl.bach1p$spearmanCoeff), decreasing=TRUE)
   tbl.bach1p <- tbl.bach1p[new.order,]
   head(tbl.bach1p, n=10)

   ddx3x.hot.targets <- sub("p$", "", subset(tbl.trena, gene=="DDX3X" & rank <= 3)$target)
   subset(tbl.discordance, gene %in% ddx3x.hot.targets & corDelta < -1)

} # klf1.as.example
#----------------------------------------------------------------------------------------------------
