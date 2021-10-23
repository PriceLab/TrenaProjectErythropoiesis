library(BSgenome.Hsapiens.UCSC.hg38)
library(RPostgreSQL)
library(memes)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
tbl.models <- get(load("trena-92-ProteinModels-32603x10.RData"))

#------------------------------------------------------------------------------------------------------------------------
# we know from earlier work that KLF1 has 7 DDX3X eCLIP binding sites in K562 cells.   search them for motifs.
quick.demo <- function()
{

   targetGenes <- sub("p$", "", unique(tbl.models$target))
   tbls <- list()
   rbp <- "DDX3X"

   for(targetGene in targetGenes){
      printf("--- %s", targetGene)
      query <- sprintf("select * from rbp where gene='%s' and target='%s'", rbp, targetGene)
      tbl.out <- dbGetQuery(geneRegDB, query)
      tbl.out <- subset(tbl.out, celltype=="K562")
      tbls[[targetGene]] <- tbl.out
      }
   length(tbls)
   tbl.out <- do.call(rbind, tbls)
   tbl.out <- subset(tbl.out, width > 6)
   dim(tbl.out)
   dna.string.set <- with(tbl.out, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, endpos))
   names(dna.string.set) <- with(tbl.out, sprintf("%s-%s:%d-%d-%s-%s", target, chrom, start, endpos, gene, celltype))
   length(dna.string.set) # [1] 261
   fasta.filename <- "targetGenes-97-ddx3x-k562.fa"
   writeXStringSet(dna.string.set, fasta.filename)

   cmd <- sprintf("~/meme/bin/meme %s -dna -oc meme.out -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun  classic -revcomp -markov_order 0", fasta.filename)

   system(cmd)
   system("open meme.out/meme.html")

} # quick.demo
#------------------------------------------------------------------------------------------------------------------------
# abandonded for now: library(memes)
use.bioc.memeSuite <- function()
{
   fasta.filename <- "targetGenes-97-ddx3x-k562.fa"
   dreme_out <- runDreme(fasta.filename, "shuffle", evalue = 1, outdir = tempdir())

} # use.bioc.memeSuite
#------------------------------------------------------------------------------------------------------------------------
runFimoOnDiscoveredMotifs <- function()
{
  source("~/github/fimoService/batchMode/fimoBatchTools.R")
  data.dir <- "~/github/TrenaProjectErythropoiesis/explore/rbp/discordants"
  meme.file <- "meme.out/motif-21.meme"
  meme.file.path <- file.path(data.dir, meme.file)
     #HLTF has strong corDelta (-1.56) and strong motif hits.  can we find them?
  file.exists(meme.file.path)
  tbls <- list()
  targets <- unique(tbl.out$target)
  for(TARGET in targets){
     tbl.target <- subset(tbl.out, target==TARGET)[, c("chrom", "start", "endpos")]
     rownames(tbl.target) <- NULL
     colnames(tbl.target) <- c("chrom", "start", "end")
     tbl.target$sequence_name <- with(tbl.target, sprintf("%s:%d-%d", chrom, start, end))
     tbl.fimo <- fimoBatch(tbl.target, matchThreshold=1e-3, genomeName="hg38",
                           pwmFile=meme.file.path, expandResultsWithMotifDb=FALSE)
     tbls[[TARGET]] <- tbl.fimo
     }


  tbl.all <- tbl.out
  rownames(tbl.all) <- NULL
  colnames(tbl.all) <- c("chrom", "start", "end")
  tbl.all$sequence_name <- with(tbl.all, sprintf("%s:%d-%d", chrom, start, end))

  tbl.hltf <- subset(tbl.out, target=="HLTF")[, c("chrom", "start", "endpos")]
  rownames(tbl.hltf) <- NULL
  colnames(tbl.hltf) <- c("chrom", "start", "end")
  tbl.hltf$sequence_name <- with(tbl.hltf, sprintf("%s:%d-%d", chrom, start, end))
  alt.meme.file <- "~/github/TrenaProjectErythropoiesis/misc/mreg/topTen.meme"
  tbl.fimo <- fimoBatch(tbl.all, matchThreshold=1e-3, genomeName="hg38",
                        pwmFile=meme.file.path, expandResultsWithMotifDb=FALSE)
#                        pwmFile=alt.meme.file)
  dim(tbl.fimo)
  }
#----------------------------------------------------------------------------------------------------
compareCorAndScore <- function()
{
   tbl.motifs <- get(load("tbl.motifs3.regions52.RData"))
   tbl.cor <- get(load("tbl.srm.rna.cor.early.late.all.delta.100genes.RData"))

   high.cor.late.genes <- subset(tbl.cor, cor.late < -0.5)$gene
   length(high.cor.late.genes) # 28

   high.cor.delta.genes <- subset(tbl.cor, corDelta < -0.7)$gene
   length(high.cor.delta.genes) # 26

   length(intersect(high.cor.late.genes, high.cor.delta.genes)) # 20

   dim(subset(tbl.cor, cor.late > 0.5))
   genes.agree <- subset(tbl.cor, cor.late > 0.5)$gene
   length(genes.agree)  # 38

   length(intersect(high.cor.late.genes, tbl.motifs$gene))/length(high.cor.late.genes)  # 36%
   length(intersect(genes.agree, tbl.motifs$gene))/length(genes.agree)                  # 18%

   mean(subset(tbl.motifs, gene %in% high.cor.late.genes)$score)                        # 11.67
   mean(subset(tbl.motifs, gene %in% genes.agree)$score)                                #  9.9

   ddx3x.regulated.genes <- unique(sub("p$", "", subset(tbl.models, gene=="DDX3X" & rank < 10)$target))
   length(ddx3x.regulated.genes)  # 35
   length(high.cor.late.genes)       # 28
   top.hits <- intersect(high.cor.late.genes, ddx3x.regulated.genes) # 15
   length(top.hits) # 15
   length(intersect(top.hits, tbl.motifs$gene))/length(top.hits)  # 40%


} # compareCorAndScore
#----------------------------------------------------------------------------------------------------
