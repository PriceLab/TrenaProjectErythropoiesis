library(TrenaProjectErythropoiesis)
library(igvR)
library(later)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
library (RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
currentColorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
if(!exists("state"))
   state <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")){
   tp <- TrenaProjectErythropoiesis()
   setTargetGene(tp, "GATA2")
   tbl.enhancers <- getEnhancers(tp)
   }

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   later(function() showGenomicRegion(igv, "GATA2"), 3)
   later(function(){
            track <- DataFrameQuantitativeTrack("gh", tbl.enhancers[, c("chrom", "start", "end", "combinedScore")],
                                                autoscale=FALSE, color="blue", min=0, max=50)
            displayTrack(igv, track)
            }, 4)
   } # igv
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function(chromosome, start.loc, end.loc)
{
   directory <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
display_ATACseq_in_enhancers <- function()
{
   chrom <- "chr3"
   start.loc <- 128476488
   end.loc <- 128502543
   tbl.atac <- getATACseq(chrom, start.loc, end.loc)
   samples <- unique(tbl.atac$sample)
    # if the area above has enough span, then open chromatin is found in 11 samples:
    # "d04_rep1" "d04_rep2" "d08_rep1" "d10_rep1" "d10_rep2" "d11_rep1" "d11_rep2" "d12_rep1" "d12_rep2" "d16_rep1" "d16_rep2"
   atac.regions.by.sample <- list()
   for(sample.x in samples){
      tbl.sample <- subset(tbl.atac, sample==sample.x)[, c("chrom", "start", "end")]
      atac.regions.by.sample[[sample.x]] <- tbl.sample
      dim(tbl.sample)
      currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
      color <- colors[currentColorNumber]
      track <- DataFrameAnnotationTrack(sample.x, tbl.sample, color=color)
      displayTrack(igv, track)
      #write.table(tbl.sample, file=sprintf("tbl.%s.bed", sample.x), quote=FALSE, row.names=FALSE)
      } # for sample.x

   state$atac.regions.by.sample <- atac.regions.by.sample

} # display_ATACseq_inEnhancers
#------------------------------------------------------------------------------------------------------------------------
makeAtacBasedModels <- function()
{
  stopifnot("atac.regions.by.sample" %in% names(state))
  models <- list()
  samples <- names(state$atac.regions.by.sample)
  for(sample.x in samples[1:2]){
     printf("--- build model with sample %s", sample.x)
     tbl.regions <- state$atac.regions.by.sample[[sample.x]]
     models[[sample.x]] <- buildModel(tbl.regions)
     } # for sample.x

 state$models <- models

} # makeAtacBasedModels
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(tbl.regions)
{
   genome <- "hg38"
   tss <- getTranscriptsTable(tp)$tss[1]
   mtx.name <- "brandLabDifferentiationTimeCourse-27171x28"
   stopifnot(mtx.name %in% getExpressionMatrixNames(tp))
   mtx <- getExpressionMatrix(tp, mtx.name)
   dim(mtx)

   #db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs", "gata2.gh.fimoBindingSites.sqlite")
   db.name <- "fimoResults-10e-3-chr3-128383794-128647775.sqlite"
   stopifnot(file.exists(db.name))

   recipe <- list(title="fimo.atacseq",
                  type="fimo.database",
                  regions=tbl.regions,
                  geneSymbol=getTargetGene(tp),
                  tss=tss,
                  matrix=mtx,
                  db.host="localhost",
                  db.port=NA_integer_,
                  databases=list(db.name),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="fimoDatabase",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.4,
                  maxModelSize=100,
                  fimoPValueThreshold=1e-5,,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- FimoDatabaseModelBuilder(genome, getTargetGene(tp),  recipe, quiet=TRUE)
   x <- build(builder)

   x

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
displayMotifsForTF <- function(tf="TBX15")
{
  pfms <- query(MotifDb, c("hsapiens", tf), c("swissregulon", "jaspar2018", "hocomoco")) # 537
  mm <- MotifMatcher("hg38", as.list(pfms), quiet=TRUE)
  loc <- getGenomicRegion(igv)
  tbl.loc <- with(loc, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))

  tbl.matches <- findMatchesByChromosomalRegion(mm, tbl.loc, pwmMatchMinimumAsPercentage=90)
  dim(tbl.matches)
  trackName <- sprintf("%s-jaspar", tf)
  track <- DataFrameQuantitativeTrack("TBX15-jaspar", tbl.matches[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                      color="red", autoscale=FALSE, min=0, max=1, trackHeight=25)
  displayTrack(igv, track)

} # displayMotifsForTF
#------------------------------------------------------------------------------------------------------------------------
hopeRestored <- function()
{
   tbl.hope.restored <- data.frame(chrom=rep("chr3", 6),
                                   start=c(128475449, 128486955, 128483451, 128497527, 128483002, 128461602),
                                   end=c(128475456, 128486962, 128483458, 128497534, 128483009, 128461609),
                                   stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("hopeRestored?", tbl.hope.restored, color="magenta", trackHeight=30)
   displayTrack(igv, track)


} # hopeRestored
#------------------------------------------------------------------------------------------------------------------------
# buildModel <- function(tbl.regions)
# {
#
# } # buildModel
# see the correlation which strengthens the "TBX15 represses GATA2" hypothesis:
#   plot(mtx["GATA2",], ylim=c(0,12), col="blue")
#   points(mtx["TBX15",], col="darkgreen")
# reproduce_TBX15_represses_GATA2_hypothesis <- function()
# {
#    printf("--- reproduce_TBX15_represses_GATA2_hypothesis")
#
#    genome <- "hg38"
#    targetGene <- "GATA2"
#    setTargetGene(tp, targetGene)
#
#    chromosome <- "chr3"
#    tss <- 128493185
#    upstream <- 5000
#    downstream <- 5000
#       # strand-aware start and end: GATA2 is on the minus strand
#    start <- tss - downstream
#    end   <- tss + upstream
#
#    tbl.enhancers <- getEnhancers(tp)
#    gr.enhancers <- GRanges(tbl.enhancers)
#    atac.directory <- "../../../prep/import/atacPeaks"
#    tbl.atac <- read.table(file.path(atac.directory, "ATAC_Cord_d04_rep1_hg38_macs2_default_peaks.narrowPeak"),
#                           sep="\t", as.is=TRUE, nrow=-1)[, 1:3]
#    colnames(tbl.atac) <- c("chrom", "start", "end") #, "name", "wdith", "strand", "score1", "score2", "score3", "score4")
#
#    gr.atac <- GRanges(tbl.atac)
#    tbl.ov <- as.data.frame(findOverlaps(gr.atac, gr.enhancers))
#    tbl.ov
#    tbl.regions <- tbl.atac[unique(tbl.ov[,1]),]
#    rownames(tbl.regions) <- NULL
#
#    mtx.name <- "brandLabDifferentiationTimeCourse-filtered"
#    checkTrue(mtx.name %in% getExpressionMatrixNames(tp))
#    mtx <- getExpressionMatrix(tp, mtx.name)
#
#    db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs", "gata2.gh.fimoBindingSites.sqlite")
#    checkTrue(file.exists(db.name))
#
#    build.spec <- list(title="fimo.5000up.5000down",
#                       type="fimo.database",
#                       regions=tbl.regions,
#                       geneSymbol=targetGene,
#                       tss=tss,
#                       matrix=mtx,
#                       db.host="localhost",
#                       db.port=NA_integer_,
#                       databases=list(db.name),
#                       annotationDbFile=dbfile(org.Hs.eg.db),
#                       motifDiscovery="fimoDatabase",
#                       tfPool=allKnownTFs(),
#                       tfMapping="MotifDB",
#                       tfPrefilterCorrelation=0.4,
#                       maxModelSize=10,
#                       orderModelByColumn="rfScore",
#                       solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
#
#    builder <- FimoDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
#    x <- build(builder)
#
#    checkTrue(all(c("model", "regulatoryRegions") %in% names(x)))
#    tbl.model <- x$model
#    checkEquals(nrow(tbl.model), 10)
#
#
# } # reproduce_TBX15_represses_GATA2_hypothesis
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   display_ATACseq_in_enhancers()
   makeAtacBasedModels()

} # run
#------------------------------------------------------------------------------------------------------------------------
