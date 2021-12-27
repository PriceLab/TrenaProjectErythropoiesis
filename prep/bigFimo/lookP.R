library(randomForest)
data.dir <- "~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints"
mtx.rna <- get(load(file.path(data.dir, "mtx-rna-pkfm-27170x14.RData")))
mtx.srm <- get(load(file.path(data.dir, "mtx-srm-copyNumber-100x13.RData")))
mtx.both <- get(load(file.path(data.dir, "mtx-rna-srm-2720x13.RData")))

early <- c("D0","D2","D4","D6","D7_5","D8")
late  <- c("D8_5", "D10","D10_5","D11","D11_5","D12","D14")



drs.p <- c("TCF3p", "CTCFp", "LDB1p", "STAT2p", "E2F4p", "STAT3p", "SUZ12p", "NR2C2p", 
           "TRIM33p", "DR1p", "TAL1p", "HLTFp", "WDHD1p", "SMC3p", "ZBTB7Ap", "TFDP1p",
           "STAT1p", "E2F8p", "GATA1p", "RAD21p", "KLF1p", "BCL11A")


#files <- system("find *p -name '*p.RData'", intern=TRUE)
files <- system("find *p -name '*p.lateTimepoints*p.RData'", intern=TRUE)
length(files)
files.good <- files[which(nchar(files) <= nchar("ZBTB7Ap/trena.p.lateTimepoints.ZBTB7Ap.RData"))]
length(files.good)
tbls <- list()
for(file in files){
    tbl <- get(load(file))
    name <- strsplit(file, "\\/")[[1]][1]
    tbls[[name]] <- tbl
    }

tbl <- do.call(rbind, tbls)
rownames(tbl) <- NULL
dim(tbl)
save(tbl, file="trena-91-ProteinModels-32603x10.RData")

tbl.rbp10 <- subset(tbl, rank < 11 & class=="rbp")   # 220 10
dim(tbl.rbp10)   # 38 10

head(as.data.frame(sort(table(tbl.rbp10$gene), decreasing=TRUE)), n=12)
#       Var1 Freq
# 1    DDX3X   36
# 2    RBM47   12
# 3    BUD13   11
# 4     FXR2    9
# 5    NCBP3    8
# 6   HNRNPC    7
# 7    DGCR8    6
# 8  IGF2BP2    5
# 9   GTF2F1    4
# 10   EIF3D    3
# 11 HNRNPA1    3

tbl.ddx3x.10 <- subset(tbl.rbp10, gene=="DDX3X")
new.order <- order(tbl.ddx3x.10$rank, decreasing=FALSE)
tbl.ddx3x.10 <- tbl.ddx3x.10[new.order,]

tbl.tf10 <- subset(tbl, rank < 11 & class=="tf")   # 220 10
dim(tbl.rbp10)   # 38 10

tbl.rbp.2 <- subset(tbl, gene %in% c("RBM47", "FXR2") & rank < 11)
tbl.rbp.2[order(tbl.rbp.2$gene),]

  #   gene betaLasso betaRidge spearmanCoeff pearsonCoeff     rfScore xgboost class rank  target
  #   FXR2     0.000  -108.790        -0.964       -0.897 26289622.85       0   rbp    7   CTCFp
  #   FXR2     0.000   -28.254        -0.821       -0.920   873123.61       0   rbp    5  GATA1p
  #   FXR2     0.000   -16.033        -0.857       -0.892  1205324.58       0   rbp    6   KLF1p
  #   FXR2     0.000   -53.384        -0.929       -0.859  6780227.04       0   rbp    6   LDB1p
  #   FXR2     0.000    -9.100        -0.929       -0.934   117906.76       0   rbp    5  NR2C2p
  #   FXR2     0.000   -23.469        -0.893       -0.901  1146141.22       0   rbp    7   SMC3p
  #   FXR2     0.000   -25.043        -0.893       -0.894  1789920.96       0   rbp    6   TAL1p
  #  RBM47     0.000    87.021         0.750        0.966   235128.10       0   rbp    8  RAD21p
  #  RBM47   174.941    15.928         0.857        0.981    14074.52       0   rbp    3  SUZ12p
  #  RBM47     0.000    81.576         0.643        0.938   670293.24       0   rbp    5 TRIM33p
  #  RBM47     0.000   319.811         0.643        0.948        0.00       0   rbp    7 ZBTB7Ap



as.data.frame(sort(table(tbl.tf10$gene), decreasing=TRUE))

    # 1  ONECUT2   12
    # 2     ETV2   10   Cooperative interaction of Etv2 and Gata2 regulates the development of endothelial and hematopoietic lineages
    # 3  BHLHE41    8
    # 4   MLXIPL    8
    # 5  CREB3L1    7
    # 6    FOXJ3    7
    # 7     TCF7    7
    # 8    VENTX    6
    # 9    RUNX3    5
    # 10   CREB1    4
    # 11    ELK4    4
    # 12  NFATC2    4
    # 13   SMAD2    4
    # 14   FOXC1    3
    # 15   OLIG1    3
    # 16    TAF1    3
    # 17   TGIF2    3
    # 18    E2F1    2
    # 19    EMX2    2
    # 20   FOXL1    2
    # 21   FOXP3    2
    # 22   GATA2    2
    # 23    IRF2    2
    # 24   KLF13    2
    # 25   MYBL2    2
    # 26   NR3C2    2
    # 27    NRF1    2



#----------------------------------------------------------------------------------------------------
explore.high.rfScore.ctcf.ddx3x <- function()
{
   load("CTCFp/trena.p.lateTimepoints.CTCFp.RData")
   tfs <- head(tbl.trena.p$gene, n=50)
   mtx.ctcf <- mtx.both[c("CTCFp", tfs),]
   round(fivenum(mtx.ctcf["CTCFp",]), digits=0)    # 33796  48298  67077  90862 113393 
   round(fivenum(mtx.ctcf["DDX3X",]), digits=0)    #  76    79    82   103   118 
   dim(mtx.ctcf)   # 51 13   
   mtx.ctcf["CTCFp", ] <- mtx.ctcf["CTCFp", ]/10000
   x <- randomForest(CTCFp ~ ., data=t(mtx.ctcf))
   head(x$importance)
   fivenum(mtx.ctcf["CTCFp",])
   fivenum(mtx.ctcf["TCF7",])
}
#----------------------------------------------------------------------------------------------------




