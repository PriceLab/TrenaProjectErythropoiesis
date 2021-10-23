onecut2.repression <- function()
{

   fivenum(tbl$corDelta[tbl$corDelta < 0]) # -1.69 -1.11 -0.73 -0.53 -0.01

   divergent.late.strong1 <- sort(unique(subset(tbl, corDelta <= -1.11)$protein))
   divergent.late.strong2 <- sort(unique(subset(tbl, corDelta > -1.11 & corDelta <= -0.71 )$protein))
   divergent.late.strong3 <- sort(unique(subset(tbl, corDelta > -0.71 & corDelta <= -0.51 )$protein))
   divergent.late.weak <- sort(unique(subset(tbl, corDelta > -0.51 & corDelta < 0)$protein))

   proteins.of.interest <- c(divergent.late.strong1, divergent.late.strong2, divergent.late.strong3)
   tbl.onecut2 <- subset(tbl.models, gene=="ONECUT2" & class=="tf" & rank <= 5 & target %in% proteins.of.interest)
   rownames(tbl.onecut2) <- NULL

#           gene   betaLasso  betaRidge spearmanCoeff pearsonCoeff      rfScore xgboost class rank  target
#  2501  ONECUT2       0.000 -80856.558        -0.893       -0.915 17188310.146       0    tf    4   CTCFp
#  3607  ONECUT2   -6098.348  -4256.789        -0.893       -0.950    60588.780       0    tf    2   E2F4p
#  13680 ONECUT2  -30259.974  -3628.359        -0.929       -0.923    39005.600       0    tf    2  KLF13p
#  14401 ONECUT2       0.000 -12694.647        -0.893       -0.871   422060.374       0    tf    5   KLF3p
#  17597 ONECUT2  -85931.740  -5531.252        -0.929       -0.981   106280.234       0    tf    1  NR2C2p
#  17933 ONECUT2    -929.791   -908.283        -0.964       -0.915     4603.984       0    tf    2  NR3C1p
#  20336 ONECUT2  -96324.656 -14718.732        -0.929       -0.973   563162.037       0    tf    5  RAD21p
#  24798 ONECUT2 -353296.600 -22589.160        -0.929       -0.969  1088384.369       0    tf    1   SMC3p
#  27995 ONECUT2  -27306.261 -15057.269        -0.929       -0.929  1033320.407       0    tf    4   TAL1p
#  28307 ONECUT2  -71149.130  -6332.839        -0.929       -0.951   226521.146       0    tf    1   TCF3p
#  29881 ONECUT2 -115512.680 -13958.887        -0.929       -0.963  1263314.250       0    tf    3 TRIM33p
#  30910 ONECUT2 -120781.618  -5162.044        -0.929       -0.976   259711.932       0    tf    1   USF1p
#  31287 ONECUT2       0.000  -9874.114        -0.893       -0.926   566378.703       0    tf    4  WDHD1p
#  31968 ONECUT2 -826000.677 -48000.351        -0.929       -0.981  3762868.834       0    tf    1 ZBTB7Ap

  #       gene   betaLasso  betaRidge spearmanCoeff pearsonCoeff     rfScore xgboost class rank  target tfbs
  #    ONECUT2       0.000 -80856.558        -0.893       -0.915 17188310.15       0    tf    4   CTCF     3
  #    ONECUT2   -6098.348  -4256.789        -0.893       -0.950    60588.78       0    tf    2   E2F4     4
  #    ONECUT2  -85931.740  -5531.252        -0.929       -0.981   106280.23       0    tf    1  NR2C2     1
  #    ONECUT2  -96324.656 -14718.732        -0.929       -0.973   563162.04       0    tf    5  RAD21     2
  #    ONECUT2 -353296.600 -22589.160        -0.929       -0.969  1088384.37       0    tf    1   SMC3     5
  #    ONECUT2  -27306.261 -15057.269        -0.929       -0.929  1033320.41       0    tf    4   TAL1     2
  #    ONECUT2  -71149.130  -6332.839        -0.929       -0.951   226521.15       0    tf    1   TCF3     3
  #    ONECUT2 -115512.680 -13958.887        -0.929       -0.963  1263314.25       0    tf    3 TRIM33     6

  #  reproduce one pearson, for example
  #   cor(mtx["CTCFp", late], mtx["ONECUT2", late], use="pairwise.complete") # -0.9146246

    data.dir <- "~/github/TrenaProjectErythropoiesis/prep/bigFimo/from-khaleesi"
    targets <- sub("p$", "", tbl.onecut2$target)
    counts <- list()
    for(target in targets){
       file <- sprintf("tbl.fimo.%s.RData", target)
       full.path <- file.path(data.dir, file)
       file.exists(full.path)
       tbl.fimo <- get(load(full.path))
       count <- nrow(subset(tbl.fimo, tf=="ONECUT2" & p.value < 1e-4))
       counts[[target]] <- count
       printf("%s: %d", target, count)
       }
   tbl.onecut2$TFBS.weak <- as.integer(counts)
   tbl.onecut2$TFBS.strong <- as.integer(counts)


   tbl.onecut2

     #        gene   betaLasso  betaRidge spearmanCoeff pearsonCoeff      rfScore xgboost class rank  target TFBS
     #  1  ONECUT2       0.000 -80856.558        -0.893       -0.915 17188310.146       0    tf    4   CTCFp    3
     #  2  ONECUT2   -6098.348  -4256.789        -0.893       -0.950    60588.780       0    tf    2   E2F4p    4
     #  3  ONECUT2  -30259.974  -3628.359        -0.929       -0.923    39005.600       0    tf    2  KLF13p    7
     #  4  ONECUT2       0.000 -12694.647        -0.893       -0.871   422060.374       0    tf    5   KLF3p    5
     #  5  ONECUT2  -85931.740  -5531.252        -0.929       -0.981   106280.234       0    tf    1  NR2C2p    1
     #  6  ONECUT2    -929.791   -908.283        -0.964       -0.915     4603.984       0    tf    2  NR3C1p   11
     #  7  ONECUT2  -96324.656 -14718.732        -0.929       -0.973   563162.037       0    tf    5  RAD21p    2
     #  8  ONECUT2 -353296.600 -22589.160        -0.929       -0.969  1088384.369       0    tf    1   SMC3p    5
     #  9  ONECUT2  -27306.261 -15057.269        -0.929       -0.929  1033320.407       0    tf    4   TAL1p    2
     #  10 ONECUT2  -71149.130  -6332.839        -0.929       -0.951   226521.146       0    tf    1   TCF3p    3
     #  11 ONECUT2 -115512.680 -13958.887        -0.929       -0.963  1263314.250       0    tf    3 TRIM33p    6
     #  12 ONECUT2 -120781.618  -5162.044        -0.929       -0.976   259711.932       0    tf    1   USF1p    1
     #  13 ONECUT2       0.000  -9874.114        -0.893       -0.926   566378.703       0    tf    4  WDHD1p    2
     #  14 ONECUT2 -826000.677 -48000.351        -0.929       -0.981  3762868.834       0    tf    1 ZBTB7Ap    0

    ctcfp.vec <- mtx["CTCFp", late]
    ctcfp.vec <- ctcfp.vec/max(ctcfp.vec)
    onecut2.vec <- mtx["ONECUT2", late]
    onecut2.vec <- onecut2.vec/max(onecut2.vec)
    plot(onecut2.vec, type="b", col="blue", main="ONECUT2 vs CTCFp, late timepoints")
    points(ctcfp.vec, type="b", col="red")

    legend(2, 0.6, c("ONECUT2 mRNA", "CTCFp srm"), c("blue", "red"))


} # onecut2.repression
