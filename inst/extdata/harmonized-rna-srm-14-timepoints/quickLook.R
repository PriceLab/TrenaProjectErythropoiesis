print(load("mtx-rna-pkfm-27170x14.RData"))
print(load("mtx-srm-copyNumber-100x13.RData"))
dim(mtx.rna)  # 27170 14
dim(mtx.srm)  # 100   13

colnames(mtx.srm)
colnames(mtx.rna)

stopifnot(all(rownames(mtx.srm) %in% rownames(mtx.rna)))
stopifnot(all(colnames(mtx.srm) %in% colnames(mtx.rna)))

#----------------------------------------------------------------------------------------------------
# executed manually (30 jul 2021)
#
fix.protein.names <- function()
{

    bad.names <- setdiff(rownames(mtx.srm), rownames(mtx.rna))
    bad.names <- list("BCL11A_XL_L" = "BCL11A",
                      "CBP" = "CREBBP",
                      "coREST" = "RCOR1",
                      "E12/E47" = "TCF3",
                      "GR" = "NR3C1",
                      "HXB4" = "HOXB4",
                      "MLL3" = "KMT2C",
                      "MLL4 (KMT2D)" = "KMT2B",
                      "NC2B" = "DR1",
                      "PO2F1" = "POU2F1",
                      "RPB1" = "POLR2A",
                      "SET1B" = "SETD1B",
                      "SETB1" = "SETDB1",
                      "SIR6" = "SIRT6",
                      "SMCA4" = "SMARCA4",
                      "SMRC1" = "SMARCC1",
                      "SNF5" = "SMARCB1",
                      "STA5A" = "STAT5A",
                      "T2FA" = "GTF2F1",
                      "TF2B" = "GTF2F2",
                      "TF3C2" = "GTF3C2",
                      "UBF1" = "UBTF",
                      "UTX" = "KDM6A",
                      "ZC11A" = "ZC3H11A")

    all(as.character(bad.names) %in% rownames(mtx.rna))
    setdiff(as.character(bad.names), rownames(mtx.rna))
    all(names(bad.names) %in% rownames(mtx.srm))
    indices <- match(names(bad.names), rownames(mtx.srm))
    rownames(mtx.srm)[indices] <- as.character(bad.names)

    stopifnot(all(rownames(mtx.srm) %in% rownames(mtx.rna)))
    save(mtx.srm, file="mtx.srm.copyNumber-100x13.RData")

} # fix
#----------------------------------------------------------------------------------------------------
