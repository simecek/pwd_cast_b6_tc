library(tidyverse)
library(qtl)


# TW ----------------------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_TW.csv", genotypes=c("A","H","P","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scantwo <- scantwo(bcp, method="hk")
bcp.perm <- scantwo(bcp, method="hk", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scantwo, bcp.perm, file = "data/scantwo_TW.rdata")

# SC  ----------------------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_SC.csv", genotypes=c("A","H","P","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scantwo <- scantwo(bcp, method="hk")
bcp.perm <- scantwo(bcp, method="hk", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scantwo, bcp.perm, file = "data/scantwo_SC.rdata")

# ASY ----------------------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_ASY.csv", genotypes=c("A","H","P","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scantwo <- scantwo(bcp, method="hk")
bcp.perm <- scantwo(bcp, method="hk", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scantwo, bcp.perm, file = "data/scantwo_ASY.rdata")

# INFERTILITY (cat.) ------------------------------------------------------

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_infertility_cat.csv", genotypes=c("A","H","P","D","C"), na="N")
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.scantwo <- scantwo(bcp, method="hk", model="binary")
bcp.perm <- scantwo(bcp, method="hk", model="binary", n.perm=1000, n.cluster = 7)

save(bcp, bcp.scantwo, bcp.perm, file = "data/scantwo_infertility_cat.rdata")


