library(qtl)

bcp <- read.cross(format = "csvsr", dir="data/", "genotypes_qtl.csv", "pheno_TW.csv", genotypes=c("A","H","P","D","C"), na="N")
bcp <- jittermap(bcp)
bcp <- calc.genoprob(bcp, step=1, error.prob=0.01)
bcp.stepwise <- stepwiseqtl(bcp, max.qtl=3, method="hk", keeptrace = TRUE)

stepwiseqtl