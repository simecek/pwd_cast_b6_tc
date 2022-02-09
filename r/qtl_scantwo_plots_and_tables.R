library(tidyverse)
library(qtl)
library(WriteXLS)
set.seed(42) # to get the same plots

alltraits_qtl_table <- NULL

# TW ----------------------------------------------------------------------

load("data/scantwo_TW.rdata")

pdf(file="outputs/scantwo_TW.pdf", width=8, height=5)
plot(bcp.scantwo, main="TW")
dev.off()

tmp <- summary(bcp.scantwo, perms=bcp.perm, pvalues=TRUE, alphas=c(0.05, 0.10, 0.10, 0.05, 0.10))
trait_table <- cbind(trait="TW", as.data.frame(tmp))
alltraits_qtl_table <- rbind(alltraits_qtl_table, trait_table)

# SC  ----------------------------------------------------------------------

load("data/scantwo_SC.rdata")

pdf(file="outputs/scantwo_SC.pdf", width=8, height=5)
plot(bcp.scantwo, main="SC")
dev.off()

tmp <- summary(bcp.scantwo, perms=bcp.perm, pvalues=TRUE, alphas=c(0.05, 0.10, 0.10, 0.05, 0.10))
trait_table <- cbind(trait="SC", as.data.frame(tmp))
alltraits_qtl_table <- rbind(alltraits_qtl_table, trait_table)


# ASY ----------------------------------------------------------------------

load("data/scantwo_ASY.rdata")

pdf(file="outputs/scantwo_ASY.pdf", width=8, height=5)
plot(bcp.scantwo, main="ASY")
dev.off()

tmp <- summary(bcp.scantwo, perms=bcp.perm, pvalues=TRUE, alphas=c(0.05, 0.10, 0.10, 0.05, 0.10))
trait_table <- cbind(trait="ASY", as.data.frame(tmp))
alltraits_qtl_table <- rbind(alltraits_qtl_table, trait_table)

# QTL table ---------------------------------------------------------------

# convert position from cM to bp
source("r/utils/cm2bp.R")
alltraits_qtl_table$pos1f <- round(cm2bp(alltraits_qtl_table$chr1, alltraits_qtl_table$pos1f))
alltraits_qtl_table$pos2f <- round(cm2bp(alltraits_qtl_table$chr2, alltraits_qtl_table$pos2f))

write_csv(alltraits_qtl_table, "outputs/scantwo_qtl_table.csv")
WriteXLS(alltraits_qtl_table, "outputs/scantwo_qtl_table.xls")

# TODOs
# 4) SEM for categorical variables  
