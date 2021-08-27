library(tidyverse)
library(qtl)

alltraits_qtl_table <- NULL

# TW ----------------------------------------------------------------------

load("data/scanone_TW.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="TW")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_TW.pdf", width=8, height=5)

significant_qtl_table <- function(trait_name, categorical = FALSE) {
  critical.lod <- summary(bcp.perm)[[1]]
  significant.chrs <- unique(bcp.scanone$chr[bcp.scanone$lod > critical.lod])
  pheno = bcp$pheno[,1]
  
  N <- length(significant.chrs)
  position <- ci_left <- ci_right <- LOD <- C_average <- C_sem <- P_average <- P_sem <- rep(NA, N)
  chrs <- rep("", N)
  
  for (i in seq_along(significant.chrs)) {
    
    chrs[i] = as.character(significant.chrs[i])
    ci = bayesint(bcp.scanone, chr=chrs[i])
    ci_left[i] = ci$pos[1]
    position[i] = ci$pos[2]
    LOD[i] = ci$lod[2]
    ci_right[i] = ci$pos[3]
    
    best_loc = sub("^c[0-9XY]*[.](loc[0-9]*)[ ]*", "\\1", rownames(ci)[2]) # remove c[CHR]. from loc name
    probs.aa = bcp[["geno"]][[as.character(significant.chrs[i])]][["prob"]][, best_loc, 1]
    gener.aa = rbinom(rep(1,length(probs.aa)), 1, prob = probs.aa)
    
    C_average[i] = mean(pheno[gener.aa==0], na.rm = TRUE)
    P_average[i] = mean(pheno[gener.aa==1], na.rm = TRUE)
    C_sem[i] = sd(pheno[gener.aa==0], na.rm = TRUE) / sqrt(sum(gener.aa==0 & !is.na(pheno)))
    P_sem[i] = sd(pheno[gener.aa==1], na.rm = TRUE) / sqrt(sum(gener.aa==1 & !is.na(pheno)))
  }
  
  tibble(trait = trait_name, chrs = chrs, position = position, LOD = LOD, ci_left=ci_left, ci_right=ci_right, 
        C_average=C_average, C_sem=C_sem, P_average=P_average, P_sem=P_sem)
}

significant_qtl_table("TW")
alltraits_qtl_table <- rbind(alltraits_qtl_table, significant_qtl_table("TW"))

# SC  ----------------------------------------------------------------------

load("data/scanone_SC.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="logSC")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_SC.pdf", width=8, height=5)

significant_qtl_table("logSC")
alltraits_qtl_table <- rbind(alltraits_qtl_table, significant_qtl_table("logSC"))


# ASY ----------------------------------------------------------------------

load("data/scanone_ASY.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="ASY%")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_ASY.pdf", width=8, height=5)

significant_qtl_table("ASY%")
alltraits_qtl_table <- rbind(alltraits_qtl_table, significant_qtl_table("ASY%"))


# INFERTILITY (cat.) ------------------------------------------------------

load("data/scanone_infertility_cat.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="Infertility (cat.)")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_infertility_cat.pdf", width=8, height=5)

significant_qtl_table("Infertility (cat.)")
alltraits_qtl_table <- rbind(alltraits_qtl_table, significant_qtl_table("Infertility (cat.)"))


# QTL table ---------------------------------------------------------------

write_csv(alltraits_qtl_table, "outputs/scanone_qtl_table.csv")
# TODOs
# 1) Make sure the P & C are in the right direction
# 2) remap cM into Mb
# 3) mean / sem summary on transformed on untransformed level?
# 4) SEM for categorical variables
