library(tidyverse)
library(qtl)
library(WriteXLS)
set.seed(42) # to get the same plots

alltraits_qtl_table <- NULL

# TW ----------------------------------------------------------------------

load("data/scanone_TW.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="TW")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_TW.pdf", width=8, height=5)

significant_qtl_table <- function(trait_name, inverse.transform = function(x) x, categorical = FALSE) {
  critical.lod <- summary(bcp.perm)[[1]]
  significant.chrs <- unique(bcp.scanone$chr[bcp.scanone$lod > critical.lod])
  pheno = inverse.transform(bcp$pheno[,1])
  
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
    
    P_average[i] = mean(pheno[gener.aa==0], na.rm = TRUE)
    C_average[i] = mean(pheno[gener.aa==1], na.rm = TRUE)
    P_sem[i] = sd(pheno[gener.aa==0], na.rm = TRUE) / sqrt(sum(gener.aa==0 & !is.na(pheno)))
    C_sem[i] = sd(pheno[gener.aa==1], na.rm = TRUE) / sqrt(sum(gener.aa==1 & !is.na(pheno)))
    
    # additional plots and tables
    chr_markers <- bcp.scanone[bcp.scanone$chr==chrs[i] & !grepl("^c[0-9XY]*[.]loc[0-9]*[ ]*", rownames(bcp.scanone)),]
    best_marker <- rownames(chr_markers)[which.max(chr_markers$lod)[1]]
    pxg_data <- plotPXG(bcp, best_marker)
    
    # pxg plot
    img_file_name = paste0(trait_name, "_chr", chrs[i], "_", best_marker)
    dev.copy2pdf(file=paste0("outputs/scanone_individual_qtls/", img_file_name, ".pdf"), width=8, height=5)
    
    # table
    marker_data = tibble(id = bcp$pheno$id, geno = c("C", "P")[as.numeric(pxg_data[,1])], 
                         pheno = inverse.transform(pxg_data$pheno), 
                         inferred = pxg_data$inferred)
    marker_data = arrange(marker_data, geno, pheno)
    names(marker_data)[2] <- best_marker
    names(marker_data)[3] <- trait_name
    write_csv(marker_data, paste0("outputs/scanone_individual_qtls/", img_file_name, ".csv"))
    WriteXLS(marker_data, paste0("outputs/scanone_individual_qtls/", img_file_name, ".xls"))  
  }
  
  tibble(trait = trait_name, Chr = chrs, Position = position, LOD = LOD, CI_left=ci_left, CI_right=ci_right, 
        C_mean=C_average, C_sem=C_sem, P_mean=P_average, P_sem=P_sem)
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

significant_qtl_table("SC", inverse.transform = function(x) exp(x)-1)
alltraits_qtl_table <- rbind(alltraits_qtl_table, 
                             significant_qtl_table("SC", inverse.transform = function(x) exp(x)-1))


# ASY ----------------------------------------------------------------------

load("data/scanone_ASY.rdata")

summary(bcp.perm)
plot(bcp.scanone, ylim = c(0, max(bcp.scanone$lod, summary(bcp.perm)[[1]])), main="ASY%")
abline(h=summary(bcp.perm)[[1]], col="red", lty=2)
abline(h=summary(bcp.perm)[[2]], col="blue", lty=2)

dev.copy2pdf(file="outputs/scanone_ASY.pdf", width=8, height=5)

asy.tranform.inverse = function(x) round(exp(x) * 100 / (1 + exp(x)), 3)
significant_qtl_table("ASY", inverse.transform = asy.tranform.inverse)
alltraits_qtl_table <- rbind(alltraits_qtl_table, 
                             significant_qtl_table("ASY", inverse.transform = asy.tranform.inverse))


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

# convert position from cM to bp
source("r/utils/cm2bp.R")
alltraits_qtl_table$Position <- round(cm2bp(alltraits_qtl_table$Chr, alltraits_qtl_table$Position))

write_csv(alltraits_qtl_table, "outputs/scanone_qtl_table.csv")
WriteXLS(alltraits_qtl_table, "outputs/scanone_qtl_table.xls")

# TODOs
# 4) SEM for categorical variables  
