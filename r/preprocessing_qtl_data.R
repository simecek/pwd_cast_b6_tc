library(tidyverse)
library(qtl)


# Joining MM a GM genotypes -----------------------------------------------

mm_data <- read_csv("data/genotypes_mm.csv", col_types=cols(chr = col_character())) %>%
  filter(chr %in% c(as.character(1:19), "X")) %>%
  select(-`CAST/EiJ`, -`PWD/PhJ`, -`C57BL/6J`)
dim(mm_data)

gm_data <- read_csv("data/genotypes_gm.csv", col_types=cols(chr = col_character())) %>%
  filter(chr %in% c(as.character(1:19), "X")) %>%
  select(-`CAST/EiJ`, -`PWD/PhJ`, -`C57BL/6J`)
dim(gm_data)

# revise GigaMUGA &	MiniMUGA cols in raw_data/SU
info <- read_csv("raw_data/SUMMARY_TW_SC_ASY_3_7_8K.csv")
info$MiniMUGA <- c('', 'miniMUGA')[as.numeric(as.character(info$ID) %in% names(mm_data))+1] 
info$GigaMUGA <- c('', 'gigaMUGA')[as.numeric(as.character(info$ID) %in% names(gm_data))+1] 
info$Hstx2 <- toupper(info$Hstx2)
info$Prdm9 <- toupper(info$Prdm9)
write_csv(info, "data/SUMMARY_TW_SC_ASY_3_7_8K_MUGAedited.csv")

# remove 7 samples that are duplicated with MiniMuga
gm_data <- select(gm_data, -one_of(c("78825", "76999", "76865", "76867", "76882", "76892", "76894")))

# keep only markers on both arrays
small_data <- inner_join(mm_data, gm_data, by = c("chr" = "chr", "bp_mm10" = "bp_mm10")) %>%
  select(-bp_mm10, -Marker.y, -cM_g2f1.y)
dim(small_data)
names(small_data)[1] <- "Marker"
names(small_data)
table(small_data$chr)

# reshape data into qtl csvsr format
names(small_data)[1:3] <- c("id", "", "")
small_data <- small_data[!duplicated(small_data$id),]  # exclude 24 cases of duplication

sample_cols = names(small_data)[grep("^[0-9]+$", names(small_data))]
small_data[sample_cols][small_data[sample_cols] == "PB"] <- "H"
small_data[sample_cols][small_data[sample_cols] == "CB"] <- "A"
small_data[sample_cols][small_data[sample_cols] == "P"] <- "P"
small_data[sample_cols][small_data[sample_cols] == "C"] <- "A"
table(small_data[sample_cols][,1])

# Checking phenotypes -----------------------------------------------------

pheno <- read_csv("raw_data/SUMMARY_TW_SC_ASY_3_7_8K.csv")

# those are missing
missing_pheno_ids <- names(small_data)[-(1:3)][!(as.numeric(names(small_data)[-(1:3)]) %in% pheno$ID)]
missing_pheno_ids
# mouse 36377 has been excluded due to an extremly weird numbers (TW>200, ASY>40, SC>...)

# samples with missing phenotypes should be excluded
dim(small_data)
small_data <- small_data[,!(names(small_data) %in% missing_pheno_ids)]
dim(small_data)
write_csv(small_data, "data/genotypes_qtl.csv")

pheno_sel <- pheno[match(names(small_data)[-(1:3)], pheno$ID),]
names(pheno_sel)[1] <- "id"
names(pheno_sel)[names(pheno_sel) == "ASY%"] <- "ASY"
dim(pheno_sel)
pheno_sel$infertility_cat <- as.numeric((pheno_sel$TW < 80) & (pheno_sel$SC < 5))

write_csv(pheno_sel, "data/phenotypes.csv")

# Transforming phenotypes for QTL mapping ---------------------------------

# TW
rpheno <- t(pheno_sel %>%
              select(TW, id)) 
write.table(as.data.frame(rpheno), "data/pheno_TW.csv", col.names = FALSE, row.names=TRUE, sep=",")

#SC
rpheno <- pheno_sel  %>%
  select(SC, id) 
rpheno$SC <- log(1+rpheno$SC)
write.table(as.data.frame(t(rpheno)), "data/pheno_SC.csv", col.names = FALSE, row.names=TRUE, sep=",")

#Infertility (cat.)
rpheno <- pheno_sel  %>%
  select(infertility_cat, id) 
write.table(as.data.frame(t(rpheno)), "data/pheno_infertility_cat.csv", col.names = FALSE, row.names=TRUE, sep=",")

# ASY
rpheno <- pheno_sel %>%
  select(ASY, id) 
rpheno$ASY <- log(rpheno$ASY / (100 - rpheno$ASY))
write.table(as.data.frame(t(rpheno)), "data/pheno_ASY.csv", col.names = FALSE, row.names=TRUE, sep=",", na='')


