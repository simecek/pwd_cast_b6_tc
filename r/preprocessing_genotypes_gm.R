library(tidyverse)

cols = cols(
  .default = col_character(),
  Chromosome = col_character(),
  Position = col_character()
)
gm_data <- read_csv("raw_data/G_MUGA_geno.csv", col_types = cols)
  dim(gm_data)

# 2nd dataset excluded due to multiple probles (missing chrX, mice with CB genotype on chr17)
# used only for info about founders
    
gm_data2 <- read_csv("raw_data/GigaMUGA_Cast_Bara_07_15.csv", col_types=cols)
dim(gm_data2)
stopifnot(all(gm_data[,1] == gm_data2[,1]))
# 
# # there is missing info about chrX in G_MUGA1603_BF_JF.csv (most likely due to error)
# # but we can impute it from additional file
# gm_data2b <- read_csv("raw_data/G_MUGA_1506_BHP_BF.csv", col_types=cols)
# stopifnot(names(gm_data2b) == names(gm_data2))
# # columns to be replaced
# replace_cols = names(gm_data2)[grep("^[0-9]+$", names(gm_data2))]
# replace_cols
# # rows to be replaced
# replace_rows = which(gm_data2$Chromosome == "X")
# replace_from =  match(gm_data2$Name[replace_rows], gm_data2b$Name)
# summary(gm_data2$Name[replace_rows] == gm_data2b$Name[replace_from])
# gm_data2[replace_rows, replace_cols] <- gm_data2b[replace_from, replace_cols]
# TODO: recode imputed genotypes to match the pattern

# we have no preprocessed data for the last set of arrays
gm_data3 <- read_delim("raw_data/Wellcome_Trust_Centre_MURGIGV01_20161113/Wellcome_Trust_Centre_MURGIGV01_20161113_FinalReport.txt",
                       skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype) %>%
  select(-starts_with("FATW"))  %>%  # remove FATW columns (samples that are not ours)
  filter(`SNP Name` %in% gm_data$Marker)
dim(gm_data3)

# arrange data to be in the same order as gm_data
gm_data3 <- gm_data3[match(gm_data[[1]], gm_data3[[1]]),]
stopifnot(all(gm_data[,1] == gm_data3[,1]))
dim(gm_data3)

# Adding map and founder genotypes ----------------------------------------

col_types = cols(
  marker = col_character(),
  chr = col_character(),
  bp_mm10 = col_double(),
  cM_cox = col_double(),
  cM_g2f1 = col_double(),
  strand = col_character(),
  snp = col_character(),
  unique = col_logical(),
  multi = col_logical(),
  unmapped = col_logical(),
  n_blast_hits = col_double(),
  n_blast_chr = col_double(),
  probe = col_character()
)

snp_map <- read_csv("raw_data/gm_uwisc_v1.csv", col_types = col_types) %>%
  select(marker, chr, bp_mm10, cM_g2f1, strand, unique, unmapped)
all(gm_data$Marker %in% snp_map$marker)
all(!duplicated(gm_data$Marker))

founders <- gm_data2[, c("Name", "PWD", "B6", "CAST")] %>%
  transmute(Marker = Name,
         `CAST/EiJ` = toupper(`CAST`),
         `PWD/PhJ` = toupper(`PWD`),
         `C57BL/6J` = toupper(`B6`))
table(gm_data3$`SNP Name` %in% founders$Marker)
table(founders$Marker %in% gm_data3$`SNP Name`)

mapped_data = founders %>%
  filter(`PWD/PhJ` != `CAST/EiJ`) %>%  #informative
  filter(`PWD/PhJ` != 'N', `CAST/EiJ` != 'N', `C57BL/6J` != "N") %>% #non-mising
  left_join(snp_map, by=c("Marker" = "marker")) %>%
  left_join(gm_data3, by=c("Marker" = "SNP Name"))

# Recoding ----------------------------------------------------------------

CB1 = paste0(mapped_data$`CAST/EiJ`, mapped_data$`C57BL/6J`)
CB2 = paste0(mapped_data$`C57BL/6J`, mapped_data$`CAST/EiJ`)
PB1 = paste0(mapped_data$`PWD/PhJ`, mapped_data$`C57BL/6J`)
PB2 = paste0(mapped_data$`C57BL/6J`, mapped_data$`PWD/PhJ`)
BB = paste0(mapped_data$`C57BL/6J`, mapped_data$`C57BL/6J`)
PP = paste0(mapped_data$`PWD/PhJ`, mapped_data$`PWD/PhJ`)
CC = paste0(mapped_data$`CAST/EiJ`, mapped_data$`CAST/EiJ`)

genotype_data = select(mapped_data, -unique, -unmapped, -strand)
male_cols = names(genotype_data)[grep("^[0-9]+$", names(genotype_data))]
autosomes = as.character(1:19)

for (c in male_cols) {
  new_geno = rep("N", nrow(genotype_data))
  
  cb_idx = (mapped_data[[c]] == CB1 | mapped_data[[c]] == CB2) & (mapped_data$chr %in% autosomes)
  new_geno[cb_idx] <- "CB"
  pb_idx = (mapped_data[[c]] == PB1 | mapped_data[[c]] == PB2) & (mapped_data$chr %in% autosomes)
  new_geno[pb_idx] <- "PB"
  
  c_idx = (mapped_data[[c]] == CC) & (mapped_data$chr == "X")
  new_geno[c_idx] <- "C"
  p_idx = (mapped_data[[c]] == PP) & (mapped_data$chr == "X")
  new_geno[p_idx] <- "P"
  
  # Y should be B6
  b_idx = (mapped_data[[c]] == BB) & (mapped_data$chr == "Y")
  new_geno[b_idx] <- "B"
  # M shoud be PWD
  p_idx = (mapped_data[[c]] == PP) & (mapped_data$chr == "M")
  new_geno[p_idx] <- "P"
  
  genotype_data[c] = new_geno
}

head(genotype_data)
sort(apply(genotype_data[,male_cols] == "N", 2, mean))  # 76879 is clearly off and should be removed
genotype_data <- select(genotype_data, -`76879`)
names(genotype_data)
stopifnot(all(gm_data[,1] == genotype_data[,1]))


# Join all GM data into one table -----------------------------------------

dim(gm_data)
dim(gm_data2)
dim(genotype_data)
#gm_joined <- cbind(gm_data, select(gm_data2, matches("^[0-9]+$")), select(genotype_data, matches("^[0-9]+$")))
gm_joined <- cbind(gm_data, select(genotype_data, matches("^[0-9]+$")))
dim(gm_joined)



# Final preprocessing -----------------------------------------------------

`%notin%` <- Negate(`%in%`)
output_data = gm_joined %>%
  select(Marker, `CAST/EiJ`, `PWD/PhJ`, `C57BL/6J`) %>%
  left_join(snp_map, by=c("Marker" = "marker")) %>%
  filter(!is.na(chr) & !is.na(bp_mm10) & (chr!="PAR")) %>%
  left_join(gm_joined[,names(gm_joined) %notin% c("CAST/EiJ", "PWD/PhJ", "C57BL/6J", "Chromosome", "Position")], by=c("Marker" = "Marker"))

# reorder
output_data[['chr2']] <- as.numeric(output_data$chr)  # warning expected
output_data$chr2[output_data$chr == "X"] = 20
output_data$chr2[output_data$chr == "Y"] = 21
output_data$chr2[output_data$chr == "M"] = 22

output_data = output_data %>% arrange(chr2, bp_mm10)

# save

output_data = select(output_data, -unique, -unmapped, -strand, -chr2)

sample_cols = names(output_data)[grep("^[0-9]+$", names(output_data))]
pct_n = apply(output_data[,sample_cols] == "N", 1, mean)
hist(pct_n)
table(pct_n > 0.33)
output_data <- output_data[pct_n <= 0.33, ]
table(output_data$chr)

write_csv(output_data, "data/genotypes_gm.csv")
