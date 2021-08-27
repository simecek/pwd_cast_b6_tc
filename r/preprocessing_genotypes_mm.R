library(tidyverse)

# joining MiniMuga experiments --------------------------------------------------------

array1 <- read_delim("raw_data/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201117/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201117/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201117_FinalReport.txt",
                     skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype)
dim(array1)

array2 <- read_delim("raw_data/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201130/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201130/Inst_of_Mol_Genetics_Forejt_MURCOMV02_20201130_FinalReport.txt",
                     skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype)
dim(array2)

array3 <- read_delim("raw_data/IMG_Forejt_MURCOMV02_20201202/IMG_Forejt_MURCOMV02_20201202/IMG_Forejt_MURCOMV02_20201202_FinalReport.txt",
                     skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype)
dim(array3)

array4 <- read_delim("raw_data/Inst_of_Molecular_Forejt_MURCOMV02_20210624/Inst_of_Molecular_Forejt_MURCOMV02_20210624/Inst_of_Molecular_Forejt_MURCOMV02_20210624_FinalReport.txt",
                     skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype)
dim(array4)

array5 <- read_delim("raw_data/Inst_of_Molecular_Forejt_MURCOMV02_20210803/Inst_of_Molecular_Forejt_MURCOMV02_20210803_FinalReport.txt",
                     skip=9, delim="\t", col_types=cols(`Sample ID` = col_character())) %>%
  mutate(genotype = paste0(`Allele1 - Forward`, `Allele2 - Forward`)) %>%
  select(`SNP Name`, `Sample ID`, genotype) %>%
  pivot_wider(names_from = `Sample ID`, values_from = genotype)

stopifnot(all(array1$`SNP Name` == array2$`SNP Name`))
stopifnot(all(array1$`SNP Name` == array3$`SNP Name`))
stopifnot(all(array1$`SNP Name` == array4$`SNP Name`))
stopifnot(all(array1$`SNP Name` == array5$`SNP Name`))

raw_data = array1 %>% 
  left_join(array2, by="SNP Name") %>%
  left_join(array3, by="SNP Name") %>%
  left_join(array4, by="SNP Name") %>%
  left_join(array5, by="SNP Name")

dim(raw_data)
names(raw_data)

raw_data = select(raw_data, !contains("fem"))
dim(raw_data)
names(raw_data)
sort(apply(raw_data == "--", 2, mean))  # fraction of missing

# Adding map and founder genotypes ----------------------------------------

# downloaded from https://github.com/kbroman/MUGAarrays
snp_map = read_csv("raw_data/mini_uwisc_v2.csv", col_types=cols(chr = col_character())) %>%
  select(marker, chr, bp_mm10, cM_g2f1, strand, unique, unmapped)
stopifnot(all(raw_data$`SNP Name` %in% snp_map$marker))

# downloaded from https://www.med.unc.edu/mmrrc/genotypes/resources/
founders = read_csv("raw_data/miniMUGA-Consensus-Genotypes.csv", col_types=cols(Chromosome = col_character()))
table(raw_data$`SNP Name` %in% founders$Marker)
stopifnot(all(founders$Marker %in% raw_data$`SNP Name`))

founders$`CAST/EiJ` = toupper(founders$`CAST/EiJ`)
founders$`PWD/PhJ` = toupper(founders$`PWD/PhJ`)
founders$`C57BL/6J` = toupper(founders$`C57BL/6J`)

mapped_data = founders %>%
  select(Marker, `CAST/EiJ`, `PWD/PhJ`, `C57BL/6J`) %>%
  filter(`PWD/PhJ` != `CAST/EiJ`) %>%  #informative
  filter(`PWD/PhJ` != 'N', `CAST/EiJ` != 'N', `C57BL/6J` != "N") %>% #non-mising
  left_join(snp_map, by=c("Marker" = "marker")) %>%
  filter(!is.na(chr) & !is.na(bp_mm10) & (chr!="PAR")) %>%
  left_join(raw_data, by=c("Marker" = "SNP Name"))

mapped_data$chr2 <- as.numeric(mapped_data$chr)  # warning expected
mapped_data$chr2[mapped_data$chr == "X"] = 20
mapped_data$chr2[mapped_data$chr == "Y"] = 21
mapped_data$chr2[mapped_data$chr == "M"] = 22

mapped_data = mapped_data %>% arrange(chr2, bp_mm10)


# Recoding ----------------------------------------------------------------

CB1 = paste0(mapped_data$`CAST/EiJ`, mapped_data$`C57BL/6J`)
CB2 = paste0(mapped_data$`C57BL/6J`, mapped_data$`CAST/EiJ`)
PB1 = paste0(mapped_data$`PWD/PhJ`, mapped_data$`C57BL/6J`)
PB2 = paste0(mapped_data$`C57BL/6J`, mapped_data$`PWD/PhJ`)
BB = paste0(mapped_data$`C57BL/6J`, mapped_data$`C57BL/6J`)
PP = paste0(mapped_data$`PWD/PhJ`, mapped_data$`PWD/PhJ`)
CC = paste0(mapped_data$`CAST/EiJ`, mapped_data$`CAST/EiJ`)

genotype_data = select(mapped_data, -unique, -unmapped, -strand, -chr2)
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

# Filter probes -----------------------------------------------------------

# filter non-functional markers
pct_n = apply(genotype_data[,male_cols] == "N", 1, mean)
hist(pct_n)
table(pct_n > 0.1)
genotype_data <- genotype_data[pct_n <= 0.1, ]
dim(genotype_data)


# Saving MiniMUGA genotypes -----------------------------------------------

write_csv(genotype_data, "data/genotypes_mm.csv")
