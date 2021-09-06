# pwd_cast_b6_tc
QTL analysis for (PWDxCAST)xB6 TC1 paper

## How to reproduce the analysis

 *  
 * Install `R` and its packages `qtl`, `tidyverse` and `WriteXLS`
 * Download raw data and save them into `raw_data` folder (raw data will be available before the paper submission)
 * Either `cd` into the project folder and call `make all` or run the R scripts in the following order:
 
    1. `preprocessing_genotypes_mm.R`
    1. `preprocessing_genotypes_gm.R`

##  