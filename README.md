# pwd_cast_b6_tc
QTL analysis for (PWDxCAST)xB6 TC1 paper

## Subfolders

  * [raw_data](raw_data/): original data obtained from Jiri Forejt (**raw data subfolder will be available after the paper submission**).
  * [r](r/): R scripts used for data preprocessing and analysis
  * [data](data/): processed data, summary tables and RData files
  * [outputs](outputs/): figures and tables used in the publication

## How to reproduce the analysis

 * Clone this repository. 
 * Install `R` and its packages `qtl`, `tidyverse` and `WriteXLS`.
 * Download raw data and save them into `raw_data` folder (**raw data will be available after the paper submission**).
 * Either `cd` into the project folder and call `make all` or run the R scripts in the following order:
 
    1. `preprocessing_genotypes_mm.R`
    1. `preprocessing_genotypes_gm.R`
    1. `preprocessing_qtl_data.R`
    1. `qtl_scanone.R`
    1. `qtl_scanone_plots_and_tables.R`
    
    
