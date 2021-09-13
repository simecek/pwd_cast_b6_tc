R_OPTS=--no-save --no-restore --no-init-file --no-site-file

all: preprocessing scanone scanone_outputs

preprocessing:
	Rscript $(R_OPTS) r/preprocessing_genotypes_mm.R
	Rscript $(R_OPTS) r/preprocessing_genotypes_gm.R
	Rscript $(R_OPTS) r/preprocessing_qtl_data.R
	rm Rplots.pdf
	
scanone:
	Rscript $(R_OPTS) r/qtl_scanone.R

scanone_outputs:
	Rscript $(R_OPTS) r/qtl_scanone_plots_and_tables.R
	
clean:
	rm -r data outputs
	mkdir data outputs outputs/scanone_individual_qtls
	