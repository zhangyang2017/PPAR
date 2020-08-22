# Makefile
# Zhang-2019-TXN
# Yang Zhang
# 2019 November

include config.mk

CSV_FILES=$(wildcard $(DATADIR)/phenome_*.csv)
RDA_FILES=$(wildcard ./*.rda)

##justFigure housekeeping.R plots.R combinePlots.R
justFigure: figure1.rda plots combinePlots

##result1:		housekeeping.R plots.R combinePlots.R stats.R
result1: figure1.rda plots combinePlots stats #dependencyGraph.png

##figure1.rda:		Copy related raw excel files for result 1, clean up raw data and store cleaned data in rda file for downstream analyses. Copied raw files will be removed.
figure1.rda: $(HOUSEKEEPING_SRC) $(CSV_FILES)
	cp $(DATADIR)/*.csv .
	$(HOUSEKEEPING_EXE)
	rm -f *.csv


##plots:			Generate individual plot for figure 1.	
plots: $(PLOTS_SRC) figure1.rda
	$(PLOTS_EXE)


##combinePlots:		Combine individual plots into final figure 1.	
combinePlots: $(COMBINE_SRC) $(RDA_FILES)
	$(COMBINE_EXE)
	rm -f Rplots.pdf


##analysis_output.txt:	Summary statistics and contrast comparison result table.	
stats: $(STAT_SRC) figure1.rda
	$(STAT_EXE)


##dependencyGraph.png:	Dependency overview.
dependencyGraph.png: Makefile
	make -Bnd | make2graph | dot -Tpng -o $@

.PHONY: clean
clean:
	rm -f $(RDA_FILES) figure1.svg figure1.pdf analysis_output.txt dependencyGraph.png

.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: variables
variables:
	@echo CSV_FILES: $(CSV_FILES)
