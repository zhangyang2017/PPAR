include config.mk

#CSV_FILES=$(wildcard $(DATADIR)/transcriptome*.csv)
RDA_FILES=$(wildcard ./*.rda)

##liverDGE:		detail
#liverDGE: figure1.rda plots combinePlots stats dependencyGraph.png

##>>transcriptome_LV.rda   :	Contains 7 objects - y, raw DGEList; y_filtered, genes removed; y_subset, leftover DGEList; y_keep, leftover DGEList after low expressed genes removal; y_keep_norm, normalized DGEList; y_keep_norm_disp, dispersion estimated DGEList, also the final DGEList; filtered_counts, unnormalized count table, can be exported as a csv file.
transcriptome_LV.rda: $(ALIGNMENT_SRC)
	$(ALIGNMENT_EXE)

##>>transcriptome_IM.rda   :	See transcriptome_LV.rda.
transcriptome_IM.rda: $(ALIGNMENT_SRC)
	$(ALIGNMENT_EXE)

##>>transcriptome_EWAT.rda :	See transcriptome_LV.rda.
transcriptome_EWAT.rda: $(ALIGNMENT_SRC)
	$(ALIGNMENT_EXE)	

##>>LVplots                :	Exploratory plots, includes library size distribution, normalization boxplots, and PCA/MD plots.	
LVplots: $(LVPLOTS_SRC) #transcriptome_LV.rda
	$(LVPLOTS_EXE)

##>>IMplots                :	Exploratory plots, includes library size distribution, normalization boxplots, and PCA/MD plots.	
IMplots: $(IMPLOTS_SRC) #transcriptome_LV.rda
	$(IMPLOTS_EXE)

##>>EWATplots              :	Exploratory plots, includes library size distribution, normalization boxplots, and PCA/MD plots.	
EWATplots: $(EWATPLOTS_SRC) #transcriptome_LV.rda
	$(EWATPLOTS_EXE)


##>>LV_DGEs                :	Detail.	
LV_DEGs: $(LV_DGEs_SRC) #transcriptome_LV.rda
	$(LV_DGEs_EXE)

##>>IM_DGEs                :	Detail.	
IM_DGEs: $(IM_DGEs_SRC) #transcriptome_LV.rda
	$(IM_DGEs_EXE)

##>>EWAT_DGEs                :	Detail.	
EWAT_DGEs: $(EWAT_DGEs_SRC) #transcriptome_LV.rda
	$(EWAT_DGEs_EXE)

##stdout                   :	detail.	
#stdout: $(STAT_SRC) figure1.rda
#	$(STAT_EXE)


##>>dependencyGraph.png    :	Dependency overview.
dependencyGraph.png: Makefile
	make -Bnd | make2graph | dot -Tpng -o $@

.PHONY: clean
clean:
	rm -f $(RDA_FILES)

.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<

#.PHONY: variables
#variables:
#	@echo CSV_FILES: $(CSV_FILES)