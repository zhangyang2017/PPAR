################################################################################
### Liver
### Input: BAM files produced by STAR, sample metadata file, GTF file
### Output: DGEList, count table 
################################################################################
library(knitr)
.cran_packages <- c("Rsamtools", "GenomicFeatures", "GenomicAlignments", "edgeR", "Mus.musculus")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
################################################################################
set.seed(20191103)

sampleTable <- read.csv("colData.csv", row.names=1)
filenames <- file.path("BAMs", paste0(rownames(sampleTable), ".bam"))
bamfiles <- BamFileList(filenames, yieldSize=2000000)
gtffile <- file.path("./gencode.vM22.primary_assembly.annotation.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")
se <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE)


colData(se) <- DataFrame(sampleTable)
row.names(se) = gsub("\\..*", "", row.names(se))
genetable <- data.frame(gene.id=rownames(se))
countdata <- assay(se)
coldata <- colData(se)

seGroups <- c(
    "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", 
    "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", 
    "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN",
    "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN")

y <- DGEList(counts=countdata, 
             samples=coldata, 
             group=factor(seGroups),
             genes=genetable)

y$genes$entrez <- mapIds(Mus.musculus,
                     keys=row.names(y$genes),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$symbol <- mapIds(Mus.musculus,
                     keys=row.names(y$genes),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$gene.name <- mapIds(Mus.musculus,
                     keys=row.names(y$genes),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

toRemove <- c("ribosomal protein", "RIKEN", "predicted gene", "pseudogene",
			  "cDNA sequence", "nuclear encoded rRNA", "microRNA", 
			  "DNA segment", "uncharacterized", "pre-RNA", "20S rRNA",
			  "non-protein coding")

yy <- y$genes[!is.na(y$genes$gene.name), ]
yy <- yy[!grepl(paste(toRemove, collapse = "|"), yy$gene.name), ]
y_subset <- y[rownames(yy), ]

y_filtered <- y[ !(y$genes$gene.name %in% y_subset$genes$gene.name), ]


keep <- filterByExpr(y, group=y$samples$Diet, min.total.count=9)
y_keep <- y[keep, , keep.lib.sizes=FALSE]

counts <- data.frame(y_keep$counts)
counts$entrez <- mapIds(Mus.musculus,
                     keys=row.names(counts),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
counts$symbol <- mapIds(Mus.musculus,
                     keys=row.names(counts),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
counts$gene.name <- mapIds(Mus.musculus,
                     keys=row.names(counts),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

write.csv(counts, "./countTable_LV.csv")
saveRDS(y_keep, "y_keep.rds")

