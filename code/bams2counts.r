########################################################################
### Rscript #1
### Input: BAM files produced by STAR, sample metadata file, GTF file
### Output: DGEList, raw count table (pre-filter, pre-normalization)
########################################################################

library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("edgeR")
library("Mus.musculus")

set.seed(20191021)

sampleTable <- read.csv("colData.csv", row.names=1)
filenames <- file.path("BAMs", paste0(rownames(sampleTable), ".bam"))
bamfiles <- BamFileList(filenames, yieldSize=2000000)
gtffile <- file.path("gencode.vM22.primary_assembly.annotation.gtf")
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
             genes=genetable,
             lib.size=coldata$seqDepth)


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

count <- data.frame(y$counts)
count$entrez <- mapIds(Mus.musculus,
                     keys=row.names(count),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

count$symbol <- mapIds(Mus.musculus,
                     keys=row.names(count),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

count$gene.name <- mapIds(Mus.musculus,
                     keys=row.names(count),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

write.csv(count, "countTable_rawLV.csv")

print("The original DGEList has been successfully generated.")
saveRDS(y, "bams2counts_raw.rds")


cutoff <- cpm(10, mean(y$samples$lib.size))
keep <- rowSums(cpm(y) > c(cutoff)) > 11
y.keep <- y[keep, ,keep.lib.sizes=FALSE]

count_filt <- data.frame(y.keep$counts)
count_filt$entrez <- mapIds(Mus.musculus,
                     keys=row.names(count_filt),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

count_filt$symbol <- mapIds(Mus.musculus,
                     keys=row.names(count_filt),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

count_filt$gene.name <- mapIds(Mus.musculus,
                     keys=row.names(count_filt),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

write.csv(count_filt, "countTable_filteredLV.csv")

print("Lowly expressed genes have been filtered.")
saveRDS(y.keep, "bams2counts_keep.rds")

sessionInfo()