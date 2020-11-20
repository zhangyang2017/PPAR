library("edgeR")
library("Mus.musculus")

set.seed(20191021)

mydata <- read.csv("filtered_LV.csv", row.names=1, header=TRUE)
sampleTable <- read.csv("colData.csv", row.names=1)
coldata <- data.frame(sampleTable)
genetable <- data.frame(gene.id=rownames(mydata))
countdata <- mydata


seGroups <- c(
    "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD", "LFD",
    "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", "HFD", 
    "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", "HXN", 
    "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN", "TXN")

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

keep <- filterByExpr(y, group=y$samples$Diet, 
    min.count=11, min.total.count=15)
y.keep <- y[keep, , keep.lib.sizes=FALSE]
dim(y.keep)

saveRDS(y.keep, "y.keep.rds")


count <- data.frame(y.keep$counts)
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

write.csv(count, "countTable_filteredLV.csv")

sessionInfo()
