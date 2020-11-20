library("edgeR")
library("Mus.musculus")

set.seed(20191021)

mydata <- read.csv("countMatrix_final.csv", row.names=1, header=TRUE)
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

saveRDS(y, "finalCount.rds")

print("The final DGEList has been successfully generated.")

sessionInfo()