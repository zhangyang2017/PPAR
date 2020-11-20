## test for DE genes
###############################################################################
options(java.parameters = "-Xmx8g")
library(knitr)
.cran_packages <- c("ggplot2", "xlsx", "RColorBrewer", "edgeR", "pryr")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
###############################################################################

load("transcriptome_LV.rda")

## Setting robust=TRUE in glmQLFit is usually recommended (Phipson et al. 2016). 

design <- model.matrix(~ 0 + y_keep_norm_disp$samples$Diet)
colnames(design) <- levels(y_keep_norm_disp$samples$Diet)

fit <- glmQLFit(y_keep_norm_disp, design, robust=TRUE)
con <- makeContrasts(TXN - HFD, levels=design)
con2 <- makeContrasts(HXN - HFD, levels=design)
con3 <- makeContrasts(LFD - HFD, levels=design)
con4 <- makeContrasts(TXN - HXN, levels=design)
con5 <- makeContrasts(TXN - LFD, levels=design)

TXN.vs.HFD <- glmQLFTest(fit, contrast=con)
TXN.vs.LFD <- glmQLFTest(fit, contrast=con5)
HXN.vs.HFD <- glmQLFTest(fit, contrast=con2)
LFD.vs.HFD <- glmQLFTest(fit, contrast=con3)
TXN.vs.HXN <- glmQLFTest(fit, contrast=con4)

stat_txnHFD_LV <- topTags(TXN.vs.HFD, n=nrow(y_keep_norm_disp))$table
stat_txnLFD_LV <- topTags(TXN.vs.LFD, n=nrow(y_keep_norm_disp))$table
stat_hxnHFD_LV <- topTags(HXN.vs.HFD, n=nrow(y_keep_norm_disp))$table
stat_lfdHFD_LV <- topTags(LFD.vs.HFD, n=nrow(y_keep_norm_disp))$table
stat_txnHXN_LV <- topTags(TXN.vs.HXN, n=nrow(y_keep_norm_disp))$table

save(stat_txnHFD_LV, stat_txnLFD_LV, stat_hxnHFD_LV, 
	stat_lfdHFD_LV, stat_txnHXN_LV,
	file="DEGs_LV.rda")


is.de <- decideTestsDGE(TXN.vs.HFD)
is.de2 <- decideTestsDGE(HXN.vs.HFD)
is.de3 <- decideTestsDGE(LFD.vs.HFD)
is.de4 <- decideTestsDGE(TXN.vs.HXN)
is.de5 <- decideTestsDGE(TXN.vs.LFD)
#summary(decideTests(is.de))


p1 %<a-% {
	plotBCV(y_keep_norm_disp)
}

p2 %<a-% {
	plotQLDisp(fit)
}

p3 %<a-% {
	plotMD(TXN.vs.HFD, status=is.de, values=c(1,-1), col=c("red","green"),
	legend="topright")
}

p4 %<a-% {
	plotMD(HXN.vs.HFD, status=is.de2, values=c(1,-1), col=c("red","green"),
	legend="topright")
}

p5 %<a-% {
	plotMD(LFD.vs.HFD, status=is.de3, values=c(1,-1), col=c("red","green"),
	legend="topright")
}


p6 %<a-% {
	plotMD(TXN.vs.LFD, status=is.de5, values=c(1,-1), col=c("red","green"),
	legend="topright")
}



pdf("DEGs_LV.pdf", width=12, height=15)
p.all %<a-% {
	split.screen(c(3, 2))
	screen(1)
	p1
	screen(2)
	p2
	screen(3)
	p3
	screen(4)
	p4
	screen(5)
	p5
	screen(6)
	p6
	close.screen(all=TRUE)
}
p.all
dev.off()

setEPS()
postscript("DEGs_LV.eps", width=12, height=15)
p.all %<a-% {
	split.screen(c(3, 2))
	screen(1)
	p1
	screen(2)
	p2
	screen(3)
	p3
	screen(4)
	p4
	screen(5)
	p5
	screen(6)
	p6
	close.screen(all=TRUE)
}
p.all
dev.off()

print("writing to excel spreadsheet, this may take a while...")
write.xlsx(stat_lfdHFD_LV, file="DEGs_LV.xlsx",
	sheetName="LFDvs.HFD", append=FALSE)
write.xlsx(stat_txnHFD_LV, file="DEGs_LV.xlsx",
	sheetName="TXNvs.HFD", append=TRUE)
write.xlsx(stat_txnLFD_LV, file="DEGs_LV.xlsx",
	sheetName="TXNvs.LFD", append=TRUE)
write.xlsx(stat_hxnHFD_LV, file="DEGs_LV.xlsx",
	sheetName="HXNvs.HFD", append=TRUE)
write.xlsx(stat_txnHXN_LV, file="DEGs_LV.xlsx",
	sheetName="TXNvs.HXN", append=TRUE)

sessionInfo()