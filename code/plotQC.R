###############################################################################
library(knitr)
.cran_packages <- c("ggplot2", "emmeans", "ggpubr", "RColorBrewer", "edgeR", "pryr")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
###############################################################################

load("transcriptome_LV.rda")

#y y_keep
mycolor <- c("#999999", "lightblue", "#E69F00", "lawngreen")
plot_y <- y$samples %>%
	dplyr::mutate(vars = factor(Diet, levels = c("LFD", "HFD", "HXN", "TXN")))
rownames(plot_y) <- rownames(y$samples)

rownames(plot_y_keep) <- rownames(y_keep$samples)
plot_y_keep <- y_keep$samples %>%
	dplyr::mutate(vars = factor(Diet, levels = c("LFD", "HFD", "HXN", "TXN")))
rownames(plot_y_keep) <- rownames(y_keep$samples)

A <- ggplot(plot_y, aes(x=row.names(plot_y), y=lib.size/1000000, fill=vars))+
	geom_bar(stat = 'identity', color="black")+
	geom_hline(yintercept=mean((plot_y$lib.size)/1000000), col="black", lty = 5) +
	geom_hline(yintercept = median((plot_y$lib.size)/1000000), col="red", lty = 5) +
	ggtitle("SPF Liver Library Sizes: Pre-filter") +
	xlab("") +
	ylab("million")+
	scale_fill_manual(values=alpha(mycolor, .6))+
	theme_bw() +
	theme(axis.title.y = element_text(size = 15, margin = margin(t=0, r=10, b=0, l=10)),
		  axis.text.y = element_text(size = 12),
		  axis.text.x = element_text(size = 8, angle = 90, hjust=1),
		  legend.position = "right",
		  legend.title = element_text(colour = "white"),
		  plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
		  panel.grid = element_blank())

B <- ggplot(plot_y_keep, aes(x=row.names(plot_y_keep), y=lib.size/1000000, fill=vars))+
	geom_bar(stat = 'identity', color="black")+
	geom_hline(yintercept=mean((plot_y_keep$lib.size)/1000000), col="black", lty = 5) +
	geom_hline(yintercept = median((plot_y_keep$lib.size)/1000000), col="red", lty = 5) +
	ggtitle("SPF Liver Library Sizes: Post-filter") +
	xlab("") +
	ylab("million")+
	scale_fill_manual(values=alpha(mycolor, .6))+
	theme_bw() +
	theme(axis.title.y = element_text(size = 15, margin = margin(t=0, r=10, b=0, l=10)),
		  axis.text.y = element_text(size = 12),
		  axis.text.x = element_text(size = 8, angle = 90, hjust=1),
		  legend.position = "right",
		  legend.title = element_text(colour = "white"),
		  plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
		  panel.grid = element_blank())

par1 <- ggarrange(A, B, align = "hv", 
				  font.label = list(size = 20, color = "black"),
				  nrow = 1, ncol = 2,
				  common.legend = TRUE, legend="bottom")


cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")

p1.pryr %<a-% {
	plot(density(lcpm[,1]), col=col[1], lwd=2, 
		 ylim=c(0,0.26), las=1, main="", xlab="")
	title(main="LogCPMs: Pre-filter", xlab="logCPM")
	abline(v=lcpm.cutoff, lty=3)
	for (i in 2:nsamples){
		den <- density(lcpm[,i])
		lines(den$x, den$y, col=col[i], lwd=2)
	}
}

cpm_keep <- cpm(y_keep)
lcpm_keep <- cpm(y_keep, log=TRUE)

L_keep <- mean(y_keep$samples$lib.size) * 1e-6
M_keep <- median(y_keep$samples$lib.size) * 1e-6
lcpm.cutoff_keep <- log2(10/M_keep + 2/L_keep)


p2.pryr %<a-% {
	plot(density(lcpm_keep[,1]), col=col[1], lwd=2, 
		 ylim=c(0,0.26), las=1, main="", xlab="")
	title(main="LogCPMs: Post-filter", xlab="logCPM")
	abline(v=lcpm.cutoff_keep, lty=3)
	for (i in 2:nsamples){
		den <- density(lcpm_keep[,i])
		lines(den$x, den$y, col=col[i], lwd=2)
	}
}

pdf("blah_base_both_pryr.pdf", width=8, height=5)
p.both %<a-% {
	split.screen(c(1, 2))
	
	screen(1)
	p1.pryr
	
	screen(2)
	p2.pryr
	
	close.screen(all=TRUE)
}
p.both
dev.off()


mycolors2 <- c(
	"#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999", "#999999",
	"lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", 
	"#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00", "#E69F00",  
	"lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen", "lawngreen")

logcounts_y <- cpm(y, log=TRUE)


ggplot(logcounts_y, aes(x = row.names(plot_y), y = logcounts_y)) +
	geom_boxplot()

boxplot(logcounts_y, xlab="", ylab="Log2 counts per million", las=2,
    font.axis = 1,
    col = mycolors2)
abline(h=median(logcounts), col="blue")
title("LogCPMs: Pre-filter")


pdf("logcpm_dist_raw.pdf", width=8, height=5)
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)
#summary(lcpm)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")

plot(density(lcpm[,1]), col=col[1], lwd=2, 
    ylim=c(0,0.26), las=1, main="", xlab="")
title(main="LogCPMs: Pre-filter", xlab="logCPM")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()

table(rowSums(y$counts==0)==11)






pdf("library.size.filter.pdf", width=8, height=5)
par(mfrow=c(1,1), mai=c(0.7,1,0.5,0.3), #bottom, right, top, left
    cex.main=1.5, cex.lab=1.5, cex.axis=0.5)

barplot((y_keep$samples$lib.size)/1000000, 
    ylim=c(0, 9),
    names=colnames(y_keep), las=2, font.axis = 1,
    xlab = "", ylab = "million",
    col = mycolors2)
abline(h = mean((y_keep$samples$lib.size)/1000000), col="black", lty = 5)
abline(h = median((y_keep$samples$lib.size)/1000000), col="red", lty = 5)
title("Library Sizes: Post-filter")
dev.off()



pdf("logcpm_filter.pdf", width=8, height=5)
par(mfrow=c(1,1), mai=c(0.7,1,0.5,0.3), #bottom, right, top, left
    cex.main=1.5, cex.lab=1.5, cex.axis=0.5)
logcounts <- cpm(y.keep, log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2,
    font.axis = 1,
    col = mycolors2)
abline(h=median(logcounts), col="blue")
title("LogCPMs: Post-filter")
dev.off()


pdf("logcpm_dist_filter.pdf", width=8, height=5)
cpm <- cpm(y.keep)
lcpm <- cpm(y.keep, log=TRUE)
#summary(lcpm)
L <- mean(y.keep$samples$lib.size) * 1e-6
M <- median(y.keep$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y.keep)
col <- brewer.pal(nsamples, "Paired")

plot(density(lcpm[,1]), col=col[1], lwd=2, 
    ylim=c(0,0.26), las=1, main="", xlab="")
title(main="LogCPMs: Post-filter", xlab="logCPM")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()


pdf("MDS_raw.pdf")
y$samples$Diet <- relevel(y$samples$Diet, "LFD")
points <- c(0, 1, 2, 15)
colors <- c("black", "blue", "orange", "darkgreen")
plotMDS(y, top = 100, dim.plot = c(1,2), 
    col=colors[y$samples$Diet], 
    pch=points[y$samples$Diet], 
    gene.selection="pairwise", cex=1)
legend("topright", legend=levels(y$samples$Diet), pch=points, col=colors, ncol=1)
dev.off()







### post-filter
pdf("post-filter.pdf", width=30, height=9)
par(mfrow=c(1,3), mai=c(0.7,1,0.5,0.3), #bottom, right, top, left
    cex.main=2, cex.lab=2, cex.axis=1.5)

barplot((y.keep$samples$lib.size)/1000000, 
    ylim=c(0, 8),
    names=colnames(y.keep), las=2,
    font.axis = 1,
    xlab = "", ylab = "million",
    col = mycolors2)
abline(h = mean((y.keep$samples$lib.size)/1000000), col="black", lty = 5)
abline(h = median((y.keep$samples$lib.size)/1000000), col="red", lty = 5)
title("A. Library Sizes: Post-filter")


logcounts <- cpm(y.keep, log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2,
    font.axis = 1,
    col = mycolors2)
abline(h=median(logcounts), col="blue")
title("B. LogCPMs: Post-filter")

cpm <- cpm(y.keep)
lcpm <- cpm(y.keep, log=TRUE)
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y.keep)
col <- brewer.pal(nsamples, "Paired")

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), main="", xlab="")
title(main="C. LogCPMs: Post-filter", xlab="logCPM")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()


sessionInfo()
### the end



y.keep.norm <- calcNormFactors(y.keep, method = "TMM")
#y.keep.norm$samples$norm.factors


y.keep.norm2 <- y.keep.norm   #x2 <- x
y.keep.norm2$samples$norm.factors <- 1
y.keep.norm2$counts[,1] <- ceiling(y.keep.norm2$counts[,1]*0.05)
y.keep.norm2$counts[,2] <- y.keep.norm2$counts[,2]*5


pdf("normalization.pdf", width=18, height=9)
par(mfrow=c(1,2))
lcpm <- cpm(y.keep.norm2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")

y.keep.norm3 <- calcNormFactors(y.keep.norm2) 
lcpm <- cpm(y.keep.norm3, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")
dev.off()


pdf("normalization_subset.pdf", width=18, height=9)
par(mfrow=c(1,2))
nsamples <- ncol(y.keep)
col <- brewer.pal(nsamples, "Paired")
lcpm <- cpm(y.keep, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data (subset)",ylab="Log-cpm")

lcpm <- cpm(y.keep.norm, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data (subset)",ylab="Log-cpm")
dev.off()



lcpm <- cpm(y.keep.norm, log=TRUE)
pdf("normalization2.pdf", width=9, height=9)
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")


pdf("MDS_filter.pdf")
y.keep.norm$samples$Diet <- relevel(y.keep.norm$samples$Diet, "LFD")
points <- c(0, 1, 2, 15)
colors <- c("black", "blue", "orange", "darkgreen")
plotMDS(y.keep.norm, top = 100, dim.plot = c(1,2), 
    col=colors[y.keep.norm$samples$Diet], 
    pch=points[y.keep.norm$samples$Diet], 
    gene.selection="pairwise", cex=1)
legend("topright", legend=levels(y.keep.norm$samples$Diet), pch=points, col=colors, ncol=1)
dev.off()



design <- model.matrix(~ 0 + y.keep.norm$samples$Diet)
colnames(design) <- levels(y.keep.norm$samples$Diet)


#design <- model.matrix(~ 0 + y.norm.keep$samples$group)
#colnames(design) <- levels(y.norm.keep$samples$group)

## estimate dispersions
y.keep.norm.disp <- estimateDisp(y.keep.norm, design, robust=TRUE)


saveRDS(y.keep.norm.disp, file = "subset39.keep.norm.disp.rds")

design <- model.matrix(~ 0 + y.keep.norm.disp$samples$Diet)
colnames(design) <- levels(y.keep.norm.disp$samples$Diet)

#design <- model.matrix(~ 0 + y.norm.keep.disp$samples$group)
#colnames(design) <- levels(y.norm.keep.disp$samples$group)


## test for DE genes

#design <- model.matrix(~ 0 + y.norm.keep.disp$samples$kind)
#colnames(design) <- levels(y.norm.keep.disp$samples$kind)
### GF vs. conv
fit <- glmQLFit(y.keep.norm.disp, design, robust=TRUE)
con <- makeContrasts(TXN - HFD, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
summary(decideTests(qlf))
stat_txnHFD_qlf <- topTags(qlf, n=nrow(y.keep.norm.disp))$table

write.csv(stat_txnHFD_qlf, "stat_txnHFD_qlf.csv")


pdf("pvalues_TXN_vs_HFD_subset39.pdf")
theme_set(theme_bw())

alpha = binw = 0.025
pi0 = 2*mean(stat_txnHFD_qlf$PValue > 0.5)
ggplot(as(stat_txnHFD_qlf, "data.frame"), aes(x=PValue)) + 
geom_histogram(binwidth=0.01, fill="Royalblue", boundary = 0) +
geom_hline(yintercept = pi0 * binw * nrow(stat_txnHFD_qlf), col = "blue") +
geom_vline(xintercept = alpha, col = "red") +
ggtitle("Histogram of p-values: TXN vs. HFD")
dev.off()


pdf("libSize_RIN.pdf", width=5, height=5)
plot(y$samples$RIN, (y$samples$lib.size)/1000000, 
	 xlim=c(5, 9), ylim=c(2, 11),
	 main="Liver (conventional)",
	 xlab="RIN", ylab="Lib Size in Millions", pch=19)
dev.off()

pdf("libSize_RNA.pdf", width=5, height=5)
plot(y$samples$RNA_input, (y$samples$lib.size)/1000000, 
	 ylim=c(2, 11),
	 main="Liver (conventional)",
	 xlab="RNA Input (µl)", ylab="Lib Size in Millions", pch=19)
dev.off()