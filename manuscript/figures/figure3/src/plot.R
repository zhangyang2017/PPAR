suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(export))
#####################################################

mycolor_alpha <- c("#8AB5F9", "#F8CB89", "#F08581", "#87A57E", "#797979")
theme_legend_free = function(){
	theme(axis.title.y = element_text(size = 17, 
	margin = margin(t=0, r=5, b=0, l=10)),
	axis.text.y = element_text(size = 17, face = 'bold'),
	axis.text.x = element_text(size = 17, face= 'bold'),
	axis.ticks.x = element_blank(),
	legend.position = "none",
    plot.title = element_text(size = 24, face = 'bold'),
	plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
	panel.grid = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank())
}

theme_legend_free2 = function(){
	theme(axis.title.y = element_text(size = 17, 
	margin = margin(t=0, r=5, b=0, l=10)),
	axis.text.y = element_text(size = 13),
	axis.text.x = element_text(size = 6, angle=45, hjust=1, face= 'bold'),
	axis.ticks.x = element_blank(),
	legend.position = "none",
	plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
	panel.grid = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank())
}

theme_legend = function(){
	theme(axis.title.y = element_text(size = 15, 
	margin = margin(t=0, r=5, b=0, l=10)),
	axis.text.y = element_text(size = 13),
	axis.text.x = element_blank(),
	axis.ticks.x = element_blank(),
	legend.title = element_text(colour = "white"),
	legend.text = element_text(size = 16),
	legend.position = "none",
	plot.margin = unit(c(1,0.5,0,0.5), "cm"),
	panel.grid = element_blank())
}



data <- read.csv('/Users/yangzhang/Box/ppar/PPAR/manuscript/figures/figure3/data/phenome_dexa_update.csv')
data$Treatment <- as.factor(data$Treatment)
data$vars <- factor(data$Treatment, levels = c("HFD", "HFD+LXN", "HFD+HXN", "HFD+TXN", "LFD"))


lvmass <- ggplot(data = data, aes(x = vars, y = lv_pct, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 4, width = 0.2, height = 0.1) +
	xlab(" ") + ylab("%") +
	ggtitle("Liver Mass Relative to Body Weight") +
    geom_hline(yintercept=4, linetype="dashed", color = "black")+
	theme_bw() +
	theme_legend_free()

label.df <- data.frame(vars = c("HFD+HXN", "HFD+TXN", "LFD"),
                       lv_pct = c(5.6, 3.7, 4.4))
lvmass <- lvmass + geom_text(data = label.df, label = c("**", "***", "***"), color='red', size = 12)


lvtag <- ggplot(data = data, aes(x = vars, y = LV_TAG, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 4, width = 0.2, height = 0.1) +
	xlab(" ") + ylab("mg/g") +
	ggtitle("Liver Triglyceride") +
	theme_bw() +
	theme_legend_free()

label.df <- data.frame(vars = c("HFD+TXN", "LFD"),
                       LV_TAG = c(625, 520))
lvtag <- lvtag + geom_text(data = label.df, label = c("**", "**"), color='red', size = 12)


plaTAG <- ggplot(data = data, aes(x = vars, y = total_TAG, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 4, width = 0.2, height = 0.1) +
	xlab(" ") + ylab("mg/dL") +
	ggtitle("Plasma Triglyceride") +
	theme_bw() +
	theme_legend_free()

label.df <- data.frame(vars = c("HFD+TXN"),
                       total_TAG = c(42))
plaTAG <- plaTAG + geom_text(data = label.df, label = c("**"), color='red', size = 12)
plaTAG

graph2pdf(lvmass, file = "lvmass.pdf", width = 8, height = 7)
graph2pdf(plaTAG, file = "plaTAG.pdf", width = 8, height = 7)
graph2pdf(lvtag, file = "lvtag.pdf", width = 8, height = 7)
