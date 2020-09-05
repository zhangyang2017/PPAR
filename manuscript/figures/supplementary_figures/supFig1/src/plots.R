#############################################################################################
library(knitr)
.cran_packages <- c("ggplot2", "emmeans", "ggpubr", 'dplyr')
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
	install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
#############################################################################################
#theme_set(theme_bw())

load("figure1.rda")
mycolor_solid <- c("#0080ff", "#ffa500", "red", "darkgreen", "black")
mycolor_alpha <- c("#8AB5F9", "#F8CB89", "#F08581", "#87A57E", "#797979")
lines <- c("solid", "solid", "solid", "solid", "dashed")
shapes <- c(24, 4, 15, 25, 1)

### for boxplot
pars <- function(x){
	r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}
outliers <- function(x){
	subset(x, x < quantile(x, probs = 0.1) | quantile(x, probs = 0.9) < x)
}

theme_legend_free = function(){
	theme(axis.title.y = element_text(size = 15, 
	margin = margin(t=0, r=5, b=0, l=10)),
	axis.text.y = element_text(size = 13),
	axis.text.x = element_blank(),
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


#############################################################################################

F1A <- ggline(feeding, x = "Week", y = "Weight_gain_percent", group = "vars",
			  add = "mean_se", color = "vars", palette = mycolor_solid,
			  ylab = "%", 
			  linetype = "vars", shape = "vars", 
			  point.size = 2,
			  ggtheme = theme_bw())

A <- F1A + scale_linetype_manual(values=c( "solid", "solid", "solid", "solid", "dashed")) + 
	scale_shape_manual(values = shapes) +
	ggtitle("Weekly Weight Gain (%)") +
	theme(axis.title.y = element_text(size = 15, margin = margin(t=0, r=10, b=0, l=10)),
		  axis.text.y = element_text(size = 15),
		  axis.text.x = element_text(size = 18),
		  legend.position = "top",
		  legend.title = element_blank(),
		  legend.margin = margin(c(0,0,0,0)),
		  plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
		  panel.grid = element_blank())

A <- A + theme(legend.position = c(0.1, 0.6), legend.background = element_rect())
#############################################################################################

B <- ggplot(data = total_gained, aes(x = vars, y = Weight_gain_percent, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab(" ") + ylab("%") +
	ggtitle("Cumulative Weight Gained") +
	theme_bw() +
	theme_legend_free()

label.df <- data.frame(vars = c("HFD+TXN", "LFD"),
                       Weight_gain_percent = c(25, 25))
B <- B + geom_text(data = label.df, label = c("***", "***"), color='red', size = 6)
#############################################################################################

food <- ggline(feeding, x = "Week", y = "Food_Intake_weekly", group = "vars",
			   add = "mean_se", color = "vars", palette = mycolor_solid,
			   ylab = "g", 
			   linetype = "vars", shape = "vars",
			   point.size = 2,
			   ggtheme = theme_bw())

C <- food + scale_shape_manual(values = shapes) +
	scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed")) + 
	ggtitle("Weekly Food Intake") +
	theme(axis.title.y = element_text(size = 15, margin = margin(t=0, r=10, b=0, l=10)),
		  axis.text.y = element_text(size = 15),
		  axis.text.x = element_text(size = 18),
		  legend.position = "none",
		  plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
		  panel.grid = element_blank())
#############################################################################################

D <- ggplot(data = total_cal, aes(x = vars, y = total, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("kcal") +
	ggtitle("Cumulative Calorie Intake") +
	theme_bw() +
	theme_legend_free()

label.df <- data.frame(vars = c("HFD+HXN", "LFD"),
                       total = c(1800, 1820))
D <- D + geom_text(data = label.df, label = c("*", "***"), color='red', size = 6)
#############################################################################################

E <- ggplot(data = dexa, aes(x = vars, y = Fat_mass, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("g") +
	ggtitle("Total Fat Mass") +
	theme_bw() +
	theme_legend()

label.df <- data.frame(vars = c("HFD+HXN", "HFD+TXN", "LFD"),
                       Fat_mass = c(12, 10, 15))
E <- E + geom_text(data = label.df, label = c("*", "***", "***"), color='red', size = 6)
#############################################################################################

F <- ggplot(data = dexa, aes(x = vars, y = mesenteric_wt, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("g") +
	ggtitle("Mesenteric Fat mass") +
	theme_bw() +
	theme_legend()

label.df <- data.frame(vars = c("HFD+HXN", "HFD+TXN", "LFD"),
                       mesenteric_wt = c(1.5, 1.2, 1.2))
F <- F + geom_text(data = label.df, label = c("*", "***", "***"), color='red', size = 6)
#############################################################################################

G <- ggplot(data = dexa, aes(x = vars, y = subq_wt, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("g") +
	ggtitle("Subcutaneous Fat mass") +
	theme_bw() +
	theme_legend()

label.df <- data.frame(vars = c("HFD+HXN", "HFD+TXN", "LFD"),
                       subq_wt = c(1, 1.3, 1.3))
G <- G + geom_text(data = label.df, label = c("**", "***", "***"), color='red', size = 6)
#############################################################################################

H <- ggplot(data = dexa, aes(x = vars, y = epi_wat, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("g") +
	ggtitle("Epididymal Fat mass") +
	theme_bw() +
	theme_legend()

label.df <- data.frame(vars = c("HFD+HXN", "LFD"),
                       epi_wat = c(2.4, 2.3))
H <- H + geom_text(data = label.df, label = c("***", "**"), color='red', size = 6)
#############################################################################################

J <- ggplot(data = dexa, aes(x = vars, y = TG_feces, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("mg/g feces") +
	ggtitle("Fecal Triglycerides") +
	theme_bw() +
	theme_legend()
#############################################################################################

K <- ggplot(data = dexa, aes(x = vars, y = total_TAG, fill = vars)) +
	stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = mycolor_alpha) +
	geom_jitter(shape = 20, size = 1, width = 0.2, height = 0.1) +
	xlab("") + ylab("mg/dL") +
	ggtitle("Fasting Plasma Triglycerides") +
	theme_bw() +
	theme_legend()

label.df <- data.frame(vars = c("HFD+TXN"),
                       total_TAG = c(65))
K <- K + geom_text(data = label.df, label = c("**"), color='red', size = 6)
#############################################################################################
#############################################################################################
save(A, B, C, D, file = "figure1_par1.rda")
save(E, F, G, H, I, J, K, file = "figure1_par2.rda")
#############################################################################################