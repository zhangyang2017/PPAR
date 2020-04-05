library(knitr)
.cran_packages <- c("ggpubr", "gridExtra", "export")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
	install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(.cran_packages, require, character.only = TRUE)


load("figure1_par1.rda")
load("figure1_par2.rda")

par1 <- ggarrange(A, B, align = "hv", 
				  widths = c(2.5,1),
				  labels = c("A", "B"),
				  font.label = list(size = 18, color = "black"),
				  nrow = 1, ncol = 2)
par2 <- ggarrange(C, D, align = "hv", 
				  widths = c(2.5,1),
				  labels = c("C", "D"),
				  font.label = list(size = 18, color = "black"),
				  nrow = 1, ncol = 2,
				  common.legend = TRUE, legend="none")
#par3 <- ggarrange(E, Ff, align = "hv", 
#				  widths = c(2.5,1),
#				  labels = c("E", "F"),
#				  font.label = list(size = 22, color = "black"),
#				  nrow = 1, ncol = 2,
#				  common.legend = TRUE, legend="none")

partA <- ggarrange(par1, par2, align = "hv",
				   heights = c(1.5,1.5),
				   ncol = 1, nrow = 2,
				   common.legend = TRUE, legend="bottom")

partB <- ggarrange(E, F, G, H, align = "hv",
				   labels = c("E", "F", "G", "H"),
				   font.label = list(size = 18, color = "black"),
				   ncol = 4, nrow = 1, #heights = c(2,2),
				   common.legend = TRUE, legend="bottom")

final <- ggarrange(partA, partB,
		  nrow = 2, ncol=1,
		  heights = c(2, 0.8),
		  common.legend = TRUE, legend="bottom")

graph2svg(final, file = "figure1.svg", width = 12, height = 10)
graph2pdf(final, file = "figure1.pdf", width = 12, height = 10)
graph2eps(final, file = "figure1.eps", width = 12, height = 10)