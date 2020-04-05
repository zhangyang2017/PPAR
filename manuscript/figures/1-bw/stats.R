sink("analysis_output.txt")

#########################################################
library(knitr)
.cran_packages <- c("dplyr", "tidyverse", "emmeans")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
	install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
#########################################################

load("figure1.rda")
#sink("analysis_output.txt")
set.seed(12345)

Contrasts <- list(
	'LXN vs HFD' = c(-1, 1, 0, 0, 0),
	'HXN vs HFD' = c(-1, 0, 1, 0, 0),
	'TXN vs HFD' = c(-1, 0, 0, 1, 0),
	'LFD vs HFD' = c(-1, 0, 0, 0, 1),
	'TXN vs LFD' = c(0, 0, 0, 1, -1))

names <- colnames(dexa)
print(names)

cat("\n===================================================================================\n")
cat("Some summary statistics and Contrast comparison between HFD and other treatments:\n")
cat("===================================================================================\n")

stats_func = function(response){
	sums <- dexa %>%
		group_by(vars) %>%
		summarise_each(funs(mean = mean(., na.rm = TRUE), 
							median = median(., na.rm = TRUE), 
							SD = sd(., na.rm = TRUE), 
							SEM = sd(., na.rm = TRUE)/sqrt(n())), response)
	cat("Summary statistics", response, ":\n")
	print(sums)
	
	form = paste(response, "~ vars")
	model <- lm(as.formula(form), data = dexa)
	leastsquare <- lsmeans(model, "vars")
	output <- contrast(leastsquare, Contrasts, adjust = "none")
	cat("\nDifferences in", response, ":\n")
	print(output)
	cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
}

variables = names(dexa)[3:41]
variables %>%
	walk(stats_func)

sink()