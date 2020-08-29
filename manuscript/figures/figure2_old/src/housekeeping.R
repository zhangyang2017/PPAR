library(knitr)
.cran_packages <- c("ggplot2", "dplyr", "tidyverse")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
	install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)


feeding <- read.csv("/Users/yangzhang/Box/MAN/pub/data/phenome_feeding.csv", as.is = TRUE)
feeding <- feeding %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))
feeding$Date <- as.Date(feeding$Date, "%m/%d/%y")
feeding$Treatment <- as.factor(feeding$Treatment)
feeding$Week <- as.factor(feeding$Week)
feeding$ID <- ordered(feeding$ID)
feeding$ID <- as.factor(feeding$ID)

weight_int <- filter(feeding, Date == "2019-02-04")
total_gained <- filter(feeding, Date == "2019-05-28")
total_food <- feeding %>% group_by(ID, Treatment) %>% summarise(total = sum(Food_Intake_weekly)) 
total_food <- total_food %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))
total_cal <- feeding %>% group_by(ID, Treatment) %>% summarise(total = sum(Calories_weekly))
total_cal <- total_cal %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))

## measurements
dexa <- read.csv("/Users/yangzhang/Box/MAN/pub/data/phenome_dexa_update.csv", as.is = TRUE)
dexa <- dexa %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))


## metabolic cage data
profile <- read.csv("/Users/yangzhang/Box/MAN/pub/data/phenome_metabolic_assessment.csv", as.is = TRUE)
profile$Treatment <- as.factor(profile$Treatment)
profile$MouseID <- ordered(profile$MouseID)
profile$MouseID <- as.factor(profile$MouseID)
profile <- profile %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))


bars <- read.csv("/Users/yangzhang/Box/MAN/pub/data/phenome_bars.csv", as.is = TRUE)
bars$Treatment <- as.factor(bars$Treatment)
bars$Cycle <- as.factor(bars$Cycle)
bars$MouseID <- ordered(bars$MouseID)
bars$MouseID <- as.factor(bars$MouseID)

bars <- bars %>%
	dplyr::mutate(vars = factor(Treatment,
								levels = c("HFD", "HFD+LXN", "HFD+HXN",
										   "HFD+TXN", "LFD")))


save(feeding, total_gained, total_cal, total_food, dexa, profile, bars, file = "figure1.rda")