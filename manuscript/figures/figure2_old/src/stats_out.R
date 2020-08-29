#########################################################
library(knitr)
.cran_packages <- c("dplyr", "tidyverse", "lsmeans")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
	install.packages(.cran_packages[!.inst])
}
sapply(.cran_packages, require, character.only = TRUE)
#########################################################

load("figure1.rda")
sink("analysis_output.txt")
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

col5 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_BMD = sprintf("%0.2f", mean(BMD)), 
			  median_BMD = sprintf("%0.2f", median(BMD)),
			  se_BMD = sprintf("%0.2f", sd(BMD)/sqrt(n())))
cat("\n1.DEXA-BMD:\n")
print(col5)
model_BMD <- lm(BMD ~ vars, data = dexa)
leastsquare_BMD <- lsmeans(model_BMD, "vars")
output_BMD <- contrast(leastsquare_BMD, Contrasts, adjust = "none")
cat("\n1.Differences in DEXA-BMD:\n")
print(output_BMD)

col6 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_BMC = sprintf("%0.2f", mean(BMC)), 
			  median_BMC = sprintf("%0.2f", median(BMC)),
			  se_BMC = sprintf("%0.2f", sd(BMC)/sqrt(n())))
cat("\n2.DEXA-BMC:\n")
print(col6)
model_BMC <- lm(BMC ~ vars, data = dexa)
leastsquare_BMC <- lsmeans(model_BMC, "vars")
output_BMC <- contrast(leastsquare_BMC, Contrasts, adjust = "none")
cat("\n2.Differences in DEXA-BMC:\n")
print(output_BMC)

col7 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_BoneArea = sprintf("%0.2f", mean(Bone_Area)), 
			  median_BoneArea = sprintf("%0.2f", median(Bone_Area)),
			  se_BoneArea = sprintf("%0.2f", sd(Bone_Area)/sqrt(n())))
cat("\n3.DEXA-Bone area:\n")
print(col7)
model_BA <- lm(Bone_Area ~ vars, data = dexa)
leastsquare_BA <- lsmeans(model_BA, "vars")
output_BA <- contrast(leastsquare_BA, Contrasts, adjust = "none")
cat("\n3.Differences in DEXA-Bone area:\n")
print(output_BA)

col8 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_pctFat = sprintf("%0.2f", mean(Fat_percent)), 
			  median_pctFat = sprintf("%0.2f", median(Fat_percent)),
			  se_pctFat = sprintf("%0.2f", sd(Fat_percent)/sqrt(n())))
cat("\n4.DEXA-percent fat:\n")
print(col8)
model_pctFAT <- lm(Fat_percent ~ vars, data = dexa)
leastsquare_pctFAT <- lsmeans(model_pctFAT, "vars")
output_pctFAT <- contrast(leastsquare_pctFAT, Contrasts, adjust = "none")
cat("\n4.Differences in DEXA-percent fat:\n")
print(output_pctFAT)

col9 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_TTM = sprintf("%0.2f", mean(TTM)), 
			  median_TTM = sprintf("%0.2f", median(TTM)),
			  se_TTM = sprintf("%0.2f", sd(TTM)/sqrt(n())))
cat("\n5.DEXA-Total tissue mass:\n")
print(col9)
model_TTM <- lm(TTM ~ vars, data = dexa)
leastsquare_TTM <- lsmeans(model_TTM, "vars")
output_TTM <- contrast(leastsquare_TTM, Contrasts, adjust = "none")
cat("\n5.Differences in DEXA-Total tissue mass:\n")
print(output_TTM)

col10 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_TAM = sprintf("%0.2f", mean(TAM)), 
			  median_TAM = sprintf("%0.2f", median(TAM)),
			  se_TAM = sprintf("%0.2f", sd(TAM)/sqrt(n())))
cat("\n6.DEXA-Total animal mass:\n")
print(col10)
model_TAM <- lm(TAM ~ vars, data = dexa)
leastsquare_TAM <- lsmeans(model_TAM, "vars")
output_TAM <- contrast(leastsquare_TAM, Contrasts, adjust = "none")
cat("\n6.Differences in DEXA-Total animal mass:\n")
print(output_TAM)

col11 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_Fat = sprintf("%0.2f", mean(Fat_mass)), 
			  median_Fat = sprintf("%0.2f", median(Fat_mass)),
			  se_Fat = sprintf("%0.2f", sd(Fat_mass)/sqrt(n())))
cat("\n7.DEXA-Total fat mass:\n")
print(col11)
model_Fat <- lm(Fat_mass ~ vars, data = dexa)
leastsquare_Fat <- lsmeans(model_Fat, "vars")
output_Fat <- contrast(leastsquare_Fat, Contrasts, adjust = "none")
cat("\n7.Differences in DEXA-Total fat mass:\n")
print(output_Fat)

col13 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_Lean = sprintf("%0.2f", mean(Lean_mass)), 
			  median_Lean = sprintf("%0.2f", median(Lean_mass)),
			  se_Lean = sprintf("%0.2f", sd(Lean_mass)/sqrt(n())))
cat("\n8.DEXA-Total lean mass:\n")
print(col13)
model_Lean <- lm(Lean_mass ~ vars, data = dexa)
leastsquare_Lean <- lsmeans(model_Lean, "vars")
output_Lean <- contrast(leastsquare_Lean, Contrasts, adjust = "none")
cat("\n8.Differences in DEXA-Total lean mass:\n")
print(output_Lean)
	
weight <- total_gained %>% 
	group_by(vars) %>% 
	summarise(mean_wt_g_pct = sprintf("%0.2f", mean(Weight_gain_percent)), 
			  median_wt_g_pct = sprintf("%0.2f", median(Weight_gain_percent)),
			  se_wt_g_pct = sprintf("%0.2f", sd(Weight_gain_percent)/sqrt(n())),
			  mean_wt = sprintf("%0.2f",mean(Weight)), 
			  median_wt = sprintf("%0.2f",median(Weight)),
			  se_wt = sprintf("%0.2f", sd(Weight)/sqrt(n()))) 
cat("\n9.weight:\n")
print(weight)
model_wg <- lm(Weight_gain_percent ~ vars, data = total_gained)
leastsquare_wg <- lsmeans(model_wg, "vars")
output_wg <- contrast(leastsquare_wg, Contrasts, adjust = "none")
cat("\n9.Differences in % weight gained:\n")
print(output_wg)

food <- total_food %>% 
	group_by(vars) %>% 
	summarise(mean = mean(total), 
			  median = median(total),
	          se = sd(total)/sqrt(n()))
cat("\n10.food intake:\n")
print(food)
model_fd <- lm(total ~ vars, data = total_food)
leastsquare_fd <- lsmeans(model_fd, "vars")
output_fd <- contrast(leastsquare_fd, Contrasts, adjust = "none")
cat("\n10.Differences in total food consumed:\n")
print(output_fd)

calories <- total_cal %>% 
	group_by(vars) %>% 
	summarise(mean = mean(total), 
			  median = median(total),
			  se = sd(total)/sqrt(n()))
cat("\n11.calorie intake:\n")
print(calories)
model_cal <- lm(total ~ vars, data = total_cal)
leastsquare_cal <- lsmeans(model_cal, "vars")
output_cal <- contrast(leastsquare_cal, Contrasts, adjust = "none")
cat("\n11.Differences in total calories consumed:\n")
print(output_cal)

col16 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_LVwt = sprintf("%0.2f", mean(liver_wt)), 
			  median_LVwt = sprintf("%0.2f", median(liver_wt)),
			  se_LVwt = sprintf("%0.2f", sd(liver_wt)/sqrt(n())))
cat("\n12.Liver weight:\n")
print(col16)
model_lv <- lm(liver_wt ~ vars, data = dexa)
leastsquare_lv <- lsmeans(model_lv, "vars")
output_lv <- contrast(leastsquare_lv, Contrasts, adjust = "none")
cat("\n12.Differences in liver weight:\n")
print(output_lv)

col19 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_BATwt = sprintf("%0.2f", mean(bat_wt)), 
			  median_BATwt = sprintf("%0.2f", median(bat_wt)),
			  se_BATwt = sprintf("%0.2f", sd(bat_wt)/sqrt(n())))
cat("\n13.Brown adipose tissue weight:\n")
print(col19)
model_bat <- lm(bat_wt ~ vars, data = dexa)
leastsquare_bat <- lsmeans(model_bat, "vars")
output_bat <- contrast(leastsquare_bat, Contrasts, adjust = "none")
cat("\n13.Differences in BAT weight:\n")
print(output_bat)

col21 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_subqwt = sprintf("%0.2f", mean(subq_wt)), 
			  median_subqwt = sprintf("%0.2f", median(subq_wt)),
			  se_subqwt = sprintf("%0.2f", sd(subq_wt)/sqrt(n())))
cat("\n14.Subcutaneous fat weight:\n")
print(col21)
model_subq <- lm(subq_wt ~ vars, data = dexa)
leastsquare_subq <- lsmeans(model_subq, "vars")
output_subq <- contrast(leastsquare_subq, Contrasts, adjust = "none")
cat("\n14.Differences in subcutaneous fat weight:\n")
print(output_subq)

col23 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_viswt = sprintf("%0.2f", mean(visceral_wt)), 
			  median_viswt = sprintf("%0.2f", median(visceral_wt)),
			  se_viswt = sprintf("%0.2f", sd(visceral_wt)/sqrt(n())))
cat("\n15.Epididymal fat weight:\n")
print(col23)
model_vis <- lm(visceral_wt ~ vars, data = dexa)
leastsquare_vis <- lsmeans(model_vis, "vars")
output_vis <- contrast(leastsquare_vis, Contrasts, adjust = "none")
cat("\n15.Differences in epididymal fat weight:\n")
print(output_vis)

col25 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_mestwt = sprintf("%0.2f", mean(mesenteric_wt)), 
			  median_mestwt = sprintf("%0.2f", median(mesenteric_wt)),
			  se_mestwt = sprintf("%0.2f", sd(mesenteric_wt)/sqrt(n())))
cat("\n16.Mesenteric fat weight:\n")
print(col25)
model_mes <- lm(mesenteric_wt ~ vars, data = dexa)
leastsquare_mes <- lsmeans(model_mes, "vars")
output_mes <- contrast(leastsquare_mes, Contrasts, adjust = "none")
cat("\n16.Differences in mesenteric fat weight:\n")
print(output_mes)

col27 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_fecal_TAG = sprintf("%0.2f", mean(fecal_fat)), 
			  median_fecal_TAG = sprintf("%0.2f", median(fecal_fat)),
			  se_fecal_TAG = sprintf("%0.2f", sd(fecal_fat)/sqrt(n())))
cat("\n17.Fecal TAG:\n")
print(col27)
model_ftag <- lm(fecal_fat ~ vars, data = dexa)
leastsquare_ftag <- lsmeans(model_ftag, "vars")
output_ftag <- contrast(leastsquare_ftag, Contrasts, adjust = "none")
cat("\n17.Differences in fecal TAG:\n")
print(output_ftag)

col28 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_fecal_outp = sprintf("%0.2f", mean(fecal_tot_output)), 
			  median_fecal_outp = sprintf("%0.2f", median(fecal_tot_output)),
			  se_fecal_outp = sprintf("%0.2f", sd(fecal_tot_output)/sqrt(n())))
cat("\n18.3-day fecal mass:\n")
print(col28)
model_fecal_outp <- lm(fecal_tot_output ~ vars, data = dexa)
leastsquare_fecal_outp <- lsmeans(model_fecal_outp, "vars")
output_fecal_outp <- contrast(leastsquare_fecal_outp, Contrasts, adjust = "none")
cat("\n18.Differences in 3-day fecal mass:\n")
print(output_fecal_outp)

col29 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plCHL_baseline = sprintf("%0.2f", mean(CHL_plasma_baseline)), 
			  median_plCHL_baseline = sprintf("%0.2f", median(CHL_plasma_baseline)),
			  se_plCHL_baseline = sprintf("%0.2f", sd(CHL_plasma_baseline)/sqrt(n())))
cat("\n19.Baseline plasma cholesterol:\n")
print(col29)
model_plCHL_baseline <- lm(CHL_plasma_baseline ~ vars, data = dexa)
leastsquare_plCHL_baseline <- lsmeans(model_plCHL_baseline, "vars")
output_plCHL_baseline <- contrast(leastsquare_plCHL_baseline, Contrasts, adjust = "none")
cat("\n19.Differences in baseline plasma cholesterol:\n")
print(output_plCHL_baseline)

col30 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plCHL_mid = sprintf("%0.2f", mean(CHL_plasma_mid)), 
			  median_plCHL_mid = sprintf("%0.2f", median(CHL_plasma_mid)),
			  se_plCHL_mid = sprintf("%0.2f", sd(CHL_plasma_mid)/sqrt(n())))
cat("\n20.Mid plasma cholesterol:\n")
print(col30)
model_plCHL_mid <- lm(CHL_plasma_mid ~ vars, data = dexa)
leastsquare_plCHL_mid <- lsmeans(model_plCHL_mid, "vars")
output_plCHL_mid <- contrast(leastsquare_plCHL_mid, Contrasts, adjust = "none")
cat("\n20.Differences in mid plasma cholesterol:\n")
print(output_plCHL_mid)

col31 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plCHL_end = sprintf("%0.2f", mean(CHL_plasma_end)), 
			  median_plCHL_end = sprintf("%0.2f", median(CHL_plasma_end)),
			  se_plCHL_end = sprintf("%0.2f", sd(CHL_plasma_end)/sqrt(n())))
cat("\n21.End plasma cholesterol:\n")
print(col31)
model_plCHL_end <- lm(CHL_plasma_end ~ vars, data = dexa)
leastsquare_plCHL_end <- lsmeans(model_plCHL_end, "vars")
output_plCHL_end <- contrast(leastsquare_plCHL_end, Contrasts, adjust = "none")
cat("\n21.Differences in end plasma cholesterol:\n")
print(output_plCHL_end)

col33 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_fastGluc = sprintf("%0.2f", mean(glucose_plasma)), 
			  median_fastGluc = sprintf("%0.2f", median(glucose_plasma)),
			  se_fastGluc = sprintf("%0.2f", sd(glucose_plasma)/sqrt(n())))
cat("\n22.Fasting glucose:\n")
print(col33)
model_fastGluc <- lm(glucose_plasma ~ vars, data = dexa)
leastsquare_fastGluc <- lsmeans(model_fastGluc, "vars")
output_fastGluc <- contrast(leastsquare_fastGluc, Contrasts, adjust = "none")
cat("\n22.Differences in fasting glucose:\n")
print(output_fastGluc)

col34 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_fastIns = sprintf("%0.2f", mean(insulin_plasma)), 
			  median_fastIns = sprintf("%0.2f", median(insulin_plasma)),
			  se_fastIns = sprintf("%0.2f", sd(insulin_plasma)/sqrt(n())))
cat("\n23.Fasting insulin:\n")
print(col34)
model_fastIns <- lm(insulin_plasma ~ vars, data = dexa)
leastsquare_fastIns <- lsmeans(model_fastIns, "vars")
output_fastIns <- contrast(leastsquare_fastIns, Contrasts, adjust = "none")
cat("\n23.Differences in fasting insulin:\n")
print(output_fastIns)

col43 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_HOMA = sprintf("%0.2f", mean(HOMA)), 
			  median_HOMA = sprintf("%0.2f", median(HOMA)),
			  se_HOMA = sprintf("%0.2f", sd(HOMA)/sqrt(n())))
cat("\n24.HOMA-IR:\n")
print(col43)
model_HOMA <- lm(HOMA ~ vars, data = dexa)
leastsquare_HOMA <- lsmeans(model_HOMA, "vars")
output_HOMA <- contrast(leastsquare_HOMA, Contrasts, adjust = "none")
cat("\n24.Differences in HOMA-IR:\n")
print(output_HOMA)

col45 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_leptin = sprintf("%0.2f", mean(leptin_new)), 
			  median_leptin = sprintf("%0.2f", median(leptin_new)),
			  se_leptin = sprintf("%0.2f", sd(leptin_new)/sqrt(n())))
cat("\n25.Fasting leptin:\n")
print(col45)
model_leptin_new <- lm(leptin_new ~ vars, data = dexa)
leastsquare_leptin_new <- lsmeans(model_leptin_new, "vars")
output_leptin_new <- contrast(leastsquare_leptin_new, Contrasts, adjust = "none")
cat("\n25.Differences in fasting leptin:\n")
print(output_leptin_new)

col41 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plTAG = sprintf("%0.2f", mean(total_TAG)), 
			  median_plTAG = sprintf("%0.2f", median(total_TAG)),
			  se_plTAG = sprintf("%0.2f", sd(total_TAG)/sqrt(n())))
cat("\n26.Fasting circulating TAG:\n")
print(col41)
model_plTAG <- lm(total_TAG ~ vars, data = dexa)
leastsquare_plTAG <- lsmeans(model_plTAG, "vars")
output_plTAG <- contrast(leastsquare_plTAG, Contrasts, adjust = "none")
cat("\n26.Differences in fasting circulating TAG:\n")
print(output_plTAG)

col38 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plCHL = sprintf("%0.2f", mean(total_CHL)), 
			  median_plCHL = sprintf("%0.2f", median(total_CHL)),
			  se_plCHL = sprintf("%0.2f", sd(total_CHL)/sqrt(n())))
cat("\n27.Fasting circulating cholesterol:\n")
print(col38)
model_plCHL <- lm(total_CHL ~ vars, data = dexa)
leastsquare_plCHL <- lsmeans(model_plCHL, "vars")
output_plCHL <- contrast(leastsquare_plCHL, Contrasts, adjust = "none")
cat("\n27.Differences in fasting circulating cholesterol:\n")
print(output_plCHL)

col39 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plLDL = sprintf("%0.2f", mean(LDL)), 
			  median_plLDL = sprintf("%0.2f", median(LDL)),
			  se_plLDL = sprintf("%0.2f", sd(LDL)/sqrt(n())))
cat("\n28.Fasting circulating LDL:\n")
print(col39)
model_plLDL <- lm(LDL ~ vars, data = dexa)
leastsquare_plLDL <- lsmeans(model_plLDL, "vars")
output_plLDL <- contrast(leastsquare_plLDL, Contrasts, adjust = "none")
cat("\n28.Differences in fasting circulating LDL:\n")
print(output_plLDL)

col40 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_plHDL = sprintf("%0.2f", mean(HDL)), 
			  median_plHDL = sprintf("%0.2f", median(HDL)),
			  se_plHDL = sprintf("%0.2f", sd(HDL)/sqrt(n())))
cat("\n29.Fasting circulating HDL:\n")
print(col40)
model_plHDL <- lm(HDL ~ vars, data = dexa)
leastsquare_plHDL <- lsmeans(model_plHDL, "vars")
output_plHDL <- contrast(leastsquare_plHDL, Contrasts, adjust = "none")
cat("\n29.Differences in fasting circulating HDL:\n")
print(output_plHDL)

index <- dexa %>%
	group_by(vars) %>%
	summarise(mean_index = sprintf("%0.2f", mean(LDL/HDL)), 
			  median_index = sprintf("%0.2f", median(LDL/HDL)),
			  se_index = sprintf("%0.2f", sd(LDL/HDL)/sqrt(n())))
cat("\n30.Atherogenic index:\n")
print(index)
model_index <- lm(LDL/HDL ~ vars, data = dexa)
leastsquare_index <- lsmeans(model_index, "vars")
output_index <- contrast(leastsquare_index, Contrasts, adjust = "none")
cat("\n30.Differences in atherogenic index:\n")
print(output_index)

col50 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_ketone = sprintf("%0.2f", mean(ketone)), 
			  median_ketone = sprintf("%0.2f", median(ketone)),
			  se_ketone = sprintf("%0.2f", sd(ketone)/sqrt(n())))
cat("\n31.Fasting circulating ketone:\n")
print(col50)
model_ketone <- lm(ketone ~ vars, data = dexa)
leastsquare_ketone <- lsmeans(model_ketone, "vars")
output_ketone <- contrast(leastsquare_ketone, Contrasts, adjust = "none")
cat("\n31.Differences in fasting circulating ketone:\n")
print(output_ketone)

col27 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_FecalTAG = sprintf("%0.2f", mean(fecal_fat)), 
			  median_FecalTAG = sprintf("%0.2f", median(fecal_fat)),
			  se_FecalTAG = sprintf("%0.2f", sd(fecal_fat)/sqrt(n())))
cat("\n32.Fecal TAG:\n")
print(col27)
model_FecalTAG <- lm(fecal_fat ~ vars, data = dexa)
leastsquare_FecalTAG <- lsmeans(model_FecalTAG, "vars")
output_FecalTAG <- contrast(leastsquare_FecalTAG, Contrasts, adjust = "none")
cat("\n32.Differences in fecal TAG:\n")
print(output_FecalTAG)

col49 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_FecalCHL = sprintf("%0.2f", mean(CHL_feces)), 
			  median_FecalCHL = sprintf("%0.2f", median(CHL_feces)),
			  se_FecalCHL = sprintf("%0.2f", sd(CHL_feces)/sqrt(n())))
cat("\n33.Fecal cholesterol:\n")
print(col49)
model_FecalCHL <- lm(CHL_feces ~ vars, data = dexa)
leastsquare_FecalCHL <- lsmeans(model_FecalCHL, "vars")
output_FecalCHL <- contrast(leastsquare_FecalCHL, Contrasts, adjust = "none")
cat("\n33.Differences in fecal cholesterol:\n")
print(output_FecalCHL)

col46 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_LV_TC = sprintf("%0.2f", mean(LV_TC)), 
			  median_LV_TC = sprintf("%0.2f", median(LV_TC)),
			  se_LV_TC = sprintf("%0.2f", sd(LV_TC)/sqrt(n())))
cat("\n34.Liver cholesterol:\n")
print(col46)
model_LV_TC <- lm(LV_TC ~ vars, data = dexa)
leastsquare_LV_TC <- lsmeans(model_LV_TC, "vars")
output_LV_TC <- contrast(leastsquare_LV_TC, Contrasts, adjust = "none")
cat("\n34.Differences in liver cholesterol:\n")
print(output_LV_TC)

col47 <- dexa %>%
	group_by(vars) %>%
	summarise(mean_LV_TAG = sprintf("%0.2f", mean(LV_TAG)), 
			  median_LV_TAG = sprintf("%0.2f", median(LV_TAG)),
			  se_LV_TAG = sprintf("%0.2f", sd(LV_TAG)/sqrt(n())))
cat("\n35.Liver TAG:\n")
print(col47)
model_LV_TAG <- lm(LV_TAG ~ vars, data = dexa)
leastsquare_LV_TAG <- lsmeans(model_LV_TAG, "vars")
output_LV_TAG <- contrast(leastsquare_LV_TAG, Contrasts, adjust = "none")
cat("\n35.Differences in liver TAG:\n")
print(output_LV_TAG)

ee <- profile %>% 
	group_by(vars) %>% 
	summarise(mean = mean(EE), 
			  median = median(EE),
			  se = sd(EE)/sqrt(n()))
cat("\n36.Hourly energy expenditure:\n")
print(ee)
model_ee <- lm(EE ~ vars, data = profile)
leastsquare_ee <- lsmeans(model_ee, "vars")
output_ee <- contrast(leastsquare_ee, Contrasts, adjust = "none")
cat("\n36.Differences in hourly energy expenditure:\n")
print(output_ee)

light <- bars %>% filter(Cycle == "1.light")
dark <- bars %>% filter(Cycle == "2.dark")
lit <- light %>% 
	group_by(vars) %>% 
	summarise(mean = mean(AllMeters), 
			  median = median(AllMeters),
			  se = sd(AllMeters)/sqrt(n()))
cat("\n37.Total meters traveled: light cycle:\n")
print(lit)
model_lit <- lm(AllMeters ~ vars, data = light)
leastsquare_lit <- lsmeans(model_lit, "vars")
output_lit <- contrast(leastsquare_lit, Contrasts, adjust = "none")
cat("\n37.Differences in total meters traveled: light cycle:\n")
print(output_lit)

dar <- dark %>% 
	group_by(vars) %>% 
	summarise(mean = mean(AllMeters), 
			  median = median(AllMeters),
			  se = sd(AllMeters)/sqrt(n()))
cat("\n38.Total meters traveled: dark cycle:\n")
print(dar)
model_dar <- lm(AllMeters ~ vars, data = dark)
leastsquare_dar <- lsmeans(model_dar, "vars")
output_dar <- contrast(leastsquare_dar, Contrasts, adjust = "none")
cat("\n38.Differences in total meters traveled:dark cycle:\n")
print(output_dar)

sink()