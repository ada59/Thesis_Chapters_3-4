################################################################################
# Script to test for metacom-level regulation (ADF test)
# Date: January 2024
# AFE
################################################################################


# Libraries: ===================================================================
library(tseries)
library(dplyr)
library(ggplot2)

# Main source: =================================================================
# https://www.r-bloggers.com/2022/06/augmented-dickey-fuller-test-in-r/
# https://rpubs.com/JSHAH/481706

# Data: ========================================================================
load("out.rawPA.RData")
load("out.rawPA_d.RData")
load("out.rawPA_u.RData")

load("out.rawAgg.RData")
load("out.rawAgg_d.RData")
load("out.rawAgg_u.RData")

# ExtractData: =================================================================

out.rawPA <- lapply(out.rawPA, function(x) {x$TDiNextEst$size_based})
out.rawPA0 <- lapply(out.rawPA, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawPA2 <- lapply(out.rawPA, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

out.rawPA_d <- lapply(out.rawPA_d, function(x) {x$TDiNextEst$size_based})
out.rawPA_d0 <- lapply(out.rawPA_d, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawPA_d2 <- lapply(out.rawPA_d, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

out.rawPA_u <- lapply(out.rawPA_u, function(x) {x$TDiNextEst$size_based})
out.rawPA_u0 <- lapply(out.rawPA_u, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawPA_u2 <- lapply(out.rawPA_u, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

out.rawAgg <- lapply(out.rawAgg, function(x) {x$TDiNextEst$size_based})
out.rawAgg0 <- lapply(out.rawAgg, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawAgg2 <- lapply(out.rawAgg, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

out.rawAgg_d <- lapply(out.rawAgg_d, function(x) {x$TDiNextEst$size_based})
out.rawAgg_d0 <- lapply(out.rawAgg_d, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawAgg_d2 <- lapply(out.rawAgg_d, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

out.rawAgg_u <- lapply(out.rawAgg_u, function(x) {x$TDiNextEst$size_based})
out.rawAgg_u0 <- lapply(out.rawAgg_u, function(x) {x[(x$Method=="Observed" & x$Order.q=="0"),]})
out.rawAgg_u2 <- lapply(out.rawAgg_u, function(x) {x[(x$Method=="Observed" & x$Order.q=="2"),]})

# Test: ========================================================================
# H0: The time series is non-stationary.
# HA: The time series is stationary.

# seems function adf.test assumes time points in order without need for the var


# PA ###########################################################################
# Northern Range (PA)-----------------------------------------------------------
adf.test(out.rawPA0$`Fish&Crustaceans`$qTD)   # p-value = 0.6731
adf.test(out.rawPA0$BenthicInvertebrates$qTD) # p-value = 0.99
adf.test(out.rawPA0$Diatoms$qTD)              # p-value = 0.01 (stat)
adf.test(out.rawPA0$Fish$qTD)                 # p-value = 0.7531

adf.test(out.rawPA2$`Fish&Crustaceans`$qTD)   # p-value = 0.9082
adf.test(out.rawPA2$BenthicInvertebrates$qTD) # p-value = 0.8106
adf.test(out.rawPA2$Diatoms$qTD)              # p-value = 0.01 (stat)
adf.test(out.rawPA2$Fish$qTD)                 # p-value = 0.6852


# Disturbed (PA)----------------------------------------------------------------
adf.test(out.rawPA_d0$`Fish&Crustaceans`$qTD)   # p-value = 0.7362
adf.test(out.rawPA_d0$BenthicInvertebrates$qTD) # p-value = 0.2041
adf.test(out.rawPA_d0$Diatoms$qTD)              # p-value = 0.6845
adf.test(out.rawPA_d0$Fish$qTD)                 # p-value = 0.6812

adf.test(out.rawPA_d2$`Fish&Crustaceans`$qTD)   # p-value = 0.8159
adf.test(out.rawPA_d2$BenthicInvertebrates$qTD) # p-value = 0.7852
adf.test(out.rawPA_d2$Diatoms$qTD)              # p-value = 0.5423
adf.test(out.rawPA_d2$Fish$qTD)                 # p-value = 0.6772


# Undisturbed (PA)--------------------------------------------------------------
adf.test(out.rawPA_u0$`Fish&Crustaceans`$qTD)   # p-value = 0.4969
adf.test(out.rawPA_u0$BenthicInvertebrates$qTD) # p-value = 0.5425
adf.test(out.rawPA_u0$Diatoms$qTD)              # p-value = 0.01 (stat)
adf.test(out.rawPA_u0$Fish$qTD)                 # p-value = 0.4821

adf.test(out.rawPA_u2$`Fish&Crustaceans`$qTD)   # p-value = 0.6377
adf.test(out.rawPA_u2$BenthicInvertebrates$qTD) # p-value = 0.6206
adf.test(out.rawPA_u2$Diatoms$qTD)              # p-value = 0.01 (stat)
adf.test(out.rawPA_u2$Fish$qTD)                 # p-value = 0.7046


# Aggregated ###################################################################
# Northern Range (agg)-----------------------------------------------------------

adf.test(out.rawAgg2$`Fish&Crustaceans`$qTD)   # p-value = 0.09752
adf.test(out.rawAgg2$BenthicInvertebrates$qTD) # p-value = 0.07252
adf.test(out.rawAgg2$Diatoms$qTD)              # p-value = 0.3588
adf.test(out.rawAgg2$Fish$qTD)                 # p-value = 0.07507


# Disturbed (agg)----------------------------------------------------------------

adf.test(out.rawAgg_d2$`Fish&Crustaceans`$qTD)   # p-value = 0.08625
adf.test(out.rawAgg_d2$BenthicInvertebrates$qTD) # p-value = 0.3145
adf.test(out.rawAgg_d2$Diatoms$qTD)              # p-value = 0.2926
adf.test(out.rawAgg_d2$Fish$qTD)                 # p-value = 0.08078


# Undisturbed (agg)--------------------------------------------------------------

adf.test(out.rawAgg_u2$`Fish&Crustaceans`$qTD)   # p-value = 0.4748
adf.test(out.rawAgg_u2$BenthicInvertebrates$qTD) # p-value = 0.02485 (stat)
adf.test(out.rawAgg_u2$Diatoms$qTD)              # p-value = 0.2342
adf.test(out.rawAgg_u2$Fish$qTD)                 # p-value = 0.4899

# End of script#################################################################
################################################################################