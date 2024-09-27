################################################################################
# Script to test for com & metacom-level regulation (ADF test)
# Date: January 2024
# AFE
################################################################################


################################################################################
# LIBRARIES: ===================================================================
library(tseries)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(data.table)
library(patchwork)
library(tidyr)


# SOURCES: =====================================================================
# https://www.r-bloggers.com/2022/06/augmented-dickey-fuller-test-in-r/
# https://rpubs.com/JSHAH/481706

rm(list=ls())

################################################################################
# DATA: ========================================================================
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
load(paste0(Results_AlphaDivRData, "/S3_lPA_Alpha.RData"))

lPA_Alpha <- lPA_Alpha[!(lPA_Alpha$Order.q == 1),]
lPA_Alpha <- lPA_Alpha[!(lPA_Alpha$Session %in% c(23, 46, 50, 54)),]
lPA_Alpha <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Taxa, lPA_Alpha$Group, lPA_Alpha$SC2)) 
lPA_Alpha <- lPA_Alpha[which(lapply(lPA_Alpha, nrow) != 0)]     


################################################################################
# TESTS: =======================================================================
# H0: The time series is non-stationary.
# HA: The time series is stationary.
# seems function adf.test assumes time points in order without need for the var
# Create a list of dataframes


resMR <- lapply(lPA_Alpha, function(x) {adf.test(x$Estimate)}) # warnings about pvals OK

Code <- names(resMR)
restodf <- function(r){
    l <- bind_rows(map_dfr(r, ~data.frame(
    Statistic = .x$statistic,
    Pval = .x$p.value)))
    l$Code <- Code
    return(l)
}

resMR <- restodf(resMR)
resMR <- separate(resMR, Code, into = paste0("Column", 1:4), sep = "\\.")


resMR$Sig <- ifelse(resMR$Pval > 0.05, 0, 1)
resMR$TaxaOrderQ <- paste0(resMR$Column2, " ", "q", resMR$Column1)
resMR$TaxaOrderQ
levs <- c( "Fish (Trait) q2","Fish (Trait) q0",
           "Fish q2","Fish q0",
           "Invertebrates q2","Invertebrates q0",
           "Diatoms q2","Diatoms q0")
labs <- levs
resMR$TaxaOrderQ <- factor(resMR$TaxaOrderQ, levels = levs,
                       labels = labs)

S8_ADF_Results <- resMR[, !names(resMR) %in% c("Column1", "Column2", "Sig")]
rownames(S8_ADF_Results) <- 1:nrow(S8_ADF_Results)
names(S8_ADF_Results)[names(S8_ADF_Results)=="Column4"] <- "CoverageValue"
names(S8_ADF_Results)[names(S8_ADF_Results)=="Column3"] <- "Group"
S8_ADF_Results$Statistic <- round(S8_ADF_Results$Statistic, 2)
S8_ADF_Results$Pval <- round(S8_ADF_Results$Pval, 3)
rownames(S8_ADF_Results) <- 1:nrow(S8_ADF_Results)
S8_ADF_Results <- S8_ADF_Results %>% relocate(c(Group, CoverageValue, TaxaOrderQ), .before=Statistic)
unique(S8_ADF_Results$TaxaOrderQ)
orderTaxaQ <- c("Diatoms q0", "Diatoms q2",
                "Invertebrates q0",
                "Invertebrates q2",
                "Fish q0", "Fish q2", 
                "Fish (Trait) q0", "Fish (Trait) q2")


S8_ADF_Results_C5 <- subset(S8_ADF_Results, S8_ADF_Results$Group=="Northern Range")
S8_ADF_Results_C5$TaxaOrderQ <- factor(S8_ADF_Results_C5$TaxaOrderQ, levels = orderTaxaQ)
S8_ADF_Results_C5 <- S8_ADF_Results_C5 %>%
  arrange(desc(CoverageValue), TaxaOrderQ)

S8_ADF_Results_C6 <- subset(S8_ADF_Results, !S8_ADF_Results$Group=="Northern Range")
S8_ADF_Results_C6$TaxaOrderQ <- factor(S8_ADF_Results_C6$TaxaOrderQ, levels = orderTaxaQ)
S8_ADF_Results_C6 <- S8_ADF_Results_C6 %>%
  arrange(Group, desc(CoverageValue), TaxaOrderQ)

table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"
write.csv(S8_ADF_Results_C5, file=paste0(table_path, "/Chapter5/S8_C5_ADF_Results.csv"), row.names = F)
write.csv(S8_ADF_Results_C6, file=paste0(table_path, "/Chapter6/S8_C6_ADF_Results.csv"), row.names = F)



################################################################################
# PLOTS: =======================================================================

resMR_C5 <- subset(resMR, resMR$Column3 =="Northern Range" & resMR$Column4=="Cmin")     # or Cmax
resMR_C6_Dis <- subset(resMR, resMR$Column3 =="Disturbed" & resMR$Column4=="Cmin")      # or Cmax
resMR_C6_Undis <- subset(resMR, resMR$Column3 =="Undisturbed" & resMR$Column4=="Cmin")  # or Cmax

plot_reg <- function(x, main=NULL){
plot <- ggplot(x, aes(x = Statistic, y = TaxaOrderQ, color=as.factor(Sig))) +
     geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
     labs(title = main,
          x = "ADF Statistic",
          y = "") +
     xlim(-6, 0.5) +
     geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
     theme_minimal() +
     theme(
       plot.title = element_text(face = "bold", size = 14),
       axis.text.y = element_text(size = 12),
       axis.text.x = element_text(size = 12),
       axis.title.x = element_text(size = 12)
     )+
     scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")
return(plot)
}

plot_C5 <- plot_reg(resMR_C5, main = "Series Stationarity (Northern Range)")
dis_C6 <- plot_reg(resMR_C6_Dis, main = "Series Stationarity (Disturbed)")
undis_C6 <- plot_reg(resMR_C6_Undis, main = "Series Stationarity (Undisturbed)")
plot_C6 <- dis_C6 + undis_C6


################################################################################
# SAVE =========================================================================
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
ggsave(paste0(plot_path, "/Chapter5/S8_C5_ForestPlot.jpg"), plot_C5, width = 6, height = 6)
ggsave(paste0(plot_path, "/Chapter6/S8_C6_ForestPlot.jpg"), plot_C6, width = 12, height = 6)


# End of script#################################################################
################################################################################
