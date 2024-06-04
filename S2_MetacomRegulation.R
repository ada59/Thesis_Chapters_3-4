################################################################################
# Script to test for com & metacom-level regulation (ADF test)
# Date: January 2024
# AFE
################################################################################


# Libraries: ===================================================================
library(tseries)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(data.table)
library(patchwork)


# Main source: =================================================================
# https://www.r-bloggers.com/2022/06/augmented-dickey-fuller-test-in-r/
# https://rpubs.com/JSHAH/481706


rm(list=ls())


# Data: ========================================================================
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
setwd(Results_AlphaDivRData)

#load("out.rawPA_Alpha_Taxonomic.RData")
#load("out.rawAgg_Alpha_Taxonomic.RData")
#load("out.rawAbu_Alpha_Taxonomic.RData")

#load("out.rawPA_Alpha_FishTrait.RData")
#load("out.rawAgg_Alpha_FishTrait.RData")
load("lPA_Alpha.RData")
setwd("C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs")


# ExtractData: =================================================================
lPA_Alpha$Session <- as.integer(lPA_Alpha$Session)
names(lPA_Alpha)[names(lPA_Alpha)=="Alpha"] <- "Estimate"
lPA_Alpha <- lPA_Alpha[!(lPA_Alpha$Session %in% c(20, 21, 22)),]
lPA_Alpha <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Diversity, lPA_Alpha$Taxa, lPA_Alpha$Group, lPA_Alpha$SC2, lPA_Alpha$Type)) 
lPA_Alpha <- lPA_Alpha[which(lapply(lPA_Alpha, nrow) != 0)]     


#out.rawAgg <- extractdat(out.rawAgg)
#out.rawAgg_d <- extractdat(out.rawAgg_d)
#out.rawAgg_u <- extractdat(out.rawAgg_u)

#out.rawAbu <- lapply(out.rawAbu, function(x) {extractdat(x)})
#out.rawAbu <- lapply(out.rawAbu, function(x) {Map(cbind, x, "Metric" = names(x))})
#out.rawAbu <- lapply(out.rawAbu, function(x) {do.call(rbind, x)})
#out.rawAbu <- lapply(out.rawAbu, function(x) {as.numeric(x$Session <- str_split_fixed(x$Metric, " ", 2)[,1]);x
#                                              x$Order.q <- str_split_fixed(x$Metric, " ", 2)[,2];x})
#out.rawAbu <- lapply(out.rawAbu, function(x) {split(x, list(x$Assemblage, x$Order.q))})


# Test: ========================================================================
# H0: The time series is non-stationary.
# HA: The time series is stationary.
# seems function adf.test assumes time points in order without need for the var
# Create a list of dataframes


resMR <- lapply(lPA_Alpha, function(x) {adf.test(x$Estimate)})



#tout.rawAgg <- performtest(out.rawAgg)
#tout.rawAggTA <- lapply(out.rawAgg, function(x) {adf.test(x$m)})

#tout.rawAgg_d <- performtest(out.rawAgg_d)
#tout.rawAggTA_d <- lapply(out.rawAgg_d, function(x) {adf.test(x$m)})

#tout.rawAgg_u <- performtest(out.rawAgg_u)
#tout.rawAggTA_u <- lapply(out.rawAgg_u, function(x) {adf.test(x$m)})

#out.rawAbu[["fishmac"]] <- lapply(out.rawAbu[["fishmac"]], function(x) {x <- x[as.numeric(x$Session) <20,]})
#out.rawAbu[["fish"]] <- lapply(out.rawAbu[["fish"]], function(x) {x <- x[as.numeric(x$Session) <20,]})

#tout.rawAbu <- lapply(out.rawAbu, function(x) {performtest(x)})
#adf.test(out.rawAbu$fish$Caura_d.q0$qTD)
#out.rawAbu$fish$Caura_d.q0$Session   # OK
#out.rawAbu$fish$Caura_d.q0$qTD       # OK

#tout.rawAbuTA <- lapply(out.rawAbu, function(x) {lapply(x, function(y) {adf.test(y$m)})})

Code <-names(resMR)
restodf <- function(r){
    l <- bind_rows(map_dfr(r, ~data.frame(
    Statistic = .x$statistic,
    Pval = .x$p.value)))
    l$Code <- Code
    return(l)
}

resMR <- restodf(resMR)
resMR <- separate(resMR, Code, into = paste0("Column", 1:6), sep = "\\.")


#rout.rawAgg <- restodf(tout.rawAgg, group="NR")
#rout.rawAggTA <- restodf(tout.rawAggTA, group="NR")
#rout.rawAgg_d <- restodf(tout.rawAgg_d, group="Disturbed")
#rout.rawAggTA_d <- restodf(tout.rawAggTA_d, group="Disturbed")
#rout.rawAgg_u <- restodf(tout.rawAgg_u, group="Undisturbed")
#rout.rawAggTA_u <- restodf(tout.rawAggTA_u, group="Undisturbed")

#sites <- unique(names(tout.rawAbu$fish))
#rout.rawAbu <- lapply(tout.rawAbu, function(x) {bind_rows(map_dfr(x, ~data.frame(
#  Statistic = .x$statistic,
#  Pval = .x$p.value)))})
#rout.rawAbu <- lapply(rout.rawAbu, function(x) {x$Site <- sites;x})

#rout.rawAbuTA <- lapply(tout.rawAbuTA, function(x) {bind_rows(map_dfr(x, ~data.frame(
#  Statistic = .x$statistic,
#  Pval = .x$p.value)))})
#rout.rawAbuTA <- lapply(rout.rawAbuTA, function(x) {x$Site <- sites;x})
#rout.rawAbuTA <- lapply(rout.rawAbuTA, function(x) {x[x$Site %like% "q0",]})


resMR$Sig <- ifelse(resMR$Pval > 0.05, 0, 1)
resMR$TaxaOrderQ <- paste0(resMR$Column3, " ", "Order q=", resMR$Column1)
resMR <- subset(resMR, resMR$Column1 != 1 & resMR$Column4 =="Northern Range" & resMR$Column5=="Cmin")

resMR$TaxaOrderQ
levs <- c( "Fish (Trait) Order q=2","Fish (Trait) Order q=0",
           "Fish Order q=2","Fish Order q=0",
           "Invertebrates Order q=2","Invertebrates Order q=0",
           "Diatoms Order q=2","Diatoms Order q=0")
labs <- levs
resMR$TaxaOrderQ <- factor(resMR$TaxaOrderQ, levels = levs,
                       labels = labs)

# Agg: -------------------------------------------------------------------------
#agg_res <- as.data.frame(rbind(rout.rawAgg, rout.rawAgg_d, rout.rawAgg_u))
#agg_resTA <- as.data.frame(rbind(rout.rawAggTA, rout.rawAggTA_d, rout.rawAggTA_u))

#agg_res$Sig <- ifelse(agg_res$Pval > 0.05, 0, 1)
#agg_resTA$Sig <- ifelse(agg_resTA$Pval > 0.05, 0, 1)
#
#agg_res <- agg_res[agg_res$Taxa %like% "q2",]
#agg_resTA <- agg_resTA[agg_resTA$Taxa %like% "q2",]
#unique(agg_res$Taxa)
#agg_res$Taxa <- factor(agg_res$Taxa, levels = levs[c(1:5)],
#                                      labels = labs[c(1:5)])
#agg_res$Group <- factor(agg_res$Group, levels = c("NR", "Disturbed", "Undisturbed"),
#                       labels = c("NR", "Disturbed", "Undisturbed"))
#agg_resTA$Taxa <- factor(agg_resTA$Taxa, levels = levs[c(1:5)],
#                       labels = labs[c(1:5)])
#agg_resTA$Group <- factor(agg_resTA$Group, levels = c("NR", "Disturbed", "Undisturbed"),
#                        labels = c("NR", "Disturbed", "Undisturbed"))

# Community abundance: ---------------------------------------------------------
#abu_com_res <- as.data.frame(do.call(rbind, Map(cbind, rout.rawAbu, "Taxa" = names(rout.rawAbu))))
##abu_com_resTA <- as.data.frame(do.call(rbind, Map(cbind, rout.rawAbuTA, "Taxa" = names(rout.rawAbuTA))))
#abu_com_resTA$Site <- gsub("q0", "m", abu_com_resTA$Site)
#abu_com_res <- as.data.frame(rbind(abu_com_res, abu_com_resTA))
#abu_com_res$Sig <- ifelse(abu_com_res$Pval > 0.05, 0, 1)
#abu_com_res <- abu_com_res[abu_com_res$Taxa %in% c("fish", "bi", "dia"),]
#sum(abu_com_res$Sig==1)  # 22 out of 144
#13/144*100               # 9 %
#unique(abu_com_res$Site[abu_com_res$Sig==1])
#abu_com_res$DU <- ifelse(abu_com_res$Site %like% "_d", "Disturbed", "Undisturbed")
#abu_com_res$Taxa <- factor(abu_com_res$Taxa, levels = c("fish", "bi", "dia"),
#                           labels = c("Fish", "Invertebrates", "Diatoms"))

################################################################################
# Plot =========================================================================

(MR_NR <- ggplot(resMR, aes(x = Statistic, y = TaxaOrderQ, color=as.factor(Sig))) +
  geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
  labs(title = "Metacommunity Regulation (Northern Range)",
       x = "ADF Statistic",
       y = "") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0, face = "bold"),
     strip.text = element_text(hjust = 0, size = 10),
     #strip.background = element_blank()
     )+
  scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none"))


#(fp_agg <- ggplot(agg_res, aes(x = Statistic, y = Taxa, color=as.factor(Sig))) +
#    geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
#    labs(title = "Aggregated Abundance Data",
#         x = "Statistic",
#         y = "") +
#    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
#    theme_minimal() +
#    theme(
#      plot.title = element_text(hjust = 0.5, face = "bold"),
#      strip.text = element_text(hjust = 0, size = 10),
#      #strip.background = element_blank()
#    )+
#   scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")+
#    facet_wrap(~Group, ncol = 1, labeller = labeller(NULL)))

#(fp_aggTA <- ggplot(agg_resTA, aes(x = Statistic, y = Taxa, color=as.factor(Sig))) +
#    geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
#    labs(title = "Total Aggregated Abundance Data",
#         x = "Statistic",
#         y = "") +
#    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
#    theme_minimal() +
#    theme(
#      plot.title = element_text(hjust = 0.5, face = "bold"),
#      strip.text = element_text(hjust = 0, size = 10),
      #strip.background = element_blank()
#    )+
#    scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")+
#    facet_wrap(~Group, ncol = 1, labeller = labeller(NULL)))

#abu_com_res0 <- abu_com_res[abu_com_res$Site %like% "q0",]
#(fp_comabu0 <- ggplot(abu_com_res0, aes(x = Statistic, y = Taxa, color=as.factor(Sig))) +
#    geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
#    labs(title = "Species Richness",
#         x = "Statistic",
#         y = "") +
#    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
#    theme_minimal() +
#    theme(
#      plot.title = element_text(hjust = 0.5, face = "bold"),
#      strip.text = element_text(hjust = 0, size = 10),
      #strip.background = element_blank()http://127.0.0.1:28323/graphics/plot_zoom_png?width=1184&height=861
#    )+
#    scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")+
#    facet_wrap(~DU, ncol = 1, labeller = labeller(NULL)))

#abu_com_res2 <- abu_com_res[abu_com_res$Site %like% "q2",]
#(fp_comabu2 <- ggplot(abu_com_res2, aes(x = Statistic, y = Taxa, color=as.factor(Sig))) +
#    geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
#    labs(title = "Evenness (q2)",
#         x = "Statistic",
#         y = "") +
#    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
#    theme_minimal() +
#    theme(
#      plot.title = element_text(hjust = 0.5, face = "bold"),
#      strip.text = element_text(hjust = 0, size = 10),
 #     #strip.background = element_blank()http://127.0.0.1:28323/graphics/plot_zoom_png?width=1184&height=861
#    )+
##    scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")+
#    facet_wrap(~DU, ncol = 1, labeller = labeller(NULL)))

#abu_com_resm <- abu_com_res[abu_com_res$Site %like% ".m",]
#(fp_comabum <- ggplot(abu_com_resm, aes(x = Statistic, y = Taxa, color=as.factor(Sig))) +
#    geom_point(aes(x = Statistic, color=as.factor(Sig)), size = 3) +
#    labs(title = "Total Abundance",
#         x = "Statistic",
#         y = "") +
#    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
#    theme_minimal() +
#    theme(
#      plot.title = element_text(hjust = 0.5, face = "bold"),
#     strip.text = element_text(hjust = 0, size = 10),
#      #strip.background = element_blank()http://127.0.0.1:28323/graphics/plot_zoom_png?width=1184&height=861
#    )+
#    scale_color_manual(labels=c("0", "1"), breaks = c("0","1"), values=c("gray", "black"), guide = "none")+
#    facet_wrap(~DU, ncol = 1, labeller = labeller(NULL)))
#fp_comabu <- fp_comabu0 + fp_comabu2 + fp_comabum
#fp_comabu

# Save =========================================================================
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
ggsave(paste0(plot_path, "/ForestPlot_main_C5.jpg"), MR_NR, width = 6, height = 6)
#ggsave(paste0(plot_path, "/ForestPlottaxonomicAggData.jpg"), fp_agg, width = 6, height = 8)
#ggsave(paste0(plot_path, "/ForestPlottaxonomicTotalAggData.jpg"), fp_aggTA, width = 6, height = 8)
#ggsave(paste0(plot_path, "/ForestPlottaxonomicComAbuData.jpg"), fp_comabu, width = 12, height = 8)

# End of script#################################################################
################################################################################
