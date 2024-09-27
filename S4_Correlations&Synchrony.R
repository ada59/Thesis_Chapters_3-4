################################################################################
# S4: Correlations & Synchrony 
# AFE
# June 2024
################################################################################


################################################################################
# LIBRARIES ====================================================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
library(patchwork)
library(Hmisc)
library(corrplot)
library(zoo)

rm(list=ls())
getwd()


################################################################################
# DATA: ========================================================================
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
Results_BetaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_BetaDivRData"

load(paste0(Results_AlphaDivRData, "/S3_lPA_Alpha.RData"))
load(paste0(Results_BetaDivRData, "/S3_lBetaLag.RData"))


################################################################################
# CORRPLOTS: ===================================================================
PlotsPath <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"

# Chapter 5: ===================================================================
alpha_corr_C5 <- subset(lPA_Alpha, lPA_Alpha$Group =="Northern Range" & lPA_Alpha$SC2=="Cmax") # or Cmax
beta_corr_C5 <- subset(lBetaLag, lBetaLag$Group =="Northern Range" & lBetaLag$SC2=="Cmax")     # or Cmax

levsc <- c( "D q0","D q1","D q2", 
            "I q0","I q1","I q2",
            "F q0","F q1","F q2",
            "F(T) q0","F(T) q1", "F(T) q2")  # order in corrplot

format_corr <- function(x){
  x$Taxa <- substr(x$Taxa, 1, 1)
  x$Taxa[x$Diversity=="FD_tau"] <- "F(T)"
  x$Order.q <- paste0("q", x$Order.q)
  x$Label <- paste0(x$Taxa, " ", x$Order.q)
  x$Label <- factor(x$Label, levels = levsc, labels=levsc)
  newx <- within(x, rm(SC, Size, Method, s.e., Diversity,
                       LCL, UCL, Group, SC2, Tau, Type, Order.q, Taxa))
  newx <- spread(newx, key="Label", value="Estimate")
  return(newx)
}

alpha_corr_C5 <- format_corr(alpha_corr_C5) 
alpha_corr_C5 <- subset(alpha_corr_C5 , alpha_corr_C5$Session < 20) # rm addition Fish obs

beta_corr_C5 <- format_corr(beta_corr_C5) 
beta_corr_C5 <- subset(beta_corr_C5 , beta_corr_C5$s1 < 20)                     # rm addition Fish obs
beta_corr_C5 <- subset(beta_corr_C5 , beta_corr_C5$s2 < 20)                     # rm addition Fish obs
beta_corr_C5 <- within(beta_corr_C5, rm(Dataset, Metric, s1, s2, Lag, LagSqrt)) # last step format


#hist.data.frame(alpha_corr_C5[,c(2:13)]) 
#hist.data.frame(alpha_corr_C5[,c(2:13)], nclass=10) # pearson
#hist.data.frame(beta_corr_C5[,c(2:13)])             # pearson
#hist(alpha_corr_C5[,5])


# C5 Supplementary:-------------------------------------------------------------

cor_ordersq <- function(x){
  l <- list()
  diatoms <- x[,c(2:4)]
  invertebrates <- x[,c(5:7)]
  fish <- x[,c(8:10)]
  fish_trait <- x[,c(11:13)]
  
  cor_diatoms <- cor(diatoms)
  cor_invertebrates <- cor(invertebrates)
  cor_fish <- cor(fish)
  cor_fish_trait <- cor(fish_trait)
  l <- list(cor_diatoms, cor_invertebrates, cor_fish, cor_fish_trait)
  l <- lapply(l, function(x) {colnames(x) <- c("q0", "q1", "q2");x
  rownames(x) <- c("q0", "q1", "q2");x})
  return(l)
}
alpha_corr_C5_sm <- cor_ordersq(alpha_corr_C5)
beta_corr_C5_sm <- cor_ordersq(beta_corr_C5)

#alpha:
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotAlpha_smDia.png"), height=250, width=250)
corrplot(alpha_corr_C5_sm[[1]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2,  diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotAlpha_smInv.png"), height=250, width=250)
corrplot(alpha_corr_C5_sm[[2]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2,  diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotAlpha_smFish.png"), height=250, width=250)
corrplot(alpha_corr_C5_sm[[3]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2, diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotAlpha_smFishT.png"), height=250, width=250)
corrplot(alpha_corr_C5_sm[[4]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2, diag=F)
dev.off()

# beta:
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotBeta_smDia.png"), height=250, width=250)
corrplot(beta_corr_C5_sm[[1]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2,  diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotBeta_smInv.png"), height=250, width=250)
corrplot(beta_corr_C5_sm[[2]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2,  diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotBeta_smFish.png"), height=250, width=250)
corrplot(beta_corr_C5_sm[[3]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2, diag=F)
dev.off()
png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotBeta_smFishT.png"), height=250, width=250)
corrplot(beta_corr_C5_sm[[4]], method = "number",type="lower", 
         tl.col="black", number.cex = 1.5, tl.cex=1.2, diag=F)
dev.off()


# C5 Main:----------------------------------------------------------------------
rm_q1 <- c("D q1","I q1","F q1","F(T) q1")
alpha_corr_C5_main <- alpha_corr_C5[,!names(alpha_corr_C5) %in% rm_q1]
beta_corr_C5_main <- beta_corr_C5[,!names(beta_corr_C5) %in% rm_q1]

cor_alpha <- cor(alpha_corr_C5_main[,c(2:9)], method = "pearson")
cor_beta <- cor(beta_corr_C5_main[,c(2:9)], method = "pearson")

png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotAlphaCmax.png"), height=650, width=650)
corrplotAlpha <- corrplot.mixed(cor_alpha, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.2, cl.cex=1.2)
mtext(paste0("A) α diversity"), at=1.3, line=0.5, cex=2)
dev.off()

png(file = paste0(PlotsPath,"/Chapter5/S4_C5_corrplotBetaCmax.png"), height=650, width=650)
corrplotBeta <- corrplot.mixed(cor_beta, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.2, cl.cex=1.2)
mtext(paste0("B) β diversity"), at=1.3, line=0.5, cex=2)
dev.off()



# Chapter 6: ===================================================================

# Change the following 4 lines for obtaining figures for Cmax
alpha_corr_Dis <- subset(lPA_Alpha, lPA_Alpha$Group =="Disturbed" & lPA_Alpha$SC2=="Cmin")     # or Cmax 
alpha_corr_Undis <- subset(lPA_Alpha, lPA_Alpha$Group =="Undisturbed" & lPA_Alpha$SC2=="Cmin") # or Cmax 

beta_corr_Dis <- subset(lBetaLag, lBetaLag$Group == "Disturbed" & lBetaLag$SC2=="Cmin")        # or Cmax 
beta_corr_Undis <- subset(lBetaLag, lBetaLag$Group == "Undisturbed" & lBetaLag$SC2=="Cmin")    # or Cmax       

levsc <- c( "D q0","D q1","D q2", 
            "I q0","I q1","I q2",
            "F q0","F q1","F q2",
            "F(T) q0","F(T) q1", "F(T) q2")  # order in corrplot

format_corr <- function(x){
  x$Taxa <- substr(x$Taxa, 1, 1)
  x$Taxa[x$Diversity=="FD_tau"] <- "F(T)"
  x$Order.q <- paste0("q", x$Order.q)
  x$Label <- paste0(x$Taxa, " ", x$Order.q)
  x$Label <- factor(x$Label, levels = levsc, labels=levsc)
  newx <- within(x, rm(SC, Size, Method, s.e., Diversity,
                       LCL, UCL, Group, SC2, Tau, Type, Order.q, Taxa))
  newx <- spread(newx, key="Label", value="Estimate")
  return(newx)
}

alpha_corr_Dis <- format_corr(alpha_corr_Dis) 
alpha_corr_Dis <- subset(alpha_corr_Dis , alpha_corr_Dis$Session < 20) 
alpha_corr_Undis <- format_corr(alpha_corr_Undis) 
alpha_corr_Undis <- subset(alpha_corr_Undis , alpha_corr_Undis$Session < 20) 

beta_corr_Dis <- format_corr(beta_corr_Dis) 
beta_corr_Dis <- subset(beta_corr_Dis, beta_corr_Dis$s1 < 20) 
beta_corr_Dis <- subset(beta_corr_Dis, beta_corr_Dis$s2 < 20) 
beta_corr_Dis_Syn <- within(beta_corr_Dis, rm(concat_sessions, Dataset, Metric, s1, s2, LagSqrt))
beta_corr_Dis <- within(beta_corr_Dis, rm(Dataset, Metric, s1, s2, Lag, LagSqrt)) 

beta_corr_Undis <- format_corr(beta_corr_Undis) 
beta_corr_Undis <- subset(beta_corr_Undis, beta_corr_Undis$s1 < 20) 
beta_corr_Undis <- subset(beta_corr_Undis, beta_corr_Undis$s2 < 20)   
beta_corr_Undis_Syn <- within(beta_corr_Undis, rm(concat_sessions, Dataset, Metric, s1, s2, LagSqrt))
beta_corr_Undis <- within(beta_corr_Undis, rm(Dataset, Metric, s1, s2, LagSqrt)) 

#hist.data.frame(alpha_corr_Dis[,c(2:13)])     
#hist.data.frame(alpha_corr_Undis[,c(2:13)])
#hist.data.frame(beta_corr_Dis[,c(2:13)])     
#hist.data.frame(beta_corr_Undis[,c(2:13)])


# Sm, corrs cross taxa & orders of q:

#alpha_corr_Dis_main <- alpha_corr_Dis[,!names(alpha_corr_Dis) %in% rm_q1]
#alpha_corr_Undis_main <- alpha_corr_Undis[,!names(alpha_corr_Undis) %in% rm_q1]

#beta_corr_Dis_main <- beta_corr_Dis[,!names(beta_corr_Dis) %in% rm_q1]
#beta_corr_Undis_main <- beta_corr_Undis[,!names(beta_corr_Undis) %in% rm_q1]   # use if wanting to rm q1

cor_alpha_Dis <- cor(alpha_corr_Dis[,c(2:13)], method = "pearson")
cor_alpha_Undis <- cor(alpha_corr_Undis[,c(2:13)], method = "pearson")

cor_beta_Dis <- cor(beta_corr_Dis[,c(2:13)], method = "pearson")
cor_beta_Undis <- cor(beta_corr_Undis[,c(2:13)], method = "pearson")


png(file = paste0(PlotsPath,"/Chapter6/S4_C6_corrplotAlphaDis.png"), height=750, width=750) # add Cmax in name when saving
corrplotAlphaD <- corrplot.mixed(cor_alpha_Dis, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.1, cl.cex=1.2)
mtext(paste0("A) α diversity Disturbed"), at=2.7, line=0.5, cex=2)
dev.off()
png(file = paste0(PlotsPath,"/Chapter6/S4_C6_corrplotAlphaUndis.png"), height=750, width=750)
corrplotAlphaU <- corrplot.mixed(cor_alpha_Undis, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.1, cl.cex=1.2)
mtext(paste0("B) α diversity Undisturbed"), at=2.7, line=0.5, cex=2)
dev.off()
png(file = paste0(PlotsPath,"/Chapter6/S4_C6_corrplotBetaDis.png"), height=750, width=750)
corrplotBetaD <- corrplot.mixed(cor_beta_Dis, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.1, cl.cex=1.2)
mtext(paste0("C) β diversity Disturbed"), at=2.7, line=0.5, cex=2)
dev.off()
png(file = paste0(PlotsPath,"/Chapter6/S4_C6_corrplotBetaUndis.png"), height=750, width=750)
corrplotBetaU <- corrplot.mixed(cor_beta_Undis, lower = 'square', upper = 'number', tl.col="black", number.cex = 1.5, tl.cex=1.1, cl.cex=1.2)
mtext(paste0("D) β diversity Undisturbed"), at=2.7, line=0.5, cex=2)
dev.off()


################################################################################
# SYNCHRONY ====================================================================
# https://towardsdatascience.com/four-ways-to-quantify-synchrony-between-time-series-data-b99136c4a9c9
# https://www.r-bloggers.com/2021/11/how-to-perform-rolling-correlation-in-r/#:~:text=Rolling%20Correlation%20in%20R%2C%20Correlations,It%20evolves%20over%20time.
# https://www.geeksforgeeks.org/how-to-calculate-rolling-correlation-in-r/


# Chapter 6: ===================================================================
# Alpha Diversity (alpha_corr_Dis), Beta Diversity(beta)
corr_Dis <- alpha_corr_Dis
#corr_Dis <- beta_corr_Dis_Syn

corr_Undis <- alpha_corr_Undis
#corr_Undis <- beta_corr_Undis_Syn

dataDis0 <- list(as.data.frame(corr_Dis[,c(1,2,5)]),  # dia & inv, q0
                as.data.frame(corr_Dis[,c(1,2,8)]),   # dia & fish, q0
                as.data.frame(corr_Dis[,c(1,2,11)]),  # dia & fish trait, q0
                as.data.frame(corr_Dis[,c(1,5,8)]),   # inv & fish, q0
                as.data.frame(corr_Dis[,c(1,5,11)]),  # inv & fish trait, q0
                as.data.frame(corr_Dis[,c(1,8,11)]))  # fish & fish trait, q0

dataDis2 <- list(as.data.frame(corr_Dis[,c(1,4,7)]),   # dia & inv, q2
                 as.data.frame(corr_Dis[,c(1,4,10)]),  # dia & fish, q2
                 as.data.frame(corr_Dis[,c(1,4,13)]),  # dia & fish trait, q2
                 as.data.frame(corr_Dis[,c(1,7,10)]),  # inv & fish, q2
                 as.data.frame(corr_Dis[,c(1,7,13)]),  # inv & fish trait, q2
                 as.data.frame(corr_Dis[,c(1,10,13)])) # fish & fish trait, q2

dataUndis0 <- list(as.data.frame(corr_Undis[,c(1,2,5)]), # dia & inv, q0
                 as.data.frame(corr_Undis[,c(1,2,8)]),   # dia & fish, q0
                 as.data.frame(corr_Undis[,c(1,2,11)]),  # dia & fish trait, q0
                 as.data.frame(corr_Undis[,c(1,5,8)]),   # inv & fish, q0
                 as.data.frame(corr_Undis[,c(1,5,11)]),  # inv & fish trait, q0
                 as.data.frame(corr_Undis[,c(1,8,11)]))  # fish & fish trait, q0

dataUndis2 <- list(as.data.frame(corr_Undis[,c(1,4,7)]), # dia & inv, q2
                 as.data.frame(corr_Undis[,c(1,4,10)]),  # dia & fish, q2
                 as.data.frame(corr_Undis[,c(1,4,13)]),  # dia & fish trait, q2
                 as.data.frame(corr_Undis[,c(1,7,10)]),  # inv & fish, q2
                 as.data.frame(corr_Undis[,c(1,7,13)]),  # inv & fish trait, q2
                 as.data.frame(corr_Undis[,c(1,10,13)])) # fish & fish trait, q2


# Alpha:------------------------------------------------------------------------
alphaDis0 <- lapply(dataDis0, function(x) {rollapply(x, width=8, function(x) cor(as.numeric(x[,2]), as.numeric(x[,3])),
                                                    by.column=FALSE)}) # rolling cors with 8 data point windows for alpha
dis0 <- as.data.frame(do.call(cbind, alphaDis0))
alphaDis2 <- lapply(dataDis2, function(x) {rollapply(x, width=8, function(x) cor(as.numeric(x[,2]), as.numeric(x[,3])),
                                                    by.column=FALSE)})
dis2 <- as.data.frame(do.call(cbind, alphaDis2))

alphaUndis0 <- lapply(dataUndis0, function(x) {rollapply(x, width=8, function(x) cor(as.numeric(x[,2]), as.numeric(x[,3])),
                                                        by.column=FALSE)})
undis0 <- as.data.frame(do.call(cbind, alphaUndis0))
alphaUndis2 <- lapply(dataUndis2, function(x) {rollapply(x, width=8, function(x) cor(as.numeric(x[,2]), as.numeric(x[,3])),
                                                        by.column=FALSE)})
undis2 <- as.data.frame(do.call(cbind, alphaUndis2))


MovCorAlpha <- as.data.frame(rbind(dis0, undis0, dis2, undis2))
xaxis <- rep(1:12, times=4)  
dfsynAlpha <- data.frame("Disturbance"=rep(c("Disturbed", "Undisturbed",
                                             "Disturbed", "Undisturbed"), each=12), 
                         "MovCor"=MovCorAlpha,
                         "Order.q"=rep(c(0,2), each=24), "Time"=xaxis)
names(dfsynAlpha)
names(dfsynAlpha) <- c("Disturbance", 
                       "Diatoms-Invertebrates",
                       "Diatoms-Fish",
                       "Diatoms-Fish(Trait)",
                       "Invertebrates-Fish",
                       "Invertebrates-Fish(Trait)",
                       "Fish-Fish(Trait)",
                       "Order.q",
                       "Time")

dfsynAlpha <- gather(dfsynAlpha, key="Pair", value="MovCor", -c(1,8,9))
str(dfsynAlpha)
levsSyn <- c("Diatoms-Invertebrates",
             "Diatoms-Fish",
             "Diatoms-Fish(Trait)",
             "Invertebrates-Fish",
             "Invertebrates-Fish(Trait)",
             "Fish-Fish(Trait)")
dfsynAlpha$Pair <- factor(dfsynAlpha$Pair, levels = levsSyn, labels=levsSyn)

dfsynAlpha0 <- subset(dfsynAlpha, dfsynAlpha$Order.q==0)
dfsynAlpha2 <- subset(dfsynAlpha, dfsynAlpha$Order.q==2)

(S4_C6_SynAlpha <- ggplot(data=dfsynAlpha, aes(x=Time, y=MovCor, color=Disturbance, linetype=as.factor(Order.q)))+
    labs(x="Time window (8 ts)", y="rolling (r) for α diversity", linetype="Order.q", color="Disturbance",
         title="Cross-taxon synchrony in α-diversity change")+
    theme_classic()+
    theme(
      strip.text = element_text(size = 11), 
      plot.title = element_text(size = 15, hjust = 0.5, vjust = 1.5)
    )+
    ylim(-1,1)+
    geom_point(data=dfsynAlpha, aes(x=Time, y=MovCor, color=Disturbance), size=1.5)+
    geom_line(data=dfsynAlpha, aes(x=Time, y=MovCor, linetype=as.factor(Order.q), color=Disturbance), size=1)+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1, fill = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    scale_linetype_manual(values=c("dotted", "solid"), breaks=c("0", "2"), labels=c("q0", "q2"))+
    scale_color_manual(values=c("#CE780F", "#0A7D9D"))+
    facet_wrap(~Pair))

ggsave(S4_C6_SynAlpha, filename= paste0(PlotsPath, "/Chapter6/S4_C6_SynAlpha.jpg"), width = 10, height = 6) # add Cmax for saving the Cmax version



# Beta: ------------------------------------------------------------------------

lag_correlations <- function(data, var1=NULL, var2=NULL, max_lag = NULL, lag_range = NULL) {
  results <- list()
  for (start_lag in 1:(max_lag - lag_range + 1)) {
    end_lag <- start_lag + lag_range - 1
    lagged_data <- data %>%
      filter(Lag >= start_lag & Lag <= end_lag)
    correlation <- cor(lagged_data[[var1]], lagged_data[[var2]])
    results[[paste0("Lag_", start_lag, "_to_", end_lag)]] <- correlation
  }
  return(results)
} 

betaDis0 <- lapply(dataDis0, function(x) {lag_correlations(x, var1 = names(x)[2], var2 = names(x)[3], max_lag = 18, lag_range = 7)})
betaDis0 <- as.data.frame(do.call(cbind, betaDis0))

betaUndis0 <- lapply(dataUndis0, function(x) {lag_correlations(x, var1 = names(x)[2], var2 = names(x)[3], max_lag = 18, lag_range = 7)})
betaUndis0 <- as.data.frame(do.call(cbind, betaUndis0))

betaDis2 <- lapply(dataDis2, function(x) {lag_correlations(x, var1 = names(x)[2], var2 = names(x)[3], max_lag = 18, lag_range = 7)})
betaDis2 <- as.data.frame(do.call(cbind, betaDis2))

betaUndis2 <- lapply(dataUndis2, function(x) {lag_correlations(x, var1 = names(x)[2], var2 = names(x)[3], max_lag = 18, lag_range = 7)})
betaUndis2 <- as.data.frame(do.call(cbind, betaUndis2))

MovCorBeta <- as.data.frame(rbind(betaDis0, 
                                  betaUndis0, 
                                  betaDis2, 
                                  betaUndis2))  # !!!use the code above from here, just need a few re-namings
MovCorBeta <- MovCorBeta %>%
  unnest(c(V1, V2, V3, V4, V5, V6))

xaxis <- rep(1:12, times=4)  
dfsynBeta <- data.frame("Disturbance"=rep(c("Disturbed", "Undisturbed",
                                             "Disturbed", "Undisturbed"), each=12), 
                         "MovCor"=MovCorBeta,
                         "Order.q"=rep(c(0,2), each=24), "Time"=xaxis)

names(dfsynBeta)
names(dfsynBeta) <- c("Disturbance", 
                       "Diatoms-Invertebrates",
                       "Diatoms-Fish",
                       "Diatoms-Fish(Trait)",
                       "Invertebrates-Fish",
                       "Invertebrates-Fish(Trait)",
                       "Fish-Fish(Trait)",
                       "Order.q",
                       "Time")

dfsynBeta <- gather(dfsynBeta, key="Pair", value="MovCor", -c(1,8,9))
str(dfsynBeta)
levsSyn <- c("Diatoms-Invertebrates",
             "Diatoms-Fish",
             "Diatoms-Fish(Trait)",
             "Invertebrates-Fish",
             "Invertebrates-Fish(Trait)",
             "Fish-Fish(Trait)")
dfsynBeta$Pair <- factor(dfsynBeta$Pair, levels = levsSyn, labels=levsSyn)

dfsynBeta0 <- subset(dfsynBeta, dfsynBeta$Order.q==0)
dfsynBeta2 <- subset(dfsynBeta, dfsynBeta$Order.q==2)
str(dfsynBeta)


(S4_C6_SynBeta <- ggplot(data=dfsynBeta, aes(x=Time, y=MovCor, color=Disturbance, linetype=as.factor(Order.q)))+
    labs(x="Time window (7 ts of time lag)", y="rolling (r) for β diversity", linetype="Order.q", color="Disturbance",
         title="Cross-taxon synchrony in β-diversity change")+
    theme_classic()+
    theme(
      strip.text = element_text(size = 11), 
      plot.title = element_text(size = 15, hjust = 0.5, vjust = 1.5)
    )+
    ylim(-1,1)+
    geom_point(data=dfsynBeta, aes(x=Time, y=MovCor, color=Disturbance), size=1.5)+
    geom_line(data=dfsynBeta, aes(x=Time, y=MovCor, linetype=as.factor(Order.q), color=Disturbance), size=1)+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.1, fill = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    scale_linetype_manual(values=c("dotted", "solid"), breaks=c("0", "2"), labels=c("q0", "q2"))+
    scale_color_manual(values=c("#CE780F", "#0A7D9D"))+
    facet_wrap(~Pair))

ggsave(S4_C6_SynBeta, filename= paste0(PlotsPath, "/Chapter6/S4_C6_SynBeta.jpg"), width = 10, height = 6) # add Cmax for saving the Cmax version


# End of script ################################################################
################################################################################

