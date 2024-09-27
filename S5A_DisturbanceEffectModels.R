################################################################################
# S5 Linear models (D/U)
# AFE
# June 2024
################################################################################


library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(AICcmodavg)
library(broom)
library(data.table)
library(scales)
library(lme4)
library(MASS)
library(lmerTest)
library(patchwork)

rm(list=ls())


################################################################################
# NOTES: =======================================================================
# Follow code from Practical Stats BL3320 (GLMs)
# Other possible structure: lmer(diversity ~ disturbance + (1 | session), 
# data = data) [mixed model with random intercept]
# In GM: provide PCA with D vs U, and graph of av. anthropogenic litter vs time 
# grouped by D/U
################################################################################


################################################################################
# read data ====================================================================
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
Results_BetaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_BetaDivRData"

load(paste0(Results_AlphaDivRData, "/S3_lPA_Alpha.RData"))
load(paste0(Results_BetaDivRData, "/S3_lBetaLag.RData"))

lPA_Alpha <- lPA_Alpha[!lPA_Alpha$Order.q == 1,]               # rm order q1
lPA_Alpha <- lPA_Alpha[!lPA_Alpha$Group == "Northern Range",]  # rm NR data
lPA_Alpha$Group <- as.factor(lPA_Alpha$Group)                  # disturbance as factor
lPA_Alpha$Taxa[lPA_Alpha$Taxa=="Fish (Trait)"] <- "FishTrait"

lBetaLag <- lBetaLag[!lBetaLag$Order.q == 1,]                  # rm order q1
lBetaLag <- lBetaLag[!lBetaLag$Group == "Northern Range",]     # rm NR data
lBetaLag$Group <- as.factor(lBetaLag$Group)                    # disturbance as factor
lBetaLag$Taxa[lBetaLag$Taxa=="Fish (Trait)"] <- "FishTrait"

# Alpha: -----------------------------------------------------------------------
str(lPA_Alpha)
alphalist23 <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Taxa, lPA_Alpha$SC2)) 
alphalist23 <- alphalist23[which(lapply(alphalist23, nrow) != 0)]  # 16, although if focusing on Cmin & q0 & q2, it's 8.


# Beta: ------------------------------------------------------------------------
betalist23 <- split(lBetaLag, list(lBetaLag$Order.q, lBetaLag$Taxa, lBetaLag$SC2)) 
betalist23 <- betalist23[which(lapply(betalist23, nrow) != 0)]     # idem above


################################################################################
# models: ======================================================================
# Interaction: does the relationship between diversity & time depend on disturbance?

run_models <- function(dt, time_pred=NULL){
  
  m1 <- lm(Estimate ~ Group*time_pred, data = dt)
  m2 <- lm(Estimate ~ Group+time_pred, data = dt)
  m3 <- lm(Estimate ~ Group, data = dt)
  m4 <- lm(Estimate ~ time_pred, data = dt)
  
  models <- list("Interaction"=m1, "Additive"=m2, "Disturbance"=m3, "Time"=m4)

  return(models)
} 

alpha_mod23 <- lapply(alphalist23, function(x) {run_models(dt=x, time_pred = x$Session)})
beta_mod23 <- lapply(betalist23, function(x) {run_models(dt=x, time_pred = x$LagSqrt)})


################################################################################
# store aic tabs ===============================================================
table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"

R2_AIC_Tab <- function(dt){
  R2 <- lapply(dt, function(x) {lapply(x, function(y) {round(summary(y)$adj.r.squared, 2)})}) # R2
  R2 <- lapply(R2, function(x) {as.data.frame(do.call(rbind, x))})
  compare <- lapply(dt, function(x) {aictab(x)})                                              # rank by AICc
  compare <- mapply(function(x,y) {x$adjR2 <- y$V1[match(x$Modnames, rownames(y))];x}, compare, R2, SIMPLIFY = F)
  raw_tabs <- do.call(rbind, compare)
  
  raw_tabs$Code <- rownames(raw_tabs) # organise & layout for SM
  rownames(raw_tabs) <- 1:nrow(raw_tabs)
  raw_tabs <- separate(raw_tabs, Code, into = paste0("Column", 1:4), sep = "\\.")          # re-assign names
  raw_tabs <- raw_tabs[, ! names(raw_tabs) %in% c("Column4")]                              # not needed
  names(raw_tabs)[names(raw_tabs)=="Column1"] <- "Orderq"
  names(raw_tabs)[names(raw_tabs)=="Column2"] <- "Taxa"
  names(raw_tabs)[names(raw_tabs)=="Column3"] <- "Coverage"
  raw_tabs <- raw_tabs %>% relocate(c(Taxa, Coverage, Orderq), .before = Modnames)
  raw_tabs <- raw_tabs[, ! names(raw_tabs)  %in% c("K", "ModelLik", "AICcWt", "Cum.Wt")]   # unnecessary vals
  raw_tabs[,c(5:8)] <- sapply(raw_tabs[,c(5:8)], function(x) {round(x, 2)})                # all rounded
  
  raw_tabs$Taxa <- factor(raw_tabs$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                    "Fish", "FishTrait"))                  # correct order
  raw_tabs <- raw_tabs[order(raw_tabs$Taxa), ]                                             # correct order
  raw_tabs <- raw_tabs %>% arrange(Taxa, desc(Coverage), Orderq)                           # correct order
  formatted_tabs <- raw_tabs

  return(formatted_tabs)
  } ### consider model averaging if there's time

alpha_aic23 <- R2_AIC_Tab(alpha_mod23)
beta_aic23 <- R2_AIC_Tab(beta_mod23)

#write.csv(alpha_aic23, file=paste0(table_path, "/S5A_Alpha_AIC_Tables_C6.csv"), row.names = F)
#write.csv(beta_aic23, file=paste0(table_path, "/S5A_Beta_AIC_Tables_C6.csv"), row.names = F)


################################################################################
# select best model:============================================================
# NOTE: Delta_AICc =0 (yet remember DeltaAICc < 2, equally plausible)
# Alternatively use backward stepwise, and double checj AICc, LL & R2 in tables

org <- function(x){
  x <- c(
    x[grep("Diatoms.Cmin", names(x))],
    x[grep("Diatoms.Cmax", names(x))],
    x[grep("Invertebrates.Cmin", names(x))],
    x[grep("Invertebrates.Cmax", names(x))],
    x[grep("Fish.Cmin", names(x))],
    x[grep("Fish.Cmax", names(x))],
    x[grep("FishTrait.Cmin", names(x))],
    x[grep("FishTrait.Cmax", names(x))]
  )
  return(x)
} # re-organise names

alpha_best <- subset(alpha_aic23, alpha_aic23$Delta_AICc==0)
alpha_best <- split(alpha_best, f=list(alpha_best$Orderq, alpha_best$Taxa, alpha_best$Coverage)) # best according to AICc
alpha_best <- lapply(alpha_best, function(x) {x$Modnames})                                       # name of best model

alphalist23 <- org(alphalist23)
alpha_mod23 <- org(alpha_mod23)
alpha_best <- org(alpha_best)   
identical(names(alpha_mod23), names(alpha_best)) # ensure name order in different lists is the same 

alpha_best<- mapply(function(x,y) {x[y]}, alpha_mod23, alpha_best, SIMPLIFY=T)                   # best according to AICc
alpha_sum <- lapply(alpha_best, function(x) {summary(x)})                                        # summary best

alpha_int <- lapply(alpha_mod23, function(x) {x <- x$Interaction})                               # interaction term



beta_best <- subset(beta_aic23, beta_aic23$Delta_AICc==0)
beta_best <- split(beta_best, f=list(beta_best$Orderq, beta_best$Taxa, beta_best$Coverage))      # best according to AICc
beta_best <- lapply(beta_best, function(x) {x$Modnames})                                         # name of best model

betalist23 <- org(betalist23)
beta_mod23 <- org(beta_mod23)
beta_best <- org(beta_best)   
identical(names(beta_mod23), names(beta_best))   # ensure name order in different lists is the same 

beta_best <- mapply(function(x,y) {x[y]}, beta_mod23, beta_best, SIMPLIFY=T)                      # best according to AICc
beta_sum <- lapply(beta_best, function(x) {summary(x)})                                           # summary best

beta_int <- lapply(beta_mod23, function(x) {x <- x$Interaction})                                  # interaction term


# Alternative
# Automatic stepwise model selection
#alpha_best <- lapply(alphalist23, function(x) {m1 <- lm(Estimate ~ Group*Session, data = x);x
#                                              res <- stepAIC(m1, direction = "backward")}) # used step AIC to check "manual" selection
#beta_best <- lapply(betalist23, function(x) {m1 <- lm(Estimate ~ Group*LagSqrt, data = x);x
#res <- stepAIC(m1, direction = "backward")})


################################################################################
# Assumptions: =================================================================
# http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/
# NOTE: here I check normality of residuals & homoscedasticty

alpha_sresid <- lapply(alpha_best, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)}) # standardised residuals (also alpha int)
beta_sresid <- lapply(beta_best, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)})   # standardised residuals (also beta int)

namesAlpha <- names(alpha_sresid)
namesAlpha <- gsub("\\.", " ", namesAlpha)

namesBeta <- names(beta_sresid)
namesBeta <- gsub("\\.", " ", namesBeta)


plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
#jpeg(file=paste0(plot_path,"/S5A_C6_Alpha_NormalityAssumptionDirectional.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#  }, alpha_sresid, namesAlpha, SIMPLIFY = T) 
#dev.off()  # assess normality of residuals (checked for the interaction model too)
#par(mfrow=c(4,4))
#lapply(alpha_bestmod23, function(x) {plot(x, which=2)})
#dev.off()

#jpeg(file=paste0(plot_path,"/S5A_C6_Alpha_HomoscedasticityAssumptionDirectional.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
 # plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, alpha_sresid, alpha_best, namesAlpha, SIMPLIFY = T)
#dev.off()  # assess homoscedasticity (checked for the interaction model too)


#jpeg(file=paste0(plot_path,"/S5A_C6_Beta_NormalityAssumptionDirectional.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#}, beta_sresid, namesBeta, SIMPLIFY = T) 
#dev.off()  

#jpeg(file=paste0(plot_path,"/S5A_C6_Beta_HomoscedasticityAssumptionDirectional.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
#  plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, beta_sresid, beta_best, namesBeta, SIMPLIFY = T)
#dev.off()  

# NOTE: plot(x, which = 2, main=y) normality, and which = 1 for homoscedasticity


# Coefs best model: ============================================================
retrieve_coefs <- function(best){
  coefs <- lapply(best, function(x) {as.data.frame(summary(x)$coefficients)}) # none with just time
  coefs <- lapply(coefs, function(x) {x$Term <- rownames(x);x})  
  coefs <- do.call(rbind, coefs)
  coefs$Code <- rownames(coefs)
  rownames(coefs) <- 1:nrow(coefs)
  
  coefs <- separate(coefs, Code, into = paste0("Column", 1:4), sep = "\\.") # warning OK
  names(coefs)[names(coefs)=="Column1"] <- "Orderq"
  names(coefs)[names(coefs)=="Column2"] <- "Taxa"
  names(coefs)[names(coefs)=="Column3"] <- "Coverage"
  names(coefs)[names(coefs)=="Column4"] <- "Model"
  coefs <- coefs %>% relocate(c(Taxa, Coverage, Orderq, Model, Term), .before = Estimate)
  #coefs[,c(6:9)] <- sapply(coefs[,c(6:9)], function(x) {round(x, 3)})
  return(coefs)
}


alpha_coefs <- retrieve_coefs(alpha_best) # warnings OK, or alpha_int
alpha_coefs$Term <- plyr::revalue(alpha_coefs$Term, c("(Intercept)"="Intercept",
                                                      "GroupUndisturbed"="Group (Undisturbed)",
                                                      "time_pred"="Session",
                                                      "GroupUndisturbed:time_pred"="Session (Undisturbed)")) 
beta_coefs <- retrieve_coefs(beta_best)   # warnings OK, or beta_int
beta_coefs$Term <- plyr::revalue(beta_coefs$Term, c("(Intercept)"="Intercept",
                                                      "GroupUndisturbed"="Group (Undisturbed)",
                                                      "time_pred" = "LagSqrt Session",
                                                      "GroupUndisturbed:time_pred"="LagSqrt Session (Undisturbed)")) 

#write.csv(alpha_coefs, file=paste0(table_path, "/S5A_Alpha_Coef_Tables_C6.csv"), row.names = F)
#write.csv(beta_coefs, file=paste0(table_path, "/S5A_Beta_Coef_Tables_C6.csv"), row.names = F)

#write.csv(alpha_coefs, file=paste0(table_path, "/S5A_Alpha_Int_Coef_Tables_C6.csv"), row.names = F)
#write.csv(beta_coefs, file=paste0(table_path, "/S5A_Beta_Int_Coef_Tables_C6.csv"), row.names = F)   # For interaction models


################################################################################
# Predict: =====================================================================
sumTime <- function(x){
  if(x[["rank"]]==3){
    R2adj <- rep(summary(x)$adj.r.squared, times=1000)
    sum <- summary(x)$coefficients
    p_valDU <- rep(sum[2,4], times=1000)
    p_valTime <- rep(sum[3,4], times=1000)
    summa <- as.data.frame(cbind("R2"=R2adj, "pvalDU"=p_valDU, "pvalTime"=p_valTime))
  }else{
    R2adj <- rep(summary(x)$adj.r.squared, times=1000)
    sum <- summary(x)$coefficients
    p_valDU <- rep(sum[2,4], times=1000)
    p_valDis <- rep(sum[3,4], times=500)
    p_valUndis <- rep(sum[4,4], times=500)
    summa <- as.data.frame(cbind("R2"=R2adj, "pvalDU"=p_valDU, "pvalTime"=c(p_valDis, p_valUndis))) 
  }
  return(summa)
}  # coefs Time effect
sumDis <- function(x){
    R2adj <- rep(summary(x)$adj.r.squared, times=1000)
    sum <- summary(x)$coefficients
    p_valDU <- rep(sum[2,4], times=1000)
    summa <- as.data.frame(cbind("R2"=R2adj, "pvalDU"=p_valDU, "pvalTime"=1))
  return(summa)
}   # coefs Disturbance effect

alphalist23 <- lapply(alphalist23, function(x) {names(x)[names(x)=="Session"] <- "time_pred";x})
betalist23 <- lapply(betalist23, function(x) {names(x)[names(x)=="LagSqrt"] <- "time_pred";x})

# Alpha : ----------------------------------------------------------------------
disAlpha <- c("2.Diatoms.Cmin","0.Diatoms.Cmax","2.Diatoms.Cmax","0.Invertebrates.Cmin",
         "0.Invertebrates.Cmax","0.Fish.Cmin","0.Fish.Cmax")

pred_best <- function(bestmodel, rawdata, type=NULL, dis=NULL){
  if(type=="TIMEeffects"){
    
    time_best <- bestmodel[!grepl("Disturbance", names(bestmodel))] 
    time_best_data <- rawdata[! names(rawdata) %in% dis]
    time_x <- lapply(time_best_data, function(x) {data.frame("time_pred"=c(seq(from = min(x$time_pred), to = max(x$time_pred), length= 500),
                                                                           seq(from = min(x$time_pred), to = max(x$time_pred), length= 500)),
                                                             "Group"=rep(c("Disturbed", "Undisturbed"), each=500))})
    time_sum <- lapply(time_best, function(x) {sumTime(x)})
 
     }else{
       
    time_best <- bestmodel[grepl("Disturbance", names(bestmodel))] # disturbance effect only models
    time_best_data <- rawdata[names(rawdata) %in% dis]
    time_x <- lapply(time_best_data, function(x) {data.frame("Group"=rep(c("Disturbed", "Undisturbed"), each=500))})
    time_sum <- lapply(time_best, function(x) {sumDis(x)})
    
  }
  
  time_pred <- mapply(function(x,y,z) {data.frame(y, predict(x, y, interval="confidence"),z)},
                      time_best, time_x, time_sum, SIMPLIFY = F) # predict
  pred <- as.data.frame(do.call(rbind, time_pred))
  pred$Taxa <- str_split_fixed(rownames(pred),"\\.", 5)[,2]
  pred$CoverageVals <- str_split_fixed(rownames(pred),"\\.", 5)[,3]
  pred$Orderq <- str_split_fixed(rownames(pred),"\\.", 5)[,1]
  pred$Model <- str_split_fixed(rownames(pred),"\\.", 5)[,4]
  rownames(pred) <- 1:nrow(pred) 
  pred$SigDU <- ifelse(pred$pvalDU < 0.05, 1, 0)
  pred$SigTime <- ifelse(pred$pvalTime < 0.05, 1, 0)
  return(pred)
}


time_models <- pred_best(alpha_best, alphalist23, type="TIMEeffects", dis=disAlpha)
dis_models <- pred_best(alpha_best, alphalist23, type="DISeffects", dis=disAlpha)

names(time_models)
names(dis_models)

pred_directional_alpha <- bind_rows(time_models, dis_models)
pred_directional_alpha$Taxa[pred_directional_alpha$Taxa == "FishTrait"] <- "Fish (Trait)"
pred_directional_alpha$Taxa <- factor(pred_directional_alpha$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                                              "Fish", "Fish (Trait)"),
                                      labels = c("Diatoms", "Invertebrates", 
                                                 "Fish", "Fish (Trait)"))
pred_directional_alpha <- pred_directional_alpha[order(pred_directional_alpha$Taxa), ] # to use for plots


# Beta:-------------------------------------------------------------------------
disBeta <- c("0.Diatoms.Cmin","2.Diatoms.Cmin","0.Diatoms.Cmax",  
             "2.Diatoms.Cmax","0.Fish.Cmin","0.Fish.Cmax",     
             "0.FishTrait.Cmin")
time_modelsB <- pred_best(beta_best, betalist23, type="TIMEeffects", dis=disBeta)
dis_modelsB <- pred_best(beta_best, betalist23, type="DISeffects", dis=disBeta)

names(time_modelsB)
names(dis_modelsB)

pred_directional_beta <- bind_rows(time_modelsB, dis_modelsB)
pred_directional_beta$Taxa[pred_directional_beta$Taxa == "FishTrait"] <- "Fish (Trait)"
pred_directional_beta$Taxa <- factor(pred_directional_beta$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                                              "Fish", "Fish (Trait)"),
                                      labels = c("Diatoms", "Invertebrates", 
                                                 "Fish", "Fish (Trait)"))
pred_directional_beta <- pred_directional_beta[order(pred_directional_beta$Taxa), ] # to use for plots


################################################################################
# Interaction Term Models ======================================================
coefInt <- function(x){
    R2adj <- rep(summary(x)$adj.r.squared, times=1000)
    sum <- summary(x)$coefficients
    p_valDis <- rep(sum[3,4], times=500)
    p_valUndis <- rep(sum[4,4], times=500)
    summa <- as.data.frame(cbind("R2"=R2adj, "pval"=c(p_valDis, p_valUndis))) 
  return(summa)
}

names(alphalist23)

alpha_int_x <- lapply(alphalist23, function(x) {data.frame("time_pred"=c(seq(from = min(x$time_pred), to = max(x$time_pred), length= 500),
                                                                       seq(from = min(x$time_pred), to = max(x$time_pred), length= 500)),
                                                           "Group"=rep(c("Disturbed", "Undisturbed"), each=500))})

sum_alpha_int <- lapply(alpha_int, function(x) {coefInt(x)})
pred_alpha_int <- mapply(function(x,y,z) {data.frame(y, predict(x, y, interval="confidence"),z)},
                         alpha_int, alpha_int_x, sum_alpha_int, SIMPLIFY = F) 
pred_alpha_int <- do.call(rbind, pred_alpha_int)
pred_alpha_int$Taxa <- str_split_fixed(rownames(pred_alpha_int),"\\.", 4)[,2]
pred_alpha_int$CoverageVals <- str_split_fixed(rownames(pred_alpha_int),"\\.", 4)[,3]
pred_alpha_int$Orderq <- str_split_fixed(rownames(pred_alpha_int),"\\.", 4)[,1]
rownames(pred_alpha_int) <- 1:nrow(pred_alpha_int) 
pred_alpha_int$Sig <- ifelse(pred_alpha_int$pval < 0.05, 1, 0)
pred_alpha_int$Taxa[pred_alpha_int$Taxa == "FishTrait"] <- "Fish (Trait)"
pred_alpha_int$Taxa <- factor(pred_alpha_int$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                              "Fish", "Fish (Trait)"),
                                                   labels = c("Diatoms", "Invertebrates", 
                                                              "Fish", "Fish (Trait)"))
pred_alpha_int <- pred_alpha_int[order(pred_alpha_int$Taxa), ] # to use for plots on interaction term models



names(betalist23)
beta_int_x <- lapply(betalist23, function(x) {data.frame("time_pred"=c(seq(from = min(x$time_pred), to = max(x$time_pred), length= 500),
                                                                         seq(from = min(x$time_pred), to = max(x$time_pred), length= 500)),
                                                          "Group"=rep(c("Disturbed", "Undisturbed"), each=500))})

sum_beta_int <- lapply(beta_int, function(x) {coefInt(x)})
pred_beta_int <- mapply(function(x,y,z) {data.frame(y, predict(x, y, interval="confidence"),z)},
                        beta_int, beta_int_x, sum_beta_int, SIMPLIFY = F) 
pred_beta_int <- do.call(rbind, pred_beta_int)
pred_beta_int$Taxa <- str_split_fixed(rownames(pred_beta_int),"\\.", 4)[,2]
pred_beta_int$CoverageVals <- str_split_fixed(rownames(pred_beta_int),"\\.", 4)[,3]
pred_beta_int$Orderq <- str_split_fixed(rownames(pred_beta_int),"\\.", 4)[,1]
rownames(pred_beta_int) <- 1:nrow(pred_beta_int) 
pred_beta_int$Sig <- ifelse(pred_beta_int$pval < 0.05, 1, 0)
pred_beta_int$Taxa[pred_beta_int$Taxa == "FishTrait"] <- "Fish (Trait)"
pred_beta_int$Taxa <- factor(pred_beta_int$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                              "Fish", "Fish (Trait)"),
                              labels = c("Diatoms", "Invertebrates", 
                                         "Fish", "Fish (Trait)"))
pred_beta_int <- pred_beta_int[order(pred_beta_int$Taxa), ]    # to use for plots on interaction term models

#View(pred_directional_alpha)
#View(pred_directional_beta)
#View(pred_alpha_int)
#View(pred_beta_int)

# Interaction Term Plot: -------------------------------------------------------
lPA_Alpha$Taxa[lPA_Alpha$Taxa=="FishTrait"] <- "Fish (Trait)"
lPA_Alpha$Taxa <- factor(lPA_Alpha$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                    "Fish", "Fish (Trait)"),
                                         labels = c("Diatoms", "Invertebrates", 
                                                    "Fish", "Fish (Trait)"))

lBetaLag$Taxa[lBetaLag$Taxa=="FishTrait"] <- "Fish (Trait)"
lBetaLag$Taxa <- factor(lBetaLag$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                  "Fish", "Fish (Trait)"),
                                        labels = c("Diatoms", "Invertebrates", 
                                                   "Fish", "Fish (Trait)"))

plot_int <- function(raw, pred, orderq=NULL, cov=NULL, abline=NULL, type=NULL){
  
  raw_int <- subset(raw, raw$SC2==cov & raw$Order.q==orderq) 
  pred_int <- subset(pred, pred$CoverageVals==cov & pred$Orderq==orderq)
  
  if (type=="alpha"){
    names(raw_int)[names(raw_int)=="Session"] <- "Time"
  }else{
    names(raw_int)[names(raw_int)=="LagSqrt"] <- "Time"
    }
  
  (p <- ggplot(data=raw_int, aes(x=Time, y=Estimate, color=as.factor(Group)))+
      labs(x="Temporal session", y=paste0("q", orderq, " (", cov, ")"), color="Disturbance")+
      geom_point(data=raw_int, aes(x=Time, y=Estimate, color=as.factor(Group)), size=1.5, shape=16, alpha=0.2)+
      scale_color_manual(labels=c("Disturbed", "Undisturbed"), values=c("#CE780F", "#0A7D9D"), breaks =c("Disturbed", "Undisturbed"))+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Group)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free"))
  (p <- p +
      geom_line(data =  pred_int[pred_int$pvalTime < 0.1,], aes(x=time_pred, y = fit, color=as.factor(Group), linetype=as.factor(SigTime)), size=0.8)+
      geom_ribbon(data =  pred_int[pred_int$pvalTime < 0.1,], aes(x=time_pred, y=fit, ymin = lwr, ymax = upr, color=as.factor(Group)), fill="white", alpha = 0.1)+
      #scale_size_manual(labels=c("0", "1"), values = c(1, 1.2), guide = "none", breaks = c(0,1)) +
      scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
      facet_wrap(~ Taxa, scales="free", nrow = 1))
  (p <- p + geom_vline(data = raw_int[raw_int$Taxa %in% c("Fish", "Fish (Trait)"), ], 
                                     aes(xintercept = abline), linetype = "dashed", color = "black"))
  if (type=="alpha"){
    p <- p
  }else{
    p <- p + labs(x="Sqrt time lag")
  }
  
  if (orderq==0){
    p <- p + labs(x="")
  }else{
    p <- p
  }
  
  return(p)
}

# Best: ------------------------------------------------------------------------
Alpha_Cmin_q0 <- plot_int(raw=lPA_Alpha, pred=pred_directional_alpha, orderq=0, cov="Cmin", abline=19.5, type="alpha")
Alpha_Cmin_q2 <- plot_int(raw=lPA_Alpha, pred=pred_directional_alpha, orderq=2, cov="Cmin", abline=19.5, type="alpha")
Alpha_Cmax_q0 <- plot_int(raw=lPA_Alpha, pred=pred_directional_alpha, orderq=0, cov="Cmax", abline=19.5, type="alpha")
Alpha_Cmax_q2 <- plot_int(raw=lPA_Alpha, pred=pred_directional_alpha, orderq=2, cov="Cmax", abline=19.5, type="alpha")
Alpha_Best_Cmin <- Alpha_Cmin_q0/Alpha_Cmin_q2 + plot_layout(guides = "collect", axes = "collect") 
Alpha_Best_Cmax <- Alpha_Cmax_q0/Alpha_Cmax_q2 + plot_layout(guides = "collect", axes = "collect") 


Beta_Cmin_q0 <- plot_int(raw=lBetaLag, pred=pred_directional_beta, orderq=0, cov="Cmin", abline=4.5, type="beta")
Beta_Cmin_q2 <- plot_int(raw=lBetaLag, pred=pred_directional_beta, orderq=2, cov="Cmin", abline=4.5, type="beta")
Beta_Cmax_q0 <- plot_int(raw=lBetaLag, pred=pred_directional_beta, orderq=0, cov="Cmax", abline=4.5, type="beta")
Beta_Cmax_q2 <- plot_int(raw=lBetaLag, pred=pred_directional_beta, orderq=2, cov="Cmax", abline=4.5, type="beta")
Beta_Best_Cmin <- Beta_Cmin_q0/Beta_Cmin_q2  + plot_layout(guides = "collect", axes = "collect") 
Beta_Best_Cmax <- Beta_Cmax_q0/Beta_Cmax_q2 + plot_layout(guides = "collect", axes = "collect") 

ggsave(paste0(plot_path, "/SM_C6_Alpha_Best_Cmin.png"), Alpha_Best_Cmin, width = 10, height = 6.5)
ggsave(paste0(plot_path, "/SM_C6_Alpha_Best_Cmax.png"), Alpha_Best_Cmax, width = 10, height = 6.5)

ggsave(paste0(plot_path, "/SM_C6_Beta_Best_Cmin.png"), Beta_Best_Cmin, width = 10, height = 6.5)
ggsave(paste0(plot_path, "/SM_C6_Beta_Best_Cmax.png"), Beta_Best_Cmax, width = 10, height = 6.5)


# Interaction:------------------------------------------------------------------
# IN function change "pred_int$pvalTime<0.1" to "pred_int$pval<0.1"
# Change also SigTime for Sig
Alpha_Cmin_q0 <- plot_int(raw=lPA_Alpha, pred=pred_alpha_int, orderq=0, cov="Cmin", abline=19.5, type="alpha")
Alpha_Cmin_q2 <- plot_int(raw=lPA_Alpha, pred=pred_alpha_int, orderq=2, cov="Cmin", abline=19.5, type="alpha")
Alpha_Cmax_q0 <- plot_int(raw=lPA_Alpha, pred=pred_alpha_int, orderq=0, cov="Cmax", abline=19.5, type="alpha")
Alpha_Cmax_q2 <- plot_int(raw=lPA_Alpha, pred=pred_alpha_int, orderq=2, cov="Cmax", abline=19.5, type="alpha")
Alpha_Interaction_Cmin <- Alpha_Cmin_q0/Alpha_Cmin_q2 + plot_layout(guides = "collect", axes = "collect") 
Alpha_Interaction_Cmax <- Alpha_Cmax_q0/Alpha_Cmax_q2 + plot_layout(guides = "collect", axes = "collect") 


Beta_Cmin_q0 <- plot_int(raw=lBetaLag, pred=pred_beta_int, orderq=0, cov="Cmin", abline=4.5, type="beta")
Beta_Cmin_q2 <- plot_int(raw=lBetaLag, pred=pred_beta_int, orderq=2, cov="Cmin", abline=4.5, type="beta")
Beta_Cmax_q0 <- plot_int(raw=lBetaLag, pred=pred_beta_int, orderq=0, cov="Cmax", abline=4.5, type="beta")
Beta_Cmax_q2 <- plot_int(raw=lBetaLag, pred=pred_beta_int, orderq=2, cov="Cmax", abline=4.5, type="beta")
Beta_Interaction_Cmin <- Beta_Cmin_q0/Beta_Cmin_q2  + plot_layout(guides = "collect", axes = "collect") 
Beta_Interaction_Cmax <- Beta_Cmax_q0/Beta_Cmax_q2 + plot_layout(guides = "collect", axes = "collect") 

ggsave(paste0(plot_path, "/SM_C6_Alpha_Interaction_Cmin.png"), Alpha_Interaction_Cmin, width = 10, height = 6.5)
ggsave(paste0(plot_path, "/SM_C6_Alpha_Interaction_Cmax.png"), Alpha_Interaction_Cmax, width = 10, height = 6.5)

ggsave(paste0(plot_path, "/SM_C6_Beta_Interaction_Cmin.png"), Beta_Interaction_Cmin, width = 10, height = 6.5)
ggsave(paste0(plot_path, "/SM_C6_Beta_Interaction_Cmax.png"), Beta_Interaction_Cmax, width = 10, height = 6.5)


# End of script ################################################################
################################################################################