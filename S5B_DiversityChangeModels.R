################################################################################
# S4 Models
# AFE
# January 2023
################################################################################


################################################################################
# LIBRARIES: ===================================================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(AICcmodavg)
library(broom)
library(data.table)
library(scales)

rm(list=ls())


################################################################################
# DATA: =========================================================================
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
Results_BetaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_BetaDivRData"

load(paste0(Results_AlphaDivRData, "/S3_lPA_Alpha.RData"))
load(paste0(Results_BetaDivRData, "/S3_lBetaLag.RData"))
load(paste0(Results_BetaDivRData, "/S3_lBetaDU.RData"))


# Alpha: -----------------------------------------------------------------------
lPA_Alpha <- subset(lPA_Alpha, !lPA_Alpha$Order.q==1)
lPA_Alpha22 <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Taxa, lPA_Alpha$Group, lPA_Alpha$SC2)) 
lPA_Alpha22 <- lPA_Alpha22[which(lapply(lPA_Alpha22, nrow) != 0)]                      # for models, although it is 23 after May 2024 :)

lPA_AlphaFish <- lPA_Alpha[(lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)")),]
lPA_Alpha1622 <- lPA_AlphaFish[(lPA_AlphaFish$Session %in% c(23, 46, 50, 54)),]        # sensitivity (second period)
lPA_AlphaDry <- lPA_AlphaFish[(lPA_AlphaFish$Session %in% c(2,6,10,14,18,46,50,54)),]  # sensitivity (dry)

lPA_Alpha1622 <- split(lPA_Alpha1622, list(lPA_Alpha1622$Order.q, lPA_Alpha1622$Taxa, lPA_Alpha1622$Group, lPA_Alpha1622$SC2)) 
lPA_Alpha1622 <- lPA_Alpha1622[which(lapply(lPA_Alpha1622, nrow) != 0)]                # sensitivity (second period)

lPA_AlphaDry <- split(lPA_AlphaDry, list(lPA_AlphaDry$Order.q, lPA_AlphaDry$Taxa, lPA_AlphaDry$Group, lPA_AlphaDry$SC2)) 
lPA_AlphaDry <- lPA_AlphaDry[which(lapply(lPA_AlphaDry, nrow) != 0)]                   # sensitivity (dry)

lPA_Alpha <- lPA_Alpha[!(lPA_Alpha$Session %in% c(23, 46, 50, 54)),]
lPA_Alpha <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Taxa, lPA_Alpha$Group, lPA_Alpha$SC2)) 
lPA_Alpha <- lPA_Alpha[which(lapply(lPA_Alpha, nrow) != 0)]                           # to use for models (and linear fish sensitivity period 1)

ldtA22 <- lPA_Alpha22    # main
ldtA <- lPA_Alpha        # main & main fish nonlinear
ldt1622 <- lPA_Alpha1622 # sensitivity fish 2nd period
ldtdry <- lPA_AlphaDry   # sensitivity fish dry


# Beta Lag: --------------------------------------------------------------------
names(lBetaLag)
lBetaLag <- as.data.frame(lBetaLag)
str(lBetaLag)

# Check obs that might be biased:
sum(lBetaLag$Estimate < 0) # 333
333/15264
table(lBetaLag$Order.q[lBetaLag$Estimate<0])   
table(lBetaLag$Group[lBetaLag$Estimate<0])     
table(lBetaLag$Dataset[lBetaLag$Estimate<0])   

lBetaLag <- subset(lBetaLag, !lBetaLag$Order.q==1)

lBetaLagFish <- lBetaLag[(lBetaLag$Taxa %in% c("Fish", "Fish (Trait)")),]
lBetaLag1622 <- lBetaLagFish[(lBetaLagFish$s1 %in% c(23, 46, 50, 54)),]        # sensitivity (second period)
lBetaLag1622 <- lBetaLag1622[(lBetaLag1622$s2 %in% c(23, 46, 50, 54)),]        # sensitivity (second period)

lBetaLagDry <- lBetaLagFish[(lBetaLagFish$s1 %in% c(2,6,10,14,18,46,50,54)),]  # sensitivity (dry)
lBetaLagDry <- lBetaLagDry[(lBetaLagDry$s2 %in% c(2,6,10,14,18,46,50,54)),]    # sensitivity (dry)

lBetaLag1622 <- split(lBetaLag1622, list(lBetaLag1622$Order.q, lBetaLag1622$Taxa, lBetaLag1622$Group, lBetaLag1622$SC2)) 
lBetaLag1622 <- lBetaLag1622[which(lapply(lBetaLag1622, nrow) != 0)]                # sensitivity (second period)

lBetaLagDry <- split(lBetaLagDry, list(lBetaLagDry$Order.q, lBetaLagDry$Taxa, lBetaLagDry$Group, lBetaLagDry$SC2)) 
lBetaLagDry <- lBetaLagDry[which(lapply(lBetaLagDry, nrow) != 0)]                   # sensitivity (dry)

lBetaLag22 <- split(lBetaLag, list(lBetaLag$Order.q, lBetaLag$Taxa, lBetaLag$Group, lBetaLag$SC2)) 
lBetaLag22 <- lBetaLag22[which(lapply(lBetaLag22, nrow) != 0)]                      # to use for models

lBetaLag <- lBetaLag[!(lBetaLag$s1 %in% c(23, 46, 50, 54)),]
lBetaLag <- lBetaLag[!(lBetaLag$s2 %in% c(23, 46, 50, 54)),]
range(lBetaLag$s1)
range(lBetaLag$Lag)
lBetaLag <- split(lBetaLag, list(lBetaLag$Order.q, lBetaLag$Taxa, lBetaLag$Group, lBetaLag$SC2)) 
lBetaLag <- lBetaLag[which(lapply(lBetaLag, nrow) != 0)]   

ldtBLag22 <- lBetaLag22     # main
ldtBLag <- lBetaLag         # main & main fish nonlinear
ldtBLag1622 <- lBetaLag1622 # sensitivity fish 2nd period
ldtBLagDry <- lBetaLagDry   # sensitivity fish dry


# Beta DU: ---------------------------------------------------------------------
names(lBetaDU)
lBetaDU <- as.data.frame(lBetaDU)
lBetaDU <- subset(lBetaDU, !lBetaDU$Order.q==1)
str(lBetaDU)
names(lBetaDU)[names(lBetaDU)=="s1"] <- "Session"

lBetaDU22 <- split(lBetaDU, list(lBetaDU$Order.q, lBetaDU$Taxa, lBetaDU$Group, lBetaDU$SC2)) 
lBetaDU22 <- lBetaDU22[which(lapply(lBetaDU22, nrow) != 0)]             # to use for models

lBetaDU <- lBetaDU[!(lBetaDU$Session %in% c(23, 46, 50, 54)),]
lBetaDU <- lBetaDU[!(lBetaDU$s2 %in% c(23, 46, 50, 54)),]
range(lBetaDU$Session)
lBetaDU <- split(lBetaDU, list(lBetaDU$Order.q, lBetaDU$Taxa, lBetaDU$Group, lBetaDU$SC2)) 
lBetaDU <- lBetaDU[which(lapply(lBetaDU, nrow) != 0)] 

ldtBDU22 <- lBetaDU22
ldtBDU <- lBetaDU


################################################################################
# POLYNOMIALS ==================================================================

polycoefs <- function(x, var=NULL, return=NULL){

  poly_order <- c(1:4)
  mcoefs <- list()
  ms <- list()
  msn <- list()
  
  for (i in 1:length(poly_order)){
    if(var=="s1"){
      m <- lm(x$Estimate ~ poly(x$Session, i), data=x)
    }else{
      m <- lm(x$Estimate ~ poly(x$LagSqrt, i), data=x)
    }
    
    msum <- summary(m)
    adjR2 <- msum$adj.r.squared
    mse <- mean(msum$residuals^2)
    fstat <- glance(m)$statistic
    pval <- glance(m)$p.value  # for overall p-val
    ms[[i]] <- m
    msn[[i]] <- paste0("poly","_", i)
    mcoefs[[i]] <- data.frame("AdjR2"=adjR2,"MSE"=mse, "FStatistic"=fstat, "Pval"=pval)
  }
  
  msn <- as.vector(do.call(rbind, msn))
  aics <- as.data.frame(aictab(ms, modnames=msn))
  
  mcoefs <- as.data.frame(do.call(rbind, mcoefs))
  mcoefs$Modnames <- paste0("poly","_", c(1:4))
  
  mcoefsII <- merge(aics, mcoefs, by="Modnames")
  mcoefsII <- mcoefsII[order(mcoefsII$Delta_AICc),]
  
  if(return=="coefs"){
    return(mcoefsII)
  }else{
    return(ms)
    }
  
} 

mpolycoefsAlpha <- lapply(ldtA, function(x) {polycoefs(x, var="s1", return = "coefs")})
mpolycoefsAlpha22 <- lapply(ldtA22, function(x) {polycoefs(x, var="s1", return = "coefs")})
mpolycoefsBetaLag <- lapply(ldtBLag, function(x) {polycoefs(x, var="Lag", return = "coefs")})
mpolycoefsBetaLag22 <- lapply(ldtBLag22, function(x) {polycoefs(x, var="Lag", return = "coefs")})
mpolycoefsBetaDU <- lapply(ldtBDU, function(x) {polycoefs(x, var="s1", return = "coefs")})
mpolycoefsBetaDU22 <- lapply(ldtBDU22, function(x) {polycoefs(x, var="s1", return = "coefs")})


################################################################################
# POLY SELECTION: ==============================================================

table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"

orgC6 <- function(x){
  x <- c(
    x[grep("Disturbed", names(x))],
    x[grep("Undisturbed", names(x))] # and already removes NR
  )
  return(x)
} # for C6
orgDU <- function(x){
  x <- c(
    x[grep("Diatoms.DU.Cmin", names(x))],
    x[grep("Diatoms.DU.Cmax", names(x))],
    x[grep("Invertebrates.DU.Cmin", names(x))],
    x[grep("Invertebrates.DU.Cmax", names(x))],
    x[grep("Fish.DU.Cmin", names(x))],
    x[grep("Fish.DU.Cmax", names(x))],
    x[grep("FishTrait.DU.Cmin", names(x))],
    x[grep("FishTrait.DU.Cmax", names(x))]
  )
  return(x)
} # for C6
orgNR <- function(x){
  x <- c(
    x[grep("Diatoms.Northern Range.Cmin", names(x))],
    x[grep("Diatoms.Northern Range.Cmax", names(x))],
    x[grep("Invertebrates.Northern Range.Cmin", names(x))],
    x[grep("Invertebrates.Northern Range.Cmax", names(x))],
    x[grep("Fish.Northern Range.Cmin", names(x))],
    x[grep("Fish.Northern Range.Cmax", names(x))],
    x[grep("Fish \\(Trait\\).Northern Range.Cmin", names(x))],
    x[grep("Fish \\(Trait\\).Northern Range.Cmax", names(x))]
    )
  return(x)
} # for C5

rm <- c("K","ModelLik","AICcWt", "Cum.Wt", "MSE") 
model_sel <- function(dt){
  dt <- do.call(rbind, dt)
  dt <- dt[, ! names(dt) %in% rm]
  dt$Code <- str_replace(rownames(dt), "\\.[^.]*$", "")
  dt <- separate(dt, Code, into = paste0("Column", 1:4), sep = "\\.") # warning OK
  dt$Modnames <- gsub("poly_", "", dt$Modnames)
  dt <- dt %>% relocate(c(Column3, Column2,Column4, Column1), .before = Modnames)
  dt$Column2 <- factor(dt$Column2, levels = c("Diatoms", "Invertebrates", 
                                              "Fish", "Fish (Trait)"))                
  dt <- dt[order(dt$Column2), ]                                            
  dt <- dt %>% arrange(Column2, desc(Column4), Column1) 
  return(dt)
}

# Chapter 5: -------------------------------------------------------------------
mpolycoefsAlphaC5 <- mpolycoefsAlpha[grep("Northern Range", names(mpolycoefsAlpha))]
mpolycoefsBetaC5 <- mpolycoefsBetaLag[grep("Northern Range", names(mpolycoefsBetaLag))]

aicAlphaC5 <- model_sel(mpolycoefsAlphaC5)
aicBetaC5 <- model_sel(mpolycoefsBetaC5)

#write.csv(aicAlphaC5, file=paste0(table_path, "/Chapter5/S5B_C5_Coef_NonlinearAlpha.csv"), row.names = F)
#write.csv(aicBetaC5, file=paste0(table_path, "/Chapter5/S5B_C5_Coef_NonlinearBeta.csv"), row.names = F)

# Chapter 6: -------------------------------------------------------------------
mpolycoefsAlphaC6 <- orgC6(mpolycoefsAlpha)
mpolycoefsBetaC6 <- orgC6(mpolycoefsBetaLag)

aicAlphaC6 <- model_sel(mpolycoefsAlphaC6)
aicBetaC6 <- model_sel(mpolycoefsBetaC6)

#write.csv(aicAlphaC6, file=paste0(table_path, "/Chapter6/S5B_C6_AIC_NonlinearAlpha.csv"), row.names = F)
#write.csv(aicBetaC6, file=paste0(table_path, "/Chapter6/S5B_C6_AIC_NonlinearBeta.csv"), row.names = F)


################################################################################
# COEFS LINEAR =================================================================
table_coefs <- function(data, type=NULL){
  if(type=="alpha"){
    models <- lapply(data, function(x) {lm(x$Estimate ~ poly(x$Session, 1), data=x)})
  }else{
    models <- lapply(data, function(x) {lm(x$Estimate ~ poly(x$LagSqrt, 1), data=x)})
  }
  
  Table <- lapply(models , function(x) {as.data.frame(summary(x)$coefficients)})
  Table <- do.call(rbind, Table)
  Table$Code <- str_replace(rownames(Table), "\\.[^.]*$", "")
  Table <- Table[!rownames(Table) %like%"(Intercept)",]
  
  R2 <- lapply(models, function(x) {as.data.frame(summary(x)$adj.r.squared)})
  R2 <- do.call(rbind, R2)
  R2$Code <- rownames(R2)
  
  Table$adjR2 <- R2[,1][match(Table$Code, R2$Code)]
  Table <- separate(Table, Code, into = paste0("Column", 1:4), sep = "\\.") # warning OK
  Table <- Table[,!names(Table) == "Column3"]
  names(Table)[names(Table) == "t value"] <- "tvalue"
  Table <- Table %>% relocate(c(Column2,Column4, Column1), .before = Estimate)
  Table <- Table %>% relocate(adjR2, .after = tvalue)
  Table$Column2 <- factor(Table$Column2, levels = c("Diatoms", "Invertebrates", 
                                                    "Fish", "Fish (Trait)"))                
  Table <- Table[order(Table$Column2), ]                                            
  Table <- Table %>% arrange(Column2, desc(Column4), Column1) 
  rownames(Table) <- 1:nrow(Table)
  return(Table)
}

# Chapter 5 (alpha): -----------------------------------------------------------
cl_alpha <- ldtA22[grep("Northern Range", names(ldtA22))]
#lapply(cl_alpha, function(x) {is.numeric(x$Session)}) 
cl_alpha <- table_coefs(cl_alpha, type = "alpha")                  # main (dia & inv up to 19, fish up to 23 sampling events)

cl_alpha_sen1 <- ldtA[grep("Northern Range", names(ldtA))]
#lapply(cl_alpha_sen1, function(x) {is.numeric(x$Session)}) 
cl_alpha_sen1 <- cl_alpha_sen1[grep("Fish", names(cl_alpha_sen1))]
cl_alpha_sen1 <- table_coefs(cl_alpha_sen1, type = "alpha")        # sensitivity Fish 1st period

cl_alpha_sen2 <- ldt1622[grep("Northern Range", names(ldt1622))]
#lapply(cl_alpha_sen2, function(x) {is.numeric(x$Session)}) 
cl_alpha_sen2 <- cl_alpha_sen2[grep("Fish", names(cl_alpha_sen2))]
cl_alpha_sen2 <- table_coefs(cl_alpha_sen2, type = "alpha")        # sensitivity Fish 2nd period

cl_alpha_sen3 <- ldtdry[grep("Northern Range", names(ldtdry))]
#lapply(cl_alpha_sen3, function(x) {is.numeric(x$Session)}) 
cl_alpha_sen3 <- cl_alpha_sen3[grep("Fish", names(cl_alpha_sen3))]
cl_alpha_sen3 <- table_coefs(cl_alpha_sen3, type = "alpha")        # sensitivity Fish end dry season sampling events 


#write.csv(cl_alpha, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_Alpha.csv"), row.names = F)
#write.csv(cl_alpha_sen1, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_AlphaSenFish1Period.csv"), row.names = F)
#write.csv(cl_alpha_sen2, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_AlphaSenFish2Period.csv"), row.names = F)
#write.csv(cl_alpha_sen3, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_AlphaSenFish3Dry.csv"), row.names = F)


# Chapter 5 (beta): ------------------------------------------------------------
cl_betalag <- ldtBLag22[grep("Northern Range", names(ldtBLag22))]
#lapply(cl_betalag, function(x) {is.numeric(x$LagSqrt)}) 
cl_betalag <- table_coefs(cl_betalag, type = "beta")

cl_betalag_sen1 <- ldtBLag[grep("Northern Range", names(ldtBLag))]
cl_betalag_sen1 <- cl_betalag_sen1[grep("Fish", names(cl_betalag_sen1))]
cl_betalag_sen1 <- table_coefs(cl_betalag_sen1, type="beta")     # sensitivity Fish 1st period

cl_betalag_sen2 <- ldtBLag1622[grep("Northern Range", names(ldtBLag1622))]
cl_betalag_sen2 <- cl_betalag_sen2[grep("Fish", names(cl_betalag_sen2))]
cl_betalag_sen2 <- table_coefs(cl_betalag_sen2, type="beta")     # sensitivity Fish 2nd period

cl_betalag_sen3 <- ldtBLagDry[grep("Northern Range", names(ldtBLagDry))]
cl_betalag_sen3 <- cl_betalag_sen3[grep("Fish", names(cl_betalag_sen3))]
cl_betalag_sen3 <- table_coefs(cl_betalag_sen3, type="beta")     # sensitivity Fish end dry season sampling events 


#write.csv(cl_betalag, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_Beta.csv"), row.names = F)
#write.csv(cl_betalag_sen1, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_BetaSenFish1Period.csv"), row.names = F)
#write.csv(cl_betalag_sen2, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_BetaSenFish2Period.csv"), row.names = F)
#write.csv(cl_betalag_sen3, file=paste0(table_path, "/Chapter5/S5B_C5_Linear_CoefTable_BetaSenFish3Dry.csv"), row.names = F)


# Linear DU: -------------------------------------------------------------------
cl_du <- orgDU(ldtBDU22)
cl_du <- table_coefs(cl_du, type = "alpha") # beta but with alpha structure
cl_du_sen <- orgDU(ldtBDU)
cl_du_sen <- cl_du_sen[grep("Fish", names(cl_du_sen))]
cl_du_sen <- table_coefs(cl_du_sen, type = "alpha") # beta but with alpha structure

#write.csv(cl_du, file=paste0(table_path, "/Chapter6/S5B_C6_DU_CoefTable.csv"), row.names = F)
#write.csv(cl_du_sen, file=paste0(table_path, "/Chapter6/S5B_C6_DU_CoefTableSen.csv"), row.names = F)


################################################################################
# BEST POLYNOMIAL: =============================================================

# 1) AICc & R2: ----------------------------------------------------------------
bestpolyF <- function(mpolycoefs){
bestpoly <- list()
for (i in 1:length(mpolycoefs)){
  
  x <- mpolycoefs[[i]]
  
  a <- min(x$AICc)
  b <- min(x$AICc) + 2
  
  y <- x %>% filter(., between(AICc, a, b) )
  
  if (nrow(y)==1){
    y <- y
  }
  if (nrow(y)>1){
    y <- y %>% filter(., AdjR2==max(AdjR2))
  }
  
  bp <- y[1,1]
  bp <- gsub("poly_", "", bp)
  
  bestpoly[[i]] <- bp
}
return(bestpoly)
} 

bestpolyAlpha <- bestpolyF(mpolycoefsAlpha)
bestpolyAlpha22 <- bestpolyF(mpolycoefsAlpha22)
bestpolyBetaLag <- bestpolyF(mpolycoefsBetaLag)
bestpolyBetaLag22 <- bestpolyF(mpolycoefsBetaLag22)
bestpolyBetaDU <- bestpolyF(mpolycoefsBetaDU)
bestpolyBetaDU22 <- bestpolyF(mpolycoefsBetaDU22)


################################################################################
# MODEL VALIDATIONS: ===========================================================
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots/Chapter5"


ldtA_C5 <- orgNR(ldtA)
ldtA22_C5 <- orgNR(ldtA22)

alpha_models <- lapply(ldtA_C5, function(x) {polycoefs(x, var="s1", return = "models")})
alpha_models22 <- lapply(ldtA22_C5, function(x) {polycoefs(x, var="s1", return = "models")})

ldtBLag_C5 <- orgNR(ldtBLag)
ldtBLag22_C5 <- orgNR(ldtBLag22)

beta_models <- lapply(ldtBLag_C5, function(x) {polycoefs(x, var="Lag", return = "models")})
beta_models22 <- lapply(ldtBLag22_C5, function(x) {polycoefs(x, var="Lag", return = "models")})


# Chapter 5 (directional): -----------------------------------------------------
alpha_linear <- lapply(alpha_models22, function(x) {x[[1]]})
alpha_linear_sresid <- lapply(alpha_linear, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)}) # standardised residuals (also alpha int)
namesAlpha <- names(alpha_linear_sresid)
namesAlpha <- gsub("\\.", " ", namesAlpha)
namesAlpha <- gsub("Northern Range", "NR", namesAlpha)

beta_linear <- lapply(beta_models22, function(x) {x[[1]]})
beta_linear_sresid <- lapply(beta_linear, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)}) # standardised residuals (also alpha int)
namesBeta <- names(beta_linear_sresid)
namesBeta <- gsub("\\.", " ", namesBeta)
namesBeta <- gsub("Northern Range", "NR", namesBeta)


#jpeg(file=paste0(plot_path,"/S5B_C5_Alpha_Linear_Normality.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#  }, alpha_linear_sresid, namesAlpha, SIMPLIFY = T) 
#dev.off()  # normality of residuals
#jpeg(file=paste0(plot_path,"/S5B_C5_Beta_Linear_Normality.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#}, beta_linear_sresid, namesBeta, SIMPLIFY = T) 
#dev.off()  # normality of residuals


#jpeg(file=paste0(plot_path,"/S5B_C5_Alpha_Linear_Homoscedasticity.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
# plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, alpha_linear_sresid, alpha_linear, namesAlpha, SIMPLIFY = T)
#dev.off()  # homoscedasticity
#jpeg(file=paste0(plot_path,"/S5B_C5_Beta_Linear_Homoscedasticity.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
#  plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, beta_linear_sresid, beta_linear, namesBeta, SIMPLIFY = T)
#dev.off()  # homoscedasticity


# Chapter 5 (nonlinear): -------------------------------------------------------
bestA_C5 <- orgNR(mpolycoefsAlpha)
bestA_C5 <- bestpolyF(bestA_C5)
bestA_C5 <- lapply(bestA_C5, as.numeric)

alpha_poly <- mapply(function(mod,best) {mod[[best]]}, alpha_models, bestA_C5, SIMPLIFY=F)
names(alpha_poly) <- paste0(lapply(names(alpha_poly), function(x) {x}), ".", "poly.",
                                   lapply(alpha_poly, function(x) {x$rank - 1 }))
alpha_poly_sresid <- lapply(alpha_poly, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)}) # standardised residuals (also alpha int)

namesAlphaPoly <- names(alpha_poly_sresid)
namesAlphaPoly <- gsub("\\.", " ", namesAlphaPoly)
namesAlphaPoly <- gsub("Northern Range", "NR", namesAlphaPoly)


bestB_C5 <- orgNR(mpolycoefsBetaLag)
bestB_C5 <- bestpolyF(bestB_C5)
bestB_C5 <- lapply(bestB_C5, as.numeric)

beta_poly <- mapply(function(mod,best) {mod[[best]]}, beta_models, bestB_C5, SIMPLIFY=F)
names(beta_poly) <- paste0(lapply(names(beta_poly), function(x) {x}), ".", "poly.",
                            lapply(beta_poly, function(x) {x$rank - 1 }))
beta_poly_sresid <- lapply(beta_poly, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)}) # standardised residuals (also alpha int)
namesBetaPoly <- names(beta_poly_sresid)
namesBetaPoly <- gsub("\\.", " ", namesBetaPoly)
namesBetaPoly <- gsub("Northern Range", "NR", namesBetaPoly)


#jpeg(file=paste0(plot_path,"/S5B_C5_Alpha_Poly_Normality.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#  }, alpha_poly_sresid, namesAlphaPoly, SIMPLIFY = T) 
#dev.off()  # normality of residuals
#jpeg(file=paste0(plot_path,"/S5B_C5_Beta_Poly_Normality.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#}, beta_poly_sresid, namesBetaPoly, SIMPLIFY = T) 
#dev.off()  # normality of residuals


#jpeg(file=paste0(plot_path,"/S5B_C5_Alpha_Poly_Homoscedasticity.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
# plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, alpha_poly_sresid, alpha_poly, namesAlphaPoly, SIMPLIFY = T)
#dev.off()  # homoscedasticity
#jpeg(file=paste0(plot_path,"/S5B_C5_Beta_Poly_Homoscedasticity.jpg"), width = 900, height = 900)
#par(mfrow=c(4,4))
#mapply(function(x,y,z){
#  plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, beta_poly_sresid, beta_poly, namesBetaPoly, SIMPLIFY = T)
#dev.off()  # homoscedasticity


# Chapter 6 (nonlinear): -------------------------------------------------------
# TBC


################################################################################
# PREDICT: =====================================================================
polypredict <- function(x, y, var = NULL){
  
  y <- as.numeric(y)
  
  if(var=="s1")
  {regionx <- data.frame("Session"=seq(from = min(x$Session), to = max(x$Session), length= 1000))}else{
   regionx <- data.frame("LagSqrt"=seq(from = min(x$LagSqrt), to = max(x$LagSqrt), length= 1000)) 
  }
  
  labelsx <- data.frame("Order.q" = rep(unique(x$Order.q), times=1000),
                        "SC2" = rep(unique(x$SC2), times=1000),
                        "Group"=rep(unique(x$Group), times=1000),
                        "Taxa"=rep(unique(x$Taxa), times=1000),
                        "Diversity"=rep(unique(x$Diversity), times=1000),
                        "Type"=rep(unique(x$Type), times=1000))
  if(var=="s1"){
    poly1 <- lm(Estimate ~ poly(Session, 1), data=x)
    polyy <- lm(Estimate ~ poly(Session, y), data=x) 
  }else{
    poly1 <- lm(Estimate ~ poly(LagSqrt, 1), data=x)
    polyy <- lm(Estimate ~ poly(LagSqrt, y), data=x)
  }
  
  
  lR2adj <- summary(poly1)$adj.r.squared
  yR2adj <- summary(polyy)$adj.r.squared
  
  lpval <- glance(poly1)$p.value  
  ypval <- glance(polyy)$p.value  
  
  lpred <- data.frame(labelsx, regionx, predict(poly1, regionx, interval = "confidence"))
  ypred <- data.frame(labelsx, regionx, predict(polyy, regionx, interval = "confidence"))
  
  if(var=="s1"){
  pred <- merge(x = lpred, y = ypred, by = c("Session", "Order.q", "SC2", "Group","Taxa", "Diversity", "Type"))}else{
  pred <- merge(x = lpred, y = ypred, by = c("LagSqrt", "Order.q", "SC2", "Group","Taxa", "Diversity", "Type")) 
  }
  pred <- data.frame(pred, "lR2adj"=rep(lR2adj, nrow(pred)), "yR2adj"=rep(yR2adj, nrow(pred)), "lpval"=rep(lpval, nrow(pred)), "ypval"=rep(ypval, nrow(pred)))
  
  return(pred)
  
}

mpolypredAlpha <- mapply(function(x, y) polypredict(x, y, var="s1"), ldtA, bestpolyAlpha, SIMPLIFY = F)
mpolypredAlpha22 <- mapply(function(x, y) polypredict(x, y, var="s1"), ldtA22, bestpolyAlpha22, SIMPLIFY = F)
mpolypredBetaLag <- mapply(function(x, y) polypredict(x, y, var="Lag"), ldtBLag, bestpolyBetaLag, SIMPLIFY = F)
mpolypredBetaLag22 <- mapply(function(x, y) polypredict(x, y, var="Lag"), ldtBLag22, bestpolyBetaLag22, SIMPLIFY = F)
mpolypredBetaDU <- mapply(function(x, y) polypredict(x, y, var="s1"), ldtBDU, bestpolyBetaDU, SIMPLIFY = F)
mpolypredBetaDU22 <- mapply(function(x, y) polypredict(x, y, var="s1"), ldtBDU22, bestpolyBetaDU22, SIMPLIFY = F)

funPVals <- function(ldt0){
  ldt <- lapply(ldt0, function(x) {x$same <- rep(NA, nrow(x)); x
  x$same <- ifelse(x$fit.x==x$fit.y, "TRUE", "FALSE"); x}) 
  ldt <- lapply(ldt, function(x) {x$lpvalS <- ifelse(x$lpval < 0.05, 1, 0);x})
  ldt <- lapply(ldt, function(x) {x$ypvalS <- ifelse(x$ypval < 0.05, 1, 0);x})
  return(ldt)
}

mpolypredAlpha <- funPVals(mpolypredAlpha)
mpolypredAlpha <- do.call(rbind, mpolypredAlpha)

mpolypredAlpha22 <- funPVals(mpolypredAlpha22)
mpolypredAlpha22 <- do.call(rbind, mpolypredAlpha22)

mpolypredBetaLag <- funPVals(mpolypredBetaLag)
mpolypredBetaLag<- do.call(rbind, mpolypredBetaLag)
mpolypredBetaLag22 <- funPVals(mpolypredBetaLag22)
mpolypredBetaLag22 <- do.call(rbind, mpolypredBetaLag22)
mpolypredBetaDU <- funPVals(mpolypredBetaDU)
mpolypredBetaDU <- do.call(rbind, mpolypredBetaDU)
mpolypredBetaDU22 <- funPVals(mpolypredBetaDU22)
mpolypredBetaDU22 <- do.call(rbind, mpolypredBetaDU22)

mrawAlpha <- do.call(rbind, ldtA)
mrawAlpha22 <- do.call(rbind, ldtA22)
mrawBetaLag<- do.call(rbind, ldtBLag)
mrawBetaLag22 <- do.call(rbind, ldtBLag22)
mrawBetaDU <- do.call(rbind, ldtBDU)
mrawBetaDU22 <- do.call(rbind, ldtBDU22)

ResultsModels <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/ResultsModels"
save(mpolypredAlpha, mpolypredAlpha22, mpolypredBetaLag, mpolypredBetaLag22, mpolypredBetaDU, mpolypredBetaDU22, file=paste0(ResultsModels,"/mpolypred.RData"))
save(mrawAlpha, mrawAlpha22, mrawBetaLag, mrawBetaLag22, mrawBetaDU, mrawBetaDU22, file=paste0(ResultsModels,"/mraw.RData"))


# End of script ################################################################
################################################################################