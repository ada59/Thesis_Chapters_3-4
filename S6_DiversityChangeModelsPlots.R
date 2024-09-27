################################################################################
# S6 Plots (C5 & C6)
# AFE
# January 2023
################################################################################


################################################################################
# LIBRARIES: ===================================================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
library(patchwork)
library(Hmisc)


################################################################################
# PALETTES: ====================================================================
#https://r-charts.com/color-palette-generator/?utm_content=cmp-true


# Diatoms:
c("#d4ff9f", "#aacc7f", "#7f995f", "#556640", "#2a3320", "#000000")
# Invertebrates:
c("#ffc906", "#cca105", "#997904", "#665002", "#332801", "#000000")
# Fish:
c("#22f4ff", "#1bc3cc", "#149299", "#0e6266", "#073133", "#000000")
# Fish (Trait):
c("#22f4ff", "#1bc3cc", "#149299", "#0e6266", "#073133", "#000000")

c("#bed1d2", "#a3b7c0", "#889eae", "#6d849c", "#516a89", "#365177", "#1b3765")
c("#1875c1", "#3181aa", "#4b8e93", "#649a7c", "#7da666", "#96b24f", "#b0bf38", "#c9cb21")

c("#982a11", "#a5452f", "#b25f4d", "#bf7a6a", "#bf7a6a", "#cc9588", "#d8afa6")
c("#9b8346", "#a8935d", "#b4a274", "#c1b28b")


rm(list=ls())
getwd()


################################################################################
# DATA =========================================================================
ResultsModels <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/ResultsModels"
load(paste0(ResultsModels,"/mpolypred.RData"))
load(paste0(ResultsModels,"/mraw.RData"))

for_plots <- function(dt){
  dt <- dt[dt$Order.q %in% c(0,2),]
  dt$Order.q <- as.character(dt$Order.q)
  dt$Taxa <- factor(dt$Taxa, c("Diatoms"="Diatoms",
                               "Invertebrates"="Invertebrates",
                                "Fish"="Fish",
                                "Fish (Trait)"="Fish (Trait)")) # sort order
  #dt$Order.q <- plyr::revalue(dt$Order.q, c("0"="Order q0",
  #                                          "2"="Order q2"))    # sort order
  return(dt)
}

mrawAlpha <- for_plots(mrawAlpha)
mrawAlpha22 <- for_plots(mrawAlpha22)
mrawBetaLag <- for_plots(mrawBetaLag)
mrawBetaLag22 <- for_plots(mrawBetaLag22)
mrawBetaDU <- for_plots(mrawBetaDU)
mrawBetaDU22 <- for_plots(mrawBetaDU22)

mpolypredAlpha <- for_plots(mpolypredAlpha)
mpolypredAlpha22 <- for_plots(mpolypredAlpha22)
mpolypredBetaLag <- for_plots(mpolypredBetaLag)
mpolypredBetaLag22 <- for_plots(mpolypredBetaLag22) 
mpolypredBetaDU <- for_plots(mpolypredBetaDU)
mpolypredBetaDU22 <- for_plots(mpolypredBetaDU22)


################################################################################
# PLOTS: =======================================================================
PlotsPath <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"


# Chapter 5: ===================================================================

plot_α <- function(raw, raw22, poly, poly22, group=NULL, sc2=NULL){
  
  rawdata <- subset(raw, raw$Group==group & raw$SC2==sc2)
  rawdata22 <-subset(raw22, raw22$Group==group & raw22$SC2==sc2)
  
  model <- subset(poly, poly$Group==group & poly$SC2==sc2)
  model22 <- subset(poly22, poly22$Group==group & poly22$SC2==sc2)
  
  
  p <- ggplot(data=rawdata, aes(x=Session, y=Estimate, color=as.factor(Order.q)))+
    labs(x="", y=paste0("α diversity (coverage = ", sc2, ")"), color="Order.q")+
    geom_point(data=rawdata, aes(x=Session, y=Estimate, color=as.factor(Order.q)), size=1.5, shape=16, alpha=0.5)+
    scale_color_manual(labels=c("Order q0", "Order q2"), values=c("#889eae", "#1b3765"))+
    #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
    #              position=position_dodge(0.05), alpha=0.7, color="grey")+
    theme_classic()+
    facet_wrap(~ Taxa, scales = "free", nrow = 1)
  p_22 <- ggplot(data=rawdata22, aes(x=Session, y=Estimate, color=as.factor(Order.q)))+
    labs(x="Temporal session", y=paste0("α diversity (coverage = ", sc2, ")"), color="Order.q")+
    geom_point(data=rawdata22, aes(x=Session, y=Estimate, color=as.factor(Order.q)), size=1.5, shape=16, alpha=0.5)+
    scale_color_manual(labels=c("Order q0", "Order q2"), values=c("#889eae", "#1b3765"))+
    #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
    #              position=position_dodge(0.05), alpha=0.7, color="grey")+
    theme_classic()+
    facet_wrap(~ Taxa, scales = "free", nrow = 1)
  p_22 <- p_22 + geom_vline(data = rawdata22[rawdata22$Taxa %in% c("Fish", "Fish (Trait)"), ], 
                            aes(xintercept = 19.5), linetype = "dashed", color = "black")
  (p_nonlinear <- p +
      geom_line(data =  model[model$ypval < 0.1,], aes(x=Session, y = fit.y, color=as.factor(Order.q), linetype=as.factor(ypvalS)), size=0.8)+
      geom_ribbon(data =  model[model$ypval < 0.1,], aes(x=Session, y=fit.y, ymin = lwr.y, ymax = upr.y, color=as.factor(Order.q)), fill="white", alpha = 0.1)+
      scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
      facet_wrap(~ Taxa, scales="free", nrow = 1))  
  (p_linear <- p_22 +
      geom_line(data =  model22[model22$lpval < 0.1,], aes(x=Session, y = fit.x, color=as.factor(Order.q), linetype=as.factor(lpvalS)), size=0.8)+
      geom_ribbon(data =  model22[model22$lpval < 0.1,], aes(x=Session, y=fit.x, ymin = lwr.x, ymax = upr.x, color=as.factor(Order.q)), fill="white", alpha = 0.1)+
      scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
      facet_wrap(~ Taxa, scales="free", nrow = 1)) 
  
  p_mix <- p_linear/p_nonlinear + plot_layout(guides = "collect", axes = "collect")
  return(p_mix)
}

(main_Alpha_C5 <- plot_α(mrawAlpha, mrawAlpha22, mpolypredAlpha, mpolypredAlpha22, group="Northern Range", sc2="Cmin"))
(sm_Alpha_C5 <- plot_α(mrawAlpha, mrawAlpha22, mpolypredAlpha, mpolypredAlpha22, group="Northern Range", sc2="Cmax"))


ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_TemporalChange_Alpha_Cmin.png"), main_Alpha_C5, width = 10, height = 6.5) # fix legend aesthetics issue on PowerPoint
ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_TemporalChange_Alpha_Cmax.png"), sm_Alpha_C5, width = 10, height = 6.5)   # fix legend aesthetics issue on PowerPoint



plot_β <- function(raw, poly, group=NULL, sc2=NULL, polynomial=NULL){
  
  rawdata_0 <- subset(raw, raw$Group==group & raw$SC2==sc2 & raw$Order.q=="0")
  rawdata_2 <- subset(raw, raw$Group==group & raw$SC2==sc2 & raw$Order.q=="2")
  model_0 <- subset(poly, poly$Group==group & poly$SC2==sc2 & poly$Order.q=="0")
  model_2 <- subset(poly, poly$Group==group & poly$SC2==sc2 & poly$Order.q=="2")
  
  
  
  if(polynomial=="linear"){
    p_0 <- ggplot(data=rawdata_0, aes(x=LagSqrt, y=Estimate))+
      labs(x="Sqrt (Time Lag)", y=paste0("Jaccard dissimilarity (coverage = ", sc2, ")"), color="Order.q")+
      geom_point(data=rawdata_0, aes(x=LagSqrt, y=Estimate), size=1.5, shape=16, alpha=0.3, color="#889eae")+
      scale_color_manual(labels=c("Order q0"), values="#889eae")+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free", nrow = 1)
    p_0 <- p_0 + geom_vline(data = rawdata_0[rawdata_0$Taxa %in% c("Fish", "Fish (Trait)"), ], 
                            aes(xintercept = 4.25), linetype = "dashed", color = "black")
    
    p_2 <- ggplot(data=rawdata_2, aes(x=LagSqrt, y=Estimate))+
      labs(x="Sqrt (Time Lag)", y=paste0("Jaccard dissimilarity (coverage = ", sc2, ")"), color="Order.q")+
      geom_point(data=rawdata_2, aes(x=LagSqrt, y=Estimate), size=1.5, shape=16, alpha=0.3, color="#1b3765")+
      scale_color_manual(labels=c("Order q2"), values="#1b3765")+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free", nrow = 1)
    p_2 <- p_2 + geom_vline(data = rawdata_2[rawdata_2$Taxa %in% c("Fish", "Fish (Trait)"), ], 
                            aes(xintercept = 4.25), linetype = "dashed", color = "black")
    (p_beta_0 <- p_0 +
        geom_line(data =  model_0[model_0$lpval < 0.1,], aes(x=LagSqrt, y = fit.x, linetype=as.factor(lpvalS)), size=0.8, color="#889eae")+
        geom_ribbon(data =  model_0[model_0$lpval < 0.1,], aes(x=LagSqrt, y=fit.x, ymin = lwr.x, ymax = upr.x), color="#889eae", fill="white", alpha = 0.1)+
        scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
        facet_wrap(~ Taxa, scales="free_x", nrow=1))
    (p_beta_2 <- p_2 +
        geom_line(data =  model_2[model_2$lpval < 0.1,], aes(x=LagSqrt, y = fit.x, linetype=as.factor(lpvalS)), size=0.8, color="#1b3765")+
        geom_ribbon(data =  model_2[model_2$lpval < 0.1,], aes(x=LagSqrt, y=fit.x, ymin = lwr.x, ymax = upr.x), color="#1b3765", fill="white", alpha = 0.1)+
        scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
        facet_wrap(~ Taxa, scales="free_x", nrow=1))
    
  }else{
    p_0 <- ggplot(data=rawdata_0, aes(x=LagSqrt, y=Estimate))+
      labs(x="Sqrt (Time Lag)", y=paste0("Jaccard dissimilarity (coverage = ", sc2, ")"), color="Order.q")+
      geom_point(data=rawdata_0, aes(x=LagSqrt, y=Estimate), size=1.5, shape=16, alpha=0.3, color="#889eae")+
      scale_color_manual(labels=c("Order q0"), values="#889eae")+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free", nrow = 1)
    p_2 <- ggplot(data=rawdata_2, aes(x=LagSqrt, y=Estimate))+
      labs(x="Sqrt (Time Lag)", y=paste0("Jaccard dissimilarity (coverage = ", sc2, ")"), color="Order.q")+
      geom_point(data=rawdata_2, aes(x=LagSqrt, y=Estimate), size=1.5, shape=16, alpha=0.3, color="#1b3765")+
      scale_color_manual(labels=c("Order q2"), values="#1b3765")+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Order.q)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free", nrow = 1)
    
    (p_beta_0 <- p_0 +
        geom_line(data =  model_0[model_0$ypval < 0.1,], aes(x=LagSqrt, y = fit.y, linetype=as.factor(ypvalS)), size=0.8, color="#889eae")+
        geom_ribbon(data =  model_0[model_0$ypval < 0.1,], aes(x=LagSqrt, y=fit.y, ymin = lwr.y, ymax = upr.y), color="#889eae", fill="white", alpha = 0.1)+
        scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
        facet_wrap(~ Taxa, scales="free_x", nrow=1))
    (p_beta_2 <- p_2 +
        geom_line(data =  model_2[model_2$ypval < 0.1,], aes(x=LagSqrt, y = fit.y, linetype=as.factor(ypvalS)), size=0.8, color="#1b3765")+
        geom_ribbon(data =  model_2[model_2$ypval < 0.1,], aes(x=LagSqrt, y=fit.y, ymin = lwr.y, ymax = upr.y), color="#1b3765", fill="white", alpha = 0.1)+
        scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
        facet_wrap(~ Taxa, scales="free_x", nrow=1))      
  }
  
  p_mix <- p_beta_0/p_beta_2 + plot_layout(guides = "collect", axes = "collect")
  
  return(p_mix)
  
}

(main_linear_Beta_C5 <- plot_β(mrawBetaLag22, mpolypredBetaLag22, group="Northern Range", sc2="Cmin", polynomial="linear"))
(sm_linear_Beta_C5 <- plot_β(mrawBetaLag22, mpolypredBetaLag22, group="Northern Range", sc2="Cmax", polynomial="linear"))

(main_nonlinear_Beta_C5 <- plot_β(mrawBetaLag, mpolypredBetaLag, group="Northern Range", sc2="Cmin", polynomial="nonlinear"))
(sm_nonlinear_Beta_C5 <- plot_β(mrawBetaLag, mpolypredBetaLag, group="Northern Range", sc2="Cmax", polynomial="nonlinear"))

#ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_LinearChange_Beta_Cmin.png"), main_linear_Beta_C5, width = 10, height = 6.5)        # legend aes fix on PowerPoint
#ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_LinearChange_Beta_Cmax.png"), sm_linear_Beta_C5, width = 10, height = 6.5)          # legend aes fix on PowerPoint
#ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_NonlinearChange_Beta_Cmin.png"), main_nonlinear_Beta_C5, width = 10, height = 6.5)  # legend aes fix on PowerPoint
#ggsave(paste0(PlotsPath, "/Chapter5/S6_C5_NonlinearChange_Beta_Cmax.png"), sm_nonlinear_Beta_C5, width = 10, height = 6.5)    # legend aes fix on PowerPoint


# Chapter 6: ===================================================================

plot_nonlinearC6 <- function(raw, pred, orderq=NULL, cov=NULL, abline=NULL, type=NULL){

  raw_int <- subset(raw, raw$SC2==cov & raw$Order.q==orderq)
  raw_int <- raw_int[!raw_int$Group=="Northern Range",]
  pred_int <- subset(pred, pred$SC2==cov & pred$Order.q==orderq)
  pred_int <- pred_int[!pred_int$Group=="Northern Range",]
  
  if (type=="alpha"){
    names(raw_int)[names(raw_int)=="Session"] <- "Time"
    names(pred_int)[names(pred_int)=="Session"] <- "Time"
  }else{
    names(raw_int)[names(raw_int)=="LagSqrt"] <- "Time"
    names(pred_int)[names(pred_int)=="LagSqrt"] <- "Time"
  }
  
  (p <- ggplot(data=raw_int, aes(x=Time, y=Estimate, color=as.factor(Group)))+
      labs(x="Temporal session", y=paste0("q", orderq, " (", cov, ")"), color="Disturbance")+
      geom_point(data=raw_int, aes(x=Time, y=Estimate, color=as.factor(Group)), size=1.5, shape=16, alpha=0.4)+
      scale_color_manual(labels=c("Disturbed", "Undisturbed"), values=c("#CE780F", "#0A7D9D"), breaks=c("Disturbed","Undisturbed"))+
      #geom_errorbar(aes(ymin=LCL, ymax=UCL, group=as.factor(Group)), width=.5,
      #              position=position_dodge(0.05), alpha=0.7, color="grey")+
      theme_classic()+
      facet_wrap(~ Taxa, scales = "free"))
  (p <- p +
      geom_line(data =  pred_int, aes(x=Time, y = fit.y, color=as.factor(Group), size=as.factor(ypvalS), linetype=as.factor(ypvalS)), size=0.8)+
      geom_ribbon(data =  pred_int, aes(x=Time, y=fit.y, ymin = lwr.y, ymax = upr.y, color=as.factor(Group)), fill="white", alpha = 0.1)+
      scale_size_manual(labels=c("0", "1"), values = c(1, 1.2), guide = "none", breaks = c(0,1)) +
      scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none",breaks = c(0,1)) +
      scale_color_manual(labels=c("Disturbed", "Undisturbed"), values=c("#CE780F", "#0A7D9D"), breaks=c("Disturbed","Undisturbed"))+
      facet_wrap(~ Taxa, scales="free", nrow = 1))
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


AlphaCminq0 <- plot_nonlinearC6(mrawAlpha, mpolypredAlpha[mpolypredAlpha$ypval<0.1,], orderq = 0, cov="Cmin", abline=19.5, type="alpha")
AlphaCmaxq0 <- plot_nonlinearC6(mrawAlpha, mpolypredAlpha[mpolypredAlpha$ypval<0.1,], orderq = 0, cov="Cmax", abline=19.5, type="alpha")

AlphaCminq2 <- plot_nonlinearC6(mrawAlpha, mpolypredAlpha[mpolypredAlpha$ypval<0.1,], orderq = 2, cov="Cmin", abline=19.5, type="alpha")
AlphaCmaxq2 <- plot_nonlinearC6(mrawAlpha, mpolypredAlpha[mpolypredAlpha$ypval<0.1,], orderq = 2, cov="Cmax", abline=19.5, type="alpha")

BetaCminq0 <- plot_nonlinearC6(mrawBetaLag, mpolypredBetaLag[mpolypredBetaLag$ypval<0.1,], orderq = 0, cov="Cmin", abline=4.5, type="beta")
BetaCmaxq0 <- plot_nonlinearC6(mrawBetaLag, mpolypredBetaLag[mpolypredBetaLag$ypval<0.1,], orderq = 0, cov="Cmax", abline=4.5, type="beta")

BetaCminq2 <- plot_nonlinearC6(mrawBetaLag, mpolypredBetaLag[mpolypredBetaLag$ypval<0.1,], orderq = 2, cov="Cmin", abline=4.5, type="beta")
BetaCmaxq2 <- plot_nonlinearC6(mrawBetaLag, mpolypredBetaLag[mpolypredBetaLag$ypval<0.1,], orderq = 2, cov="Cmax", abline=4.5, type="beta")


AlphaCmin <- AlphaCminq0/AlphaCminq2 + plot_layout(guides = "collect", axes = "collect")
AlphaCmax <- AlphaCmaxq0/AlphaCmaxq2 + plot_layout(guides = "collect", axes = "collect")

BetaCmin <- BetaCminq0/BetaCminq2 + plot_layout(guides = "collect", axes = "collect")
BetaCmax <- BetaCmaxq0/BetaCmaxq2 + plot_layout(guides = "collect", axes = "collect")


#ggsave(paste0(PlotsPath, "/Chapter6/S6_C6_NonlinearChange_Alpha_Cmin.png"), AlphaCmin, width = 10, height = 6.5)
#ggsave(paste0(PlotsPath, "/Chapter6/S6_C6_NonlinearChange_Alpha_Cmax.png"), AlphaCmax, width = 10, height = 6.5)

#ggsave(paste0(PlotsPath, "/Chapter6/S6_C6_NonlinearChange_Beta_Cmin.png"), BetaCmin, width = 10, height = 6.5)
#ggsave(paste0(PlotsPath, "/Chapter6/S6_C6_NonlinearChange_Beta_Cmax.png"), BetaCmax, width = 10, height = 6.5) # legend aesthetics fixed in PowerPointw  when needed


# Chapter 6 (DU) ===============================================================
mrawBetaDU22 <- subset(mrawBetaDU22, mrawBetaDU22$SC2=="Cmin")                 # or Cmax
mpolypredBetaDU22 <- subset(mpolypredBetaDU22, mpolypredBetaDU22$SC2=="Cmin")  # or Cmax

(pDU <- ggplot(data=mrawBetaDU22, aes(x=Session, y=Estimate, color=Order.q))+
  labs(x="Time Session", y="Jaccard dissimilarity Cmin", color="Order.q")+
  geom_point(data=mrawBetaDU22, aes(x=Session, y=Estimate, color=Order.q), size=1.5, shape=16, alpha=0.6)+
  theme_classic()+
  scale_y_continuous(limits = c(-0.15, NA)) +
  scale_color_manual(labels=c("Order q0", "Order q2"), values=c("#889eae", "#1b3765"))+
  facet_wrap(~ Taxa, scales = "free", nrow = 1))
pDU <- pDU + geom_vline(data = mrawBetaDU22[mrawBetaDU22$Taxa %in% c("Fish", "Fish (Trait)"), ], 
                        aes(xintercept = 19.5), linetype = "dashed", color = "black")
pDU
pDU <- pDU +
       geom_line(data =  mpolypredBetaDU22[mpolypredBetaDU22$lpval<0.1,], aes(x=Session, y = fit.x, color=Order.q, linetype=as.factor(lpvalS)),size=0.8)+
       geom_ribbon(data =  mpolypredBetaDU22[mpolypredBetaDU22$lpval<0.1,], aes(x=Session, y=fit.x, ymin = lwr.x, ymax = upr.x, color=Order.q), fill="white", alpha = 0.1)+
       #scale_size_manual(labels=c("0", "1"), values = c(1, 1.2), guide = "none") +
       scale_linetype_manual(labels=c("0", "1"), values = c("dashed", "solid"), guide = "none", breaks = c(0,1)) +
       scale_color_manual(labels=c("Order q0", "Order q2"), values=c("#889eae", "#1b3765"), breaks = c(0,2))+
       facet_wrap(~ Taxa, scales="free_x", nrow=1)
pDU

#ggsave(paste0(PlotsPath, "/Chapter6/S6_C6_LinearChange_DU_Cmin.png"), pDU, width = 10, height = 3.25)  # or Cmax


# End of script ################################################################
################################################################################