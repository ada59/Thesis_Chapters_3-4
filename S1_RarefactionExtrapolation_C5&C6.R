################################################################################
# S1 Rarefaction/ Extrapolation (sample size and coverage)
# AFE
# Dec 2023
################################################################################


# libraries: ===================================================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(vegan)
library(data.table)
library(devtools)
library(iNEXT.3D)
library(patchwork)


# read data: ===================================================================
rm(list=ls())
Thesis_GeneralMethods <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods"
mainanalysisRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/mainanalysisRData"
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"

all_files <- list.files(mainanalysisRData, pattern = "\\.RData$", full.names = TRUE)

# Loop through the list and load each RData file
for (file in all_files) {
  load(file)
}

load("C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods/FormattedTrait_Trinidad/dist_Euc_av.RData")
acr <- read.csv(paste0(Thesis_GeneralMethods,"/acronymsNR.csv"), h=T)
acr$acronym2 <- paste0(substring(acr$acronym,1,1), "_", sub(".", "", acr$acronym))
identical(sort(unique(colnames(dist_Euc_av))), sort(unique(acr$acronym2)))      # TRUE
rownames(dist_Euc_av) <- acr$name[match(rownames(dist_Euc_av), acr$acronym2)]
colnames(dist_Euc_av) <- acr$name[match(colnames(dist_Euc_av), acr$acronym2)]
dist_Euc_av <- (dist_Euc_av - min(dist_Euc_av)) / (max(dist_Euc_av) - min(dist_Euc_av))            # Re-scale between 0 & 1
save(dist_Euc_av, file=paste0(mainanalysisRData, "/dist_Euc_av.RData"))         # species names formatted properly using acr


################################################################################
# curves (one x metacom) =======================================================

# PA ---------------------------------------------------------------------------
#out.rawPA <- lapply(lAllPA, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})

#out.rawPA_d <- lapply(lPA_d, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})

#out.rawPA_u <- lapply(lPA_u, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})
                                     


# Agg --------------------------------------------------------------------------
#out.rawAgg <- lapply(lAllAbuAgg, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#})

#out.rawAgg_d <- lapply(lAbuAgg_d, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#})

#out.rawAgg_u <- lapply(lAbuAgg_u, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)                    
#})


# Community abundance ==========================================================
#out.rawAbu <- lapply(lAllAbu, function(x) {lapply(x, function(y) {iNEXT3D(data = y, diversity = "TD",
#                                                                         q = c(0, 1, 2), datatype = "abundance", nboot = 100)})})


#save(out.rawPA, out.rawPA_d, out.rawPA_u, file = paste0(Results_AlphaDivRData, "/S1_out.rawPA_Alpha_Taxonomic.RData"))
#save(out.rawAgg, out.rawAgg_d, out.rawAgg_u, file = paste0(Results_AlphaDivRData, "/S1_out.rawAgg_Alpha_Taxonomic.RData"))
#save(out.rawAbu, file = paste0(Results_AlphaDivRData, "/S1_out.rawAbu_Alpha_Taxonomic.RData"))



################################################################################
# Trait Diversity:--------------------------------------------------------------

#out.rawPA_FD <- iNEXT3D(data = lAllPA$fish, diversity = "FD",
#                        FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                        q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)       # default tau is mean distance
#out.rawPA_dFD <- iNEXT3D(data = lPA_d$fish, diversity = "FD",
#                        FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                        q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#out.rawPA_uFD <- iNEXT3D(data = lPA_u$fish, diversity = "FD",
#                         FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                         q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#out.rawAgg_FD <- iNEXT3D(data = lAllAbuAgg$fish, diversity = "FD",
#                        FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                       q = c(0, 1, 2), datatype = "abundance", nboot = 100)           # default tau is mean distance
#out.rawAgg_dFD <- iNEXT3D(data = lAbuAgg_d$fish, diversity = "FD", 
#                         FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                         q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#out.rawAgg_uFD <- iNEXT3D(data = lAbuAgg_u$fish, diversity = "FD", 
#                         FDdistM = dist_Euc_av, FDtype = "tau_values", FDtau = NULL,
#                         q = c(0, 1, 2), datatype = "abundance", nboot = 100)


#save(out.rawPA_FD, out.rawPA_dFD, out.rawPA_uFD, file = paste0(Results_AlphaDivRData, "/S1_out.rawPA_Alpha_FishTrait.RData"))
#save(out.rawAgg_FD, out.rawAgg_dFD, out.rawAgg_uFD, file = paste0(Results_AlphaDivRData, "/S1_out.rawAgg_Alpha_FishTrait.RData"))


################################################################################
# data info ====================================================================

create_counts_vector <- function(df) {
  counts_vector <- c(
    Singletons = sum(rowSums(df)==1),
    Doubletons = sum(rowSums(df)==2)
  )
  return(counts_vector)
}
functionDataInfo <- function(dt, type){
  
  if (type=="PA"){
    DataInfo <-   lapply(dt, function(x) {
      DataInfo3D(data = x, diversity = "TD",
                 datatype = "incidence_raw")
    })  
    SingDoub <- lapply(dt,  function(z){lapply(z, create_counts_vector)})
  }else{
    DataInfo <- lapply(dt, function(x) {
      DataInfo3D(data = x, diversity = "TD",
                 datatype = "abundance")
    }) 
    SingDoub <- lapply(dt,  function(z){lapply(z, function(z) {c(sum(z==1), sum(z==2))})})
  }
  
  DataInfo <- as.data.frame(do.call(rbind, DataInfo))
  DataInfo$Taxa <- str_split_fixed(rownames(DataInfo), "\\.", 2)[,1]
  DataInfo$Session <- str_split_fixed(rownames(DataInfo), "\\.", 2)[,2]
  rownames(DataInfo) <- 1:nrow(DataInfo)
  
  SingDoub <- lapply(SingDoub,  function(z){as.data.frame(do.call(rbind, z))})
  SingDoub <- as.data.frame(do.call(rbind, SingDoub))
  SingDoub$Taxa <- str_split_fixed(rownames(SingDoub), "\\.", 2)[,1]
  SingDoub$Session <- str_split_fixed(rownames(SingDoub), "\\.", 2)[,2]
  rownames(SingDoub) <- 1:nrow(SingDoub)
  DataInfo_Table <- merge(DataInfo, SingDoub, by=c("Taxa", "Session"))
  DataInfo_Table$Taxa <- recode_factor(DataInfo_Table$Taxa, "fish" = "Fish",
                               "bi"="Invertebrates",
                               "dia"="Diatoms")
  DataInfo_Table$Session <- as.numeric(DataInfo_Table$Session)
  DataInfo_Table <- DataInfo_Table[order(DataInfo_Table$Taxa, DataInfo_Table$Session), ]
  DataInfo_Table$Taxa <- factor(DataInfo_Table$Taxa,
                                levels = c("Diatoms", "Invertebrates", "Fish"))
  return(DataInfo_Table)
}

# Presence-absence: ------------------------------------------------------------  
DataInfoPA_Table <- functionDataInfo(lAllPA, type="PA")
DataInfoPA_Table <- within(DataInfoPA_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))       # NR
DataInfoPA_d_Table <- functionDataInfo(lPA_d, type="PA")
DataInfoPA_d_Table <- within(DataInfoPA_d_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))   # DIS
DataInfoPA_u_Table <- functionDataInfo(lPA_u, type="PA")
DataInfoPA_u_Table <- within(DataInfoPA_u_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))   # UNDIS


# Aggregated abundances: -------------------------------------------------------
DataInfoAgg_Table <- functionDataInfo(lAllAbuAgg, type="Abu")
DataInfoAgg_Table <- within(DataInfoAgg_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V2"] <- "Doubletons"               # NR


DataInfoAgg_d_Table <- functionDataInfo(lAbuAgg_d, type="Abu")
DataInfoAgg_d_Table <- within(DataInfoAgg_d_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_d_Table)[names(DataInfoAgg_d_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_d_Table)[names(DataInfoAgg_d_Table)=="V2"] <- "Doubletons"           # DIS


DataInfoAgg_u_Table <- functionDataInfo(lAbuAgg_u, type="Abu")
DataInfoAgg_u_Table <- within(DataInfoAgg_u_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_u_Table)[names(DataInfoAgg_u_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_u_Table)[names(DataInfoAgg_u_Table)=="V2"] <- "Doubletons"           # UNDIS


################################################################################
# save: ========================================================================
table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"
write.csv(DataInfoPA_Table, file=paste0(table_path, "/DataInfoPA_Table.csv"), row.names = F)
write.csv(DataInfoPA_d_Table, file=paste0(table_path, "/DataInfoPA_d_Table.csv"), row.names = F)
write.csv(DataInfoPA_u_Table, file=paste0(table_path, "/DataInfoPA_u_Table.csv"), row.names = F)


write.csv(DataInfoAgg_Table, file=paste0(table_path, "/DataInfoAgg_Table.csv"), row.names = F)
write.csv(DataInfoAgg_d_Table, file=paste0(table_path, "/DataInfoAgg_d_Table.csv"), row.names = F)
write.csv(DataInfoAgg_u_Table, file=paste0(table_path, "/DataInfoAgg_u_Table.csv"), row.names = F)


################################################################################
# plots: =======================================================================
all_files <- list.files(Results_AlphaDivRData, pattern = "\\.RData$", full.names = TRUE)

# Loop through the list and load each RData file
for (file in all_files) {
  load(file)
}


# panel of 3 x 4 ===============================================================
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
range(out.rawPA$fish$TDiNextEst$size_based$nT)
extract_size_based <- function(dt, dimension="tax") {
  if(dimension=="tax"){
    result <- lapply(dt, function(x) {x$TDiNextEst$size_based})
    result <- as.data.frame(do.call(rbind, result))
    result$Taxa <- str_split_fixed(rownames(result), "\\.", 2)[,1]
    rownames(result) <- 1:nrow(result)
    result$Order.q <- as.factor(result$Order.q)
    result$Taxa <- recode_factor(result$Taxa, "fish" = "Fish",
                                 "bi"="Invertebrates",
                                 "dia"="Diatoms")
    result$Taxa <- factor(result$Taxa,
                          levels = c("Diatoms", "Invertebrates", "Fish"))
  }else{
    result <- dt[["FDiNextEst"]][["size_based"]]
    result <- result[,! names(result) =="Tau"]
    result$Taxa <- "Fish (Trait)"
    result$Order.q <- as.factor(result$Order.q)
    names(result)[names(result)=="qFD"] <- "qTD"
    names(result)[names(result)=="qFD.LCL"] <- "qTD.LCL"
    names(result)[names(result)=="qFD.UCL"] <- "qTD.UCL"
    
  }
  
  return(result)
}
out.rawPA_plot <- extract_size_based(out.rawPA, dimension = "tax")
out.rawPA_d_plot <- extract_size_based(out.rawPA_d, dimension = "tax")
out.rawPA_u_plot <- extract_size_based(out.rawPA_u, dimension = "tax")

identical(names(out.rawPA_plot), names(out.rawPA_d_plot)) # TRUE

out.rawAgg_plot <- extract_size_based(out.rawAgg, dimension = "tax")
out.rawAgg_d_plot <- extract_size_based(out.rawAgg_d, dimension = "tax")
out.rawAgg_u_plot <- extract_size_based(out.rawAgg_u, dimension = "tax")

out.rawPA_FD_plot <- extract_size_based(out.rawPA_FD, dimension = "trait")
out.rawPA_FD_d_plot <- extract_size_based(out.rawPA_dFD, dimension = "trait")
out.rawPA_FD_u_plot <- extract_size_based(out.rawPA_uFD, dimension = "trait")

out.rawAgg_FD_plot <- extract_size_based(out.rawAgg_FD, dimension = "trait")
out.rawAgg_FD_d_plot <- extract_size_based(out.rawAgg_dFD, dimension = "trait")
out.rawAgg_FD_u_plot <- extract_size_based(out.rawAgg_uFD, dimension = "trait")



out.rawPA_S_C <- rbind(out.rawPA_plot, out.rawPA_FD_plot)
out.rawPA_S_C$Group <- "Northern Range"
out.rawPA_S_C$Type <- "PA"
out.rawPAd_S_C <- rbind(out.rawPA_d_plot, out.rawPA_FD_d_plot)
out.rawPAd_S_C$Group <- "Disturbed"
out.rawPAd_S_C$Type <- "PA"
out.rawPAu_S_C <- rbind(out.rawPA_u_plot, out.rawPA_FD_u_plot)
out.rawPAu_S_C$Group <- "Undisturbed"
out.rawPAu_S_C$Type <- "PA"

out.rawAgg_S_C <- rbind(out.rawAgg_plot, out.rawAgg_FD_plot)
out.rawAgg_S_C$Group <- "Northern Range"
out.rawAgg_S_C$Type <- "Agg"
out.rawAggd_S_C <- rbind(out.rawAgg_d_plot, out.rawAgg_FD_d_plot)
out.rawAggd_S_C$Group <- "Disturbed"
out.rawAggd_S_C$Type <- "Agg"
out.rawAggu_S_C <- rbind(out.rawAgg_u_plot, out.rawAgg_FD_u_plot)
out.rawAggu_S_C$Group <- "Undisturbed"
out.rawAggu_S_C$Type <- "Agg"


for_observed <- bind_rows(out.rawPA_S_C, out.rawPAd_S_C, out.rawPAu_S_C,
                     out.rawAgg_S_C, out.rawAggd_S_C, out.rawAggu_S_C)
for_observed$SC <- round(for_observed$SC, digits=4)
for_observed$Concat <- paste0(for_observed$Assemblage, "_", for_observed$Taxa, "_",
                              for_observed$Group, "_", for_observed$Type, "_", for_observed$SC)
for_observed <- subset(for_observed, for_observed$Method=="Observed")
save(for_observed, file=paste0(mainanalysisRData, "/for_observed.RData"))


# Size vs coverage =============================================================

size_vs_coverage <- function(dt, type=NULL){
  plotPA <- dt
  
  dt_plot <- subset(plotPA, plotPA$Order.q=="0")
  dt_plot_ref <- subset(dt_plot, dt_plot$Method=="Observed")
  dt_plot_r <- subset(dt_plot, dt_plot$Method=="Rarefaction")
  dt_plot_e <- subset(dt_plot, dt_plot$Method=="Extrapolation")
  
  if ("nT" %in% names(plotPA)) { # presence-absence
  p <- ggplot(data=dt_plot, aes(x=nT, y=SC, group=Assemblage, color=Taxa)) +
    labs(title=type)+
    geom_point(data=dt_plot_ref, aes(x=nT, y=SC, group=Assemblage))+
    geom_line(data=dt_plot_r, aes(x=nT, y=SC, group=Assemblage), linetype="solid", alpha=0.4)+
    geom_line(data=dt_plot_e, aes(x=nT, y=SC, group=Assemblage), linetype="dotdash")+
    scale_color_manual(values=c("#aacc7f","#cca105", "#149299"))+
    theme_classic()+
    scale_y_continuous(
      breaks = c(0.5,0.6,0.7,0.8,0.9,1),
      limits = c(0.4, 1))+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))+
    geom_hline(yintercept=0.8, linetype=2, color = "black")+
    facet_wrap(~Taxa, nrow=1)}else{ # agg abus
      p <- ggplot(data=dt_plot, aes(x=m, y=SC, group=Assemblage, color=Taxa)) +
        labs(title=type)+
        geom_point(data=dt_plot_ref, aes(x=m, y=SC, group=Assemblage))+
        geom_line(data=dt_plot_r, aes(x=m, y=SC, group=Assemblage), linetype="solid", alpha=0.4)+
        geom_line(data=dt_plot_e, aes(x=m, y=SC, group=Assemblage), linetype="dotdash")+
        scale_color_manual(values=c("#aacc7f","#cca105", "#149299"))+
        theme_classic()+
        theme(legend.position = "none",
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(size = 12),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10))+
        geom_hline(yintercept=0.8, linetype=2, color = "black")+
        facet_wrap(~Taxa, nrow=1, scales = "free_x")  
    }
  return(p)
}
  

pAll <- size_vs_coverage(out.rawPA_plot, type="Northern Range")
pAll
pD <- size_vs_coverage(out.rawPA_d_plot, type="Disturbed")
pD
pU <- size_vs_coverage(out.rawPA_u_plot, type="Undisturbed")
pU

resultPA <- pAll/pD/pU
gt <- patchwork::patchworkGrob(resultPA)
gt2 <- gridExtra::grid.arrange(gt, left = "Coverage", bottom = "Size (sampling units)")
ggsave(paste0(plot_path, "/Size_vs_CoveragePA.png"), gt2, width = 8, height = 10)
  

pAggAll <- size_vs_coverage(out.rawAgg_plot, type="Northern Range")
pAggAll
pAggD <- size_vs_coverage(out.rawAgg_d_plot, type="Disturbed")
pAggD
pAggU <- size_vs_coverage(out.rawAgg_u_plot, type="Undisturbed")
pAggU

resultAgg <- pAggAll/pAggD/pAggU
gt <- patchwork::patchworkGrob(resultAgg)
gt2 <- gridExtra::grid.arrange(gt, left = "Coverage", bottom = "Size (number of individuals)")
ggsave(paste0(plot_path, "/Size_vs_CoverageAggregated.png"), gt2, width = 8, height = 10)


# Diversity vs coverage ========================================================
diversity_vs_coverage <- function(dt, order=NULL) {
  if ("m" %in% names(dt)){
    names(dt)[names(dt)=="m"] <- "nT"
  }else{names(dt) <- names(dt)}
  plotPA <- dt
  dt_plot <- subset(plotPA, plotPA$Order.q==order)
  dt_plot_ref <- subset(dt_plot, dt_plot$Method=="Observed")
  dt_plot_r <- subset(dt_plot, dt_plot$Method=="Rarefaction")
  dt_plot_e <- subset(dt_plot, dt_plot$Method=="Extrapolation")
  
  p <- ggplot(data=dt_plot, aes(x=nT, y=qTD, group=Assemblage, color=Taxa)) +
    labs(y=paste0("Order q", order))+
    geom_point(data=dt_plot_ref, aes(x=nT, y=qTD, group=Assemblage))+
    geom_line(data=dt_plot_r, aes(x=nT, y=qTD, group=Assemblage), linetype="solid", alpha=0.4)+
    geom_line(data=dt_plot_e, aes(x=nT, y=qTD, group=Assemblage), linetype="dotdash")+
    scale_color_manual(values=c("#aacc7f","#cca105", "#149299", "#9900CC"))+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))+
    #geom_vline(xintercept=0.8, linetype=2, color = "black")+
    facet_wrap(~Taxa, nrow=1, scales = "free_x")
  #if ("nT" %in% names(plotPA)) {
  #  p <- p + scale_x_continuous(breaks = c(0.5,0.8,1),
  #    limits = c(0.4, 1)) +
  #    scale_y_continuous(
  #    limits = c(5, 30))
  #}else{
  #  p <- p
  #} (use this is plotting coverage-based instead)
    
  return(p)
  
}

layout_fun <- function(x,y,z){
  result <- x/y/z
  gt <- patchwork::patchworkGrob(result)
  gt2 <- gridExtra::grid.arrange(gt, bottom = "Size")
  return(gt2)
}

# Order q0 (richness): ---------------------------------------------------------
pTD0 <- diversity_vs_coverage(out.rawPA_S_C, order="0")
pTD0
pTD0_d <- diversity_vs_coverage(out.rawPAd_S_C, order="0")
pTD0_d
pTD0_u <- diversity_vs_coverage(out.rawPAu_S_C, order="0")
pTD0_u

pTDAgg0 <- diversity_vs_coverage(out.rawAgg_S_C, order="0")
pTDAgg0
pTDAgg0_d <- diversity_vs_coverage(out.rawAggd_S_C, order="0")
pTDAgg0_d
pTDAgg0_u <- diversity_vs_coverage(out.rawAggu_S_C, order="0")
pTDAgg0_u


# Order q1 (Shannon): ----------------------------------------------------------
pTD1 <- diversity_vs_coverage(out.rawPA_S_C, order="1")
pTD1
pTD1_d <- diversity_vs_coverage(out.rawPAd_S_C, order="1")
pTD1_d
pTD1_u <- diversity_vs_coverage(out.rawPAu_S_C, order="1")
pTD1_u

pTDAgg1 <- diversity_vs_coverage(out.rawAgg_S_C, order="1")
pTDAgg1
pTDAgg1_d <- diversity_vs_coverage(out.rawAggd_S_C, order="1")
pTDAgg1_d
pTDAgg1_u <- diversity_vs_coverage(out.rawAggu_S_C, order="1")
pTDAgg1_u


# Order q=2 (Simpson): ---------------------------------------------------------
pTD2 <- diversity_vs_coverage(out.rawPA_S_C, order="2")
pTD2
pTD2_d <- diversity_vs_coverage(out.rawPAd_S_C, order="2")
pTD2_d
pTD2_u <- diversity_vs_coverage(out.rawPAu_S_C, order="2")
pTD2_u

pTDAgg2 <- diversity_vs_coverage(out.rawAgg_S_C, order="2")
pTDAgg2
pTDAgg2_d <- diversity_vs_coverage(out.rawAggd_S_C, order="2")
pTDAgg2_d
pTDAgg2_u <- diversity_vs_coverage(out.rawAggu_S_C, order="2")
pTDAgg2_u



resultTD <- layout_fun(pTD0, pTD1, pTD2)
ggsave(paste0(plot_path, "/Size_vs_Diversity_NorthernRangePA.png"), resultTD, width = 8, height = 10)

resultTD_d <- layout_fun(pTD0_d, pTD1_d, pTD2_d)
ggsave(paste0(plot_path, "/Size_vs_Diversity_DistubedPA.png"), resultTD_d, width = 8, height = 10)

resultTD_u <- layout_fun(pTD0_u,pTD1_u,pTD2_u)
ggsave(paste0(plot_path, "/Size_vs_Diversity_UndistubedPA.png"), resultTD_u, width = 8, height = 10)


resultAggTD <- layout_fun(pTDAgg0, pTDAgg1, pTDAgg2)
ggsave(paste0(plot_path, "/Size_vs_Diversity_NorthernRangeAgg.png"), resultAggTD, width = 10, height = 10)

resultAggTD_d <- layout_fun(pTDAgg0_d, pTDAgg1_d,pTDAgg2_d)
ggsave(paste0(plot_path, "/Size_vs_Diversity_DistubedAgg.png"), resultAggTD_d, width = 10, height = 10)

resultAggTD_u <- layout_fun(pTDAgg0_u, pTDAgg1_u,pTDAgg2_u)
ggsave(paste0(plot_path, "/Size_vs_Diversity_UndistubedAgg.png"), resultAggTD_u, width = 10, height = 10)


################################################################################
# curves (pooling x 5 y) =======================================================
# To Be Completed for Ms (following recommendation by Anne Chao)

# End of script#################################################################