################################################################################
# S1 Sample size & coverage
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
load ("lAllPA.RData")
load ("lAllAbuAgg.RData") 
load ("lAllAbu.RData") 

names(lAllPA) <- c("Fish&Crustaceans", "BenthicInvertebrates", "Diatoms")
names(lAllAbuAgg) <- c("Fish&Crustaceans", "BenthicInvertebrates", "Diatoms")
names(lAllAbu) <- c("Fish&Crustaceans", "BenthicInvertebrates", "Diatoms")


# create lists sensitivity (fish with and without crustaceans)==================
sen_fishPA <- lAllPA[[1]]
sen_fishPA <- lapply(sen_fishPA, function(x) {x[!rownames(x) %in% c("Macrobrachium_spp", "Eudaniela_garmani",
                                                                    "Atya_spp"),]})

elements_to_remove <- c("Macrobrachium_spp", "Eudaniela_garmani",
                        "Atya_spp")
remove_elements <- function(my_list, elements_to_remove) {
  map(my_list, ~ discard(.x, names(.x) %in% elements_to_remove))
}
sen_fishAgg <- lAllAbuAgg[[1]]
sen_fishAgg <- remove_elements(sen_fishAgg, elements_to_remove) # OK

sen_fishAbu <- lAllAbu[[1]]
sen_fishAbu <- lapply(sen_fishAbu, function(x) {x[!rownames(x) %in% c("Macrobrachium_spp", "Eudaniela_garmani",
                                                                      "Atya_spp"),]}) 
                                                                    
lAllPA[["Fish"]] <- sen_fishPA
lAllAbuAgg[["Fish"]] <- sen_fishAgg
lAllAbu[["Fish"]] <- sen_fishAbu


# create lists D & U ===========================================================
dtPA_d <- lapply(lAllPA, function(x) {lapply(x, function(y) {y[names(y) %like% "_d"]})})
dtPA_u <- lapply(lAllPA, function(x) {lapply(x, function(y) {y[names(y) %like% "_u"]})})

dtAgg_d <- lapply(lAllAbu, function(x) {lapply(x, function(y) {y[names(y) %like% "_d"]})})
dtAgg_d <- lapply(dtAgg_d, function(x) {lapply(x, function(y) {rowSums(y)})})

dtAgg_u <- lapply(lAllAbu, function(x) {lapply(x, function(y) {y[names(y) %like% "_u"]})})
dtAgg_u <- lapply(dtAgg_u, function(x) {lapply(x, function(y) {rowSums(y)})})



################################################################################
# curves (one x metacom) =======================================================

# PA ---------------------------------------------------------------------------
#out.rawPA <- lapply(lAllPA, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})
#save(out.rawPA, file="out.rawPA.RData")
#out.rawPA_d <- lapply(dtPA_d, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})
#save(out.rawPA_d, file="out.rawPA_d.RData")
#out.rawPA_u <- lapply(dtPA_u, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "incidence_raw", nboot = 100)
#})
#save(out.rawPA_u, file="out.rawPA_u.RData")

# Agg --------------------------------------------------------------------------
#out.rawAgg <- lapply(lAllAbuAgg, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#})
#save(out.rawAgg, file="out.rawAgg.RData")
#out.rawAgg_d <- lapply(dtAgg_d, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#})
#save(out.rawAgg_d, file="out.rawAgg_d.RData")
#out.rawAgg_u <- lapply(dtAgg_u, function(x) {
#  iNEXT3D(data = x, diversity = "TD",
#          q = c(0, 1, 2), datatype = "abundance", nboot = 100)
#})
#save(out.rawAgg_u, file="out.rawAgg_u.RData")


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
  return(DataInfo_Table)
}
  
DataInfoPA_Table <- functionDataInfo(lAllPA, type="PA")
DataInfoPA_Table <- within(DataInfoPA_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))
DataInfoPA_d_Table <- functionDataInfo(dtPA_d, type="PA")
DataInfoPA_d_Table <- within(DataInfoPA_d_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))
DataInfoPA_u_Table <- functionDataInfo(dtPA_u, type="PA")
DataInfoPA_u_Table <- within(DataInfoPA_u_Table, rm(Assemblage, Q1, Q2, Q3, Q4, Q5))


DataInfoAgg_Table <- functionDataInfo(lAllAbuAgg, type="Abu")
DataInfoAgg_Table <- within(DataInfoAgg_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V2"] <- "Doubletons"


DataInfoAgg_Table <- functionDataInfo(lAllAbuAgg, type="Abu")
DataInfoAgg_Table <- within(DataInfoAgg_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_Table)[names(DataInfoAgg_Table)=="V2"] <- "Doubletons"

DataInfoAgg_d_Table <- functionDataInfo(dtAgg_d, type="Abu")
DataInfoAgg_d_Table <- within(DataInfoAgg_d_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_d_Table)[names(DataInfoAgg_d_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_d_Table)[names(DataInfoAgg_d_Table)=="V2"] <- "Doubletons"


DataInfoAgg_u_Table <- functionDataInfo(dtAgg_u, type="Abu")
DataInfoAgg_u_Table <- within(DataInfoAgg_u_Table, rm(Assemblage, f1,f2,f3,f4,f5))
names(DataInfoAgg_u_Table)[names(DataInfoAgg_u_Table)=="V1"] <- "Singletons"
names(DataInfoAgg_u_Table)[names(DataInfoAgg_u_Table)=="V2"] <- "Doubletons"


# save: ========================================================================
table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"
write.csv(DataInfoPA_Table, file=paste0(table_path, "/DataInfoPA_Table.csv"), row.names = F)
write.csv(DataInfoPA_d_Table, file=paste0(table_path, "/DataInfoPA_d_Table.csv"), row.names = F)
write.csv(DataInfoPA_u_Table, file=paste0(table_path, "/DataInfoPA_u_Table.csv"), row.names = F)


write.csv(DataInfoAgg_Table, file=paste0(table_path, "/DataInfoAgg_Table.csv"), row.names = F)
write.csv(DataInfoAgg_d_Table, file=paste0(table_path, "/DataInfoAgg_d_Table.csv"), row.names = F)
write.csv(DataInfoAgg_u_Table, file=paste0(table_path, "/DataInfoAgg_u_Table.csv"), row.names = F)


# palettes: ====================================================================

# Invertebrates:
c("#ffc906", "#cca105", "#997904", "#665002", "#332801", "#000000")
# Diatoms:
c("#d4ff9f", "#aacc7f", "#7f995f", "#556640", "#2a3320", "#000000")
# Fish:
c("#22f4ff", "#1bc3cc", "#149299", "#0e6266", "#073133", "#000000")



################################################################################
# plots: =======================================================================
load("out.rawPA.RData")
load("out.rawPA_d.RData")
load("out.rawPA_u.RData")

load("out.rawAgg.RData")
load("out.rawAgg_d.RData")
load("out.rawAgg_u.RData")


# panel of 3 x 4 
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"

extract_size_based <- function(dt) {
  result <- lapply(dt, function(x) {x$TDiNextEst$size_based})
  result <- as.data.frame(do.call(rbind, result))
  result$Taxa <- str_split_fixed(rownames(result), "\\.", 2)[,1]
  rownames(result) <- 1:nrow(result)
  result$Order.q <- as.factor(result$Order.q)
  result$Taxa <- recode_factor(result$Taxa, "BenthicInvertebrates" = "Benthic Invertebrates",
                               "Fish&Crustaceans" = "Fish & Crustaceans")
  result$Taxa <- factor(result$Taxa,
                  levels = c("Benthic Invertebrates", "Diatoms", "Fish", "Fish & Crustaceans"))
  
  return(result)
}
out.rawPA_plot <- extract_size_based(out.rawPA)
out.rawPA_d_plot <- extract_size_based(out.rawPA_d)
out.rawPA_u_plot <- extract_size_based(out.rawPA_u)

identical(names(out.rawPA_plot), names(out.rawPA_d_plot)) # TRUE

out.rawAgg_plot <- extract_size_based(out.rawAgg)
out.rawAgg_d_plot <- extract_size_based(out.rawAgg_d)
out.rawAgg_u_plot <- extract_size_based(out.rawAgg_u)


# Size vs coverage =============================================================

size_vs_coverage <- function(dt, type=NULL){
  plotPA <- dt
  
  dt_plot <- subset(plotPA, plotPA$Order.q=="0")
  dt_plot_ref <- subset(dt_plot, dt_plot$Method=="Observed")
  dt_plot_r <- subset(dt_plot, dt_plot$Method=="Rarefaction")
  dt_plot_e <- subset(dt_plot, dt_plot$Method=="Extrapolation")
  
  if ("nT" %in% names(plotPA)) {
  p <- ggplot(data=dt_plot, aes(x=nT, y=SC, group=Assemblage, color=Taxa)) +
    labs(title=type)+
    geom_point(data=dt_plot_ref, aes(x=nT, y=SC, group=Assemblage))+
    geom_line(data=dt_plot_r, aes(x=nT, y=SC, group=Assemblage), linetype="solid", alpha=0.4)+
    geom_line(data=dt_plot_e, aes(x=nT, y=SC, group=Assemblage), linetype="dotdash")+
    scale_color_manual(values=c("#cca105", "#aacc7f", "#149299", "#0e6266"))+
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
    facet_wrap(~Taxa, nrow=1)}else{
      p <- ggplot(data=dt_plot, aes(x=m, y=SC, group=Assemblage, color=Taxa)) +
        labs(title=type)+
        geom_point(data=dt_plot_ref, aes(x=m, y=SC, group=Assemblage))+
        geom_line(data=dt_plot_r, aes(x=m, y=SC, group=Assemblage), linetype="solid", alpha=0.4)+
        geom_line(data=dt_plot_e, aes(x=m, y=SC, group=Assemblage), linetype="dotdash")+
        scale_color_manual(values=c("#cca105", "#aacc7f", "#149299", "#0e6266"))+
        theme_classic()+
        theme(legend.position = "none",
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_text(size = 12),
              axis.text.x = element_text(size = 10),
              axis.text.y = element_text(size = 10))+
        geom_hline(yintercept=0.8, linetype=2, color = "black")+
        facet_wrap(~Taxa, nrow=1)  
    }
  return(p)
}
  

pAll <- size_vs_coverage(out.rawPA_plot, type="a) 16 sites")
pAll
pD <- size_vs_coverage(out.rawPA_d_plot, type="b) 8 sites (disturbed)")
pD
pU <- size_vs_coverage(out.rawPA_u_plot, type="c) 8 sites (undisturbed)")
pU

result <- pAll/pD/pU
gt <- patchwork::patchworkGrob(result)
gt2 <- gridExtra::grid.arrange(gt, left = "Coverage", bottom = "Size")
ggsave(paste0(plot_path, "/Size_vs_Coverage.png"), gt2, width = 8, height = 10)
  

pAggAll <- size_vs_coverage(out.rawAgg_plot, type="a) 16 sites")
pAggAll
pAggD <- size_vs_coverage(out.rawAgg_d_plot, type="b) 8 sites (disturbed)")
pAggD
pAggU <- size_vs_coverage(out.rawAgg_u_plot, type="c) 8 sites (undisturbed)")
pAggU

result <- pAggAll/pAggD/pAggU
gt <- patchwork::patchworkGrob(result)
gt2 <- gridExtra::grid.arrange(gt, left = "Coverage", bottom = "Size (number of individuals)")
ggsave(paste0(plot_path, "/Size_vs_Coverage_Aggregated.png"), gt2, width = 8, height = 10)


# Diversity vs coverage ========================================================
diversity_vs_coverage <- function(dt, order=NULL) {
  plotPA <- dt
  
  dt_plot <- subset(plotPA, plotPA$Order.q==order)
  dt_plot_ref <- subset(dt_plot, dt_plot$Method=="Observed")
  dt_plot_r <- subset(dt_plot, dt_plot$Method=="Rarefaction")
  dt_plot_e <- subset(dt_plot, dt_plot$Method=="Extrapolation")
  
  p <- ggplot(data=dt_plot, aes(x=SC, y=qTD, group=Assemblage, color=Taxa)) +
    labs(y=paste0("Order q", order))+
    geom_point(data=dt_plot_ref, aes(x=SC, y=qTD, group=Assemblage))+
    geom_line(data=dt_plot_r, aes(x=SC, y=qTD, group=Assemblage), linetype="solid", alpha=0.4)+
    geom_line(data=dt_plot_e, aes(x=SC, y=qTD, group=Assemblage), linetype="dotdash")+
    scale_color_manual(values=c("#cca105", "#aacc7f", "#149299", "#0e6266"))+
    theme_classic()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))+
    geom_vline(xintercept=0.8, linetype=2, color = "black")+
    facet_wrap(~Taxa, nrow=1)
  if ("nT" %in% names(plotPA)) {
    p <- p + scale_x_continuous(breaks = c(0.5,0.8,1),
      limits = c(0.4, 1)) +
      scale_y_continuous(
      limits = c(5, 30))
  }else{
    p <- p
  }
    
  return(p)
  
}

layout_fun <- function(x,y,z){
  result <- x/y/z
  gt <- patchwork::patchworkGrob(result)
  gt2 <- gridExtra::grid.arrange(gt, bottom = "Coverage")
  return(gt2)
}

# Order q=0 (richness): --------------------------------------------------------
pTD0 <- diversity_vs_coverage(out.rawPA_plot, order="0")
pTD0
pTD0_d <- diversity_vs_coverage(out.rawPA_d_plot, order="0")
pTD0_d
pTD0_u <- diversity_vs_coverage(out.rawPA_u_plot, order="0")
pTD0_u

pTDAgg0 <- diversity_vs_coverage(out.rawAgg_plot, order="0")
pTDAgg0
pTDAgg0_d <- diversity_vs_coverage(out.rawAgg_d_plot, order="0")
pTDAgg0_d
pTDAgg0_u <- diversity_vs_coverage(out.rawAgg_u_plot, order="0")
pTDAgg0_u

# Order q=1 (shannon): ---------------------------------------------------------
pTD1 <- diversity_vs_coverage(out.rawPA_plot, order="1")
pTD1
pTD1_d <- diversity_vs_coverage(out.rawPA_d_plot, order="1")
pTD1_d
pTD1_u <- diversity_vs_coverage(out.rawPA_u_plot, order="1")
pTD1_u

pTDAgg1 <- diversity_vs_coverage(out.rawAgg_plot, order="1")
pTDAgg1
pTDAgg1_d <- diversity_vs_coverage(out.rawAgg_d_plot, order="1")
pTDAgg1_d
pTDAgg1_u <- diversity_vs_coverage(out.rawAgg_u_plot, order="1")
pTDAgg1_u

# Order q=2 (simpson): ---------------------------------------------------------
pTD2 <- diversity_vs_coverage(out.rawPA_plot, order="2")
pTD2
pTD2_d <- diversity_vs_coverage(out.rawPA_d_plot, order="2")
pTD2_d
pTD2_u <- diversity_vs_coverage(out.rawPA_u_plot, order="2")
pTD2_u

pTDAgg2 <- diversity_vs_coverage(out.rawAgg_plot, order="2")
pTDAgg2
pTDAgg2_d <- diversity_vs_coverage(out.rawAgg_d_plot, order="2")
pTDAgg2_d
pTDAgg2_u <- diversity_vs_coverage(out.rawAgg_u_plot, order="2")
pTDAgg2_u


resultTD <- layout_fun(pTD0,pTD1,pTD2)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_AllSites.png"), resultTD, width = 8, height = 10)

resultTD_d <- layout_fun(pTD0_d,pTD1_d,pTD2_d)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_Distubed.png"), resultTD_d, width = 8, height = 10)

resultTD_u <- layout_fun(pTD0_u,pTD1_u,pTD2_u)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_Undistubed.png"), resultTD_u, width = 8, height = 10)


resultAggTD <- layout_fun(pTDAgg0,pTDAgg1,pTDAgg2)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_AggAllSites.png"), resultAggTD, width = 8, height = 10)

resultAggTD_d <- layout_fun(pTDAgg0_d,pTDAgg1_d,pTDAgg2_d)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_AggDistubed.png"), resultAggTD_d, width = 8, height = 10)

resultAggTD_u <- layout_fun(pTDAgg0_u,pTDAgg1_u,pTDAgg2_u)
ggsave(paste0(plot_path, "/Coverage_vs_Diversity_AggUndistubed.png"), resultAggTD_u, width = 8, height = 10)


################################################################################
# curves (pooling x 5 y) =======================================================
# TBC (following recommendation by AC)

# End of script#################################################################