################################################################################
# S3 Compute Diversity Chao Numbers
# AFE
# January 2023
################################################################################


# Libraries: ===================================================================
library(rlang)
library(tidyverse)
library(ggcorset)
library(patchwork)
library(ape)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gg.gap)
library(data.table)
library(plotly)
library(reshape2)
library(lmtest)
library(future.apply)
library(abind)
library(parallel)
library(geosphere)
library(openxlsx)
library(dplyr)
library(corrplot)
library(gridExtra)
library(ggplotify)
library(grid)

library(devtools)
library(iNEXT.3D)
library(iNEXT.beta3D)


rm(list=ls())

################################################################################
# Load Data ====================================================================
mainanalysisRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/mainanalysisRData"

all_files <- list.files(mainanalysisRData, pattern = "\\.RData$", full.names = TRUE)

# Loop through the list and load each RData file
for (file in all_files) {
  load(file)
}

names(lAllPA) <- c("Fish", "Invertebrates", "Diatoms")
names(lPA_d) <- c("Fish", "Invertebrates", "Diatoms")
names(lPA_u) <- c("Fish", "Invertebrates", "Diatoms")

names(lAllAbuAgg) <- c("Fish", "Invertebrates", "Diatoms")
names(lAbuAgg_d) <- c("Fish", "Invertebrates", "Diatoms")
names(lAbuAgg_u) <- c("Fish", "Invertebrates", "Diatoms")


#load("dist_matS.RData")
#dist_matS <- dist_matS[-19, -19]
#acr <- read.csv("acronymsNR.csv", h=T)
#acr$acronym[acr$acronym=="Osp"] <- "Omossambicus"
#rownames(dist_matS) <- acr$name[match(rownames(dist_matS), acr$acronym)]
#colnames(dist_matS) <- acr$name[match(colnames(dist_matS), acr$acronym)]

Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
Results_BetaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_BetaDivRData"

################################################################################
# Diversity:====================================================================


# Alpha PA =====================================================================
# NOTE: Had to re-compute with iNext3DBeta / otherwise no Cmin & Cmax across sampling events

# Presence-absence: ------------------------------------------------------------

lAllPAalpha <- lapply(lAllPA, function(x){lapply(x, function(y) {list(as.matrix(y), as.matrix(y))})})
lPA_dalpha <- lapply(lPA_d , function(x){lapply(x, function(y) {list(as.matrix(y), as.matrix(y))})})
lPA_ualpha  <- lapply(lPA_u , function(x){lapply(x, function(y) {list(as.matrix(y), as.matrix(y))})})

lAllPATD <- lapply(lAllPAalpha, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "incidence_raw",
                                        base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Northern Range")], nboot = 100, conf = 0.95)})
lPA_dTD  <- lapply(lPA_dalpha, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "incidence_raw",
                                                  base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Disturbed")], nboot = 100, conf = 0.95)})
lPA_uTD  <- lapply(lPA_ualpha, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "incidence_raw",
                                                    base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Undisturbed")], nboot = 100, conf = 0.95)})
save(lAllPATD, file = paste0(Results_AlphaDivRData, "/lAllPATD.RData"))
save(lPA_dTD, file = paste0(Results_AlphaDivRData, "/lPA_dTD.RData"))
save(lPA_uTD, file = paste0(Results_AlphaDivRData, "/lPA_uTD.RData"))


lAllPAFD <- iNEXTbeta3D(lAllPAalpha$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "incidence_raw", FDtype = "tau_value",
            base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Northern Range")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
lPA_dFD <- iNEXTbeta3D(lPA_dalpha$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "incidence_raw", FDtype = "tau_value",
                        base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Disturbed")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
lPA_uFD <- iNEXTbeta3D(lPA_ualpha$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "incidence_raw", FDtype = "tau_value",
                       base = "coverage", level = lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Undisturbed")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
save(lAllPAFD, file = paste0(Results_AlphaDivRData, "/lAllPAFD.RData"))
save(lPA_dFD, file = paste0(Results_AlphaDivRData, "/lPA_dFD.RData"))
save(lPA_uFD, file = paste0(Results_AlphaDivRData, "/lPA_uFD.RData"))


# Format PA: -------------------------------------------------------------------
retrieveAlphaTD <- function(dt){
  dt <- lapply(dt, function(x) {lapply(x, function(y) {y$alpha})})
  dt <- lapply(dt, function(x) {as.data.frame(do.call(rbind, x))})
  dt <- Map(cbind, dt, "Taxa"=names(dt))
  dt <- as.data.frame(do.call(rbind, dt))
  names(dt)[names(dt)=="Dataset"] <- "Session"
  return(dt)
}

retrieveAlphaFD <- function(dt){
  dt <- lapply(dt, function(x) {x$alpha})
  dt <- as.data.frame(do.call(rbind, dt))
  names(dt)[names(dt)=="Dataset"] <- "Session"
  return(dt)
}

lAllPATD <- retrieveAlphaTD(lAllPATD)
lAllPATD$Group <- rep("Northern Range", nrow(lAllPATD))
lPA_dTD <- retrieveAlphaTD(lPA_dTD)
lPA_dTD$Group <- rep("Disturbed", nrow(lPA_dTD))
lPA_uTD <- retrieveAlphaTD(lPA_uTD)
lPA_uTD$Group <- rep("Undisturbed", nrow(lPA_uTD))

lAllPAFD <- retrieveAlphaFD(lAllPAFD)
lAllPAFD$Group <- rep("Northern Range", nrow(lAllPAFD))
lPA_dFD <- retrieveAlphaFD(lPA_dFD)
lPA_dFD$Group <- rep("Disturbed", nrow(lPA_dFD))
lPA_uFD <- retrieveAlphaFD(lPA_uFD)
lPA_uFD$Group <- rep("Undisturbed", nrow(lPA_uFD))
lPAFD <- rbind(lAllPAFD, lPA_dFD, lPA_uFD)
lPAFD$Taxa <- "Fish (Trait)"

lPA_Alpha <- bind_rows(lAllPATD, lPA_dTD, lPA_uTD, lPAFD)
lPA_Alpha$SC2 <- rep(NA, nrow(lPA_Alpha))
lPA_Alpha$Type <- rep("PA", nrow(lPA_Alpha))

lcov

lPA_Alpha$SC2[lPA_Alpha$SC==0.9387 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Northern Range"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9547 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Northern Range"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9940 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Northern Range"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9980 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Northern Range"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9785 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Northern Range"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9893 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Northern Range"]<- "Cmax"


lPA_Alpha$SC2[lPA_Alpha$SC==0.8995 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Disturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9307 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Disturbed"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9745 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Disturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9909 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Disturbed"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9599 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Disturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9907 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Disturbed"]<- "Cmax"


lPA_Alpha$SC2[lPA_Alpha$SC==0.8515 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Undisturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.8977 & lPA_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lPA_Alpha$Group =="Undisturbed"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9731 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Undisturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9870 & lPA_Alpha$Taxa %in% c("Invertebrates") & lPA_Alpha$Group =="Undisturbed"]<- "Cmax"

lPA_Alpha$SC2[lPA_Alpha$SC==0.9581 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Undisturbed"]<- "Cmin"
lPA_Alpha$SC2[lPA_Alpha$SC==0.9732 & lPA_Alpha$Taxa %in% c("Diatoms") & lPA_Alpha$Group =="Undisturbed"]<- "Cmax"


#lPA_Alpha$Concat <- paste0(lPA_Alpha$Session, "_", lPA_Alpha$Taxa, "_",
#                           lPA_Alpha$Group, "_", lPA_Alpha$Type, "_", lPA_Alpha$SC)

#lPA_Alpha$SC3 <- for_observed$SC[match(lPA_Alpha$Concat, for_observed$Concat)] # no matches

lPA_Alpha <- lPA_Alpha[!is.na(lPA_Alpha$SC2),]
save(lPA_Alpha, file = paste0(Results_AlphaDivRData, "/lPA_Alpha.RData"))



# Alpha Aggregated =============================================================
# Aggregated: ------------------------------------------------------------------

lAllAbuAgg <- lapply(lAllAbuAgg, function(x){lapply(x, function(y) {as.matrix(y)})})
lAllAbuAggTD <- lapply(lAllAbuAgg, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "abundance",
                                                  base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Northern Range")], nboot = 0, conf = 0.95)})
lAbuAgg_d <- lapply(lAbuAgg_d, function(x){lapply(x, function(y) {as.matrix(y)})})
lAbuAgg_dTD <- lapply(lAbuAgg_d, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "abundance",
                                                          base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Disturbed")], nboot = 0, conf = 0.95)})
lAbuAgg_u <- lapply(lAbuAgg_u, function(x){lapply(x, function(y) {as.matrix(y)})})
lAbuAgg_uTD <- lapply(lAbuAgg_u, function(x) {iNEXTbeta3D(x, diversity = "TD", q = c(0, 1, 2), datatype = "abundance",
                                                          base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Undisturbed")], nboot = 0, conf = 0.95)})


lAllAbuAggTD <- retrieveAlphaTD(lAllAbuAggTD)
lAllAbuAggTD$Group <- rep("Northern Range", nrow(lAllAbuAggTD))
lAbuAgg_dTD <- retrieveAlphaTD(lAbuAgg_dTD)
lAbuAgg_dTD$Group <- rep("Disturbed", nrow(lAbuAgg_dTD))
lAbuAgg_uTD <- retrieveAlphaTD(lAbuAgg_uTD)
lAbuAgg_uTD$Group <- rep("Undisturbed", nrow(lAbuAgg_uTD))

################################################################################
#lAllAbuAggFD <- iNEXTbeta3D(lAllAbuAgg$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "abundance", FDtype = "tau_value",
#                        base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Northern Range")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
#lAbuAgg_dFD <- iNEXTbeta3D(lAbuAgg_d$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "abundance", FDtype = "tau_value",
#                       base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Disturbed")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
#lAbuAgg_uFD <- iNEXTbeta3D(lAbuAgg_u$Fish, diversity = "FD", q = c(0, 1, 2), datatype = "abundance", FDtype = "tau_value",
#                       base = "coverage", level = lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Undisturbed")], nboot = 0, conf = 0.95, FDdistM = dist_matS)
# DOES NOT RUN!!!!!!!!!!!!!!

#lAllAbuAggFD <- retrieveAlphaFD(lAllAbuAggFD)
#lAllAbuAggFD$Group <- rep("Northern Range", nrow(lAllAbuAggFD))

#lAbuAgg_dFD <- retrieveAlphaFD(lAbuAgg_dFD)
#lAbuAgg_dFD$Group <- rep("Disturbed", nrow(lAbuAgg_dFD))

#lAbuAgg_uFD <- retrieveAlphaFD(lAbuAgg_uFD)
#lAbuAgg_uFD$Group <- rep("Undisturbed", nrow(lAbuAgg_uFD))

#lAbuAgg_FD <- rbind(lAllAbuAggFD, lAbuAgg_dFD, lAbuAgg_uFD)
#lAbuAgg_FD$Taxa <- "Fish (Trait)"
################################################################################

lAgg_Alpha <- bind_rows(lAllAbuAggTD, lAbuAgg_dTD, lAbuAgg_uTD) 
lAgg_Alpha$SC2 <- rep(NA, nrow(lAgg_Alpha))
lAgg_Alpha$Type <- rep("Agg", nrow(lAgg_Alpha))

lcov

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9985 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Northern Range"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9990 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Northern Range"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9998 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Northern Range"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==1.0000 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Northern Range"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9998 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Northern Range"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9999 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Northern Range"]<- "Cmax"


lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9984 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Disturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9993 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Disturbed"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9989 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Disturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9996 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Disturbed"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9993 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Disturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9997 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Disturbed"]<- "Cmax"


lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9967 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9983 & lAgg_Alpha$Taxa %in% c("Fish", "Fish (Trait)") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9992 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9997 & lAgg_Alpha$Taxa %in% c("Invertebrates") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmax"

lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9996 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmin"
lAgg_Alpha$SC2[lAgg_Alpha$SC==0.9999 & lAgg_Alpha$Taxa %in% c("Diatoms") & lAgg_Alpha$Group =="Undisturbed"]<- "Cmax"

lAgg_Alpha <- lAgg_Alpha[!is.na(lAgg_Alpha$SC2),]   
save(lAgg_Alpha, file = paste0(Results_AlphaDivRData, "/lAgg_Alpha.RData"))



# Final Alpha ==================================================================
lAlpha <- as.data.frame(bind_rows(lPA_Alpha, lAgg_Alpha))
save(lAlpha, file=paste0(Results_AlphaDivRData, "/lAlpha.RData"))



# Beta PA ======================================================================
combosfunA <- function(x, y){
  list_of_lists <- unlist(lapply(x, function(a) lapply(y, function (b) list(a, b))), recursive=FALSE) # needed structure to compute lag beta diversity
  return(list_of_lists)
}
combosfunB <- function(dt,dt2){
  l_1 <- dt
  l_2 <- dt2 # correct, same
  mIll <- mapply(function(x,y){combosfunA(x,y)}, l_1, l_2, SIMPLIFY = F)                                # keep the same general structure
  mIll <- lapply(mIll, function(x) {lapply(x, function(y){lapply(y, function(z){as.matrix(z)})})})      # ready for iNEXT beta
  return(mIll)
}

mPA <- combosfunB(lAllPA, lAllPA)
mPA_d <- combosfunB(lPA_d, lPA_d)
mPA_u <- combosfunB(lPA_u, lPA_u)
mPA_du <- combosfunB(lPA_d, lPA_u) 


mPA_FD <- list(mPA$Fish, mPA$Fish)
mPA_dFD <- list(mPA_d$Fish, mPA_d$Fish)
mPA_uFD <- list(mPA_u$Fish, mPA_u$Fish)
mPA_duFD <- list(mPA_du$Fish, mPA_du$Fish)


# Taxonomic: -------------------------------------------------------------------
#beta.out.rawPA <- lapply(mPA, function(x) {iNEXTbeta3D(data = x, diversity = "TD", datatype = "incidence_raw",
#                                           base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Northern Range")])})              
#beta.out.rawPA_d <- lapply(mPA_d, function(x) {iNEXTbeta3D(data = x, diversity = "TD", datatype = "incidence_raw",
#                                                        base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Disturbed")])})               
#beta.out.rawPA_u <- lapply(mPA_u, function(x) {iNEXTbeta3D(data = x, diversity = "TD", datatype = "incidence_raw",
#                                                        base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Undisturbed")])})    
#beta.out.rawPA_du <- lapply(mPA_du, function(x) {iNEXTbeta3D(data = x, diversity = "TD", datatype = "incidence_raw",
#                                                        base = "coverage", nboot = 0, level=lcov$CoverageVals)}) 
#save(beta.out.rawPA, beta.out.rawPA_d, beta.out.rawPA_u, beta.out.rawPA_du, file = paste0(Results_BetaDivRData, "/beta.out.rawPATaxonomic.RData"))

# Fish Trait:-------------------------------------------------------------------
#mPA <- list(mPA$fish, mPA$fish)
#beta.out.rawPAFD <- lapply(mPA, function(x) {iNEXTbeta3D(data = x, diversity = "FD", datatype = "incidence_raw",
#base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Northern Range")], FDdistM=dist_matS)})              
#identical(beta.out.rawPAFD[[1]], beta.out.rawPAFD[[2]]) # TRUE

#beta.out.rawPA_FD <- lapply(mPA_FD, function(x) {iNEXTbeta3D(data = x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_value",
#                                                       base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Northern Range")], FDdistM=dist_matS)})              
#beta.out.rawPA_dFD <- lapply(mPA_dFD, function(x) {iNEXTbeta3D(data = x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_value",
#                                                           base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Disturbed")], FDdistM=dist_matS)})               
#beta.out.rawPA_uFD <- lapply(mPA_uFD, function(x) {iNEXTbeta3D(data = x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_value",
#                                                          base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="PA" & lcov$Group=="Undisturbed")], FDdistM=dist_matS)})    
#beta.out.rawPA_duFD <- lapply(mPA_duFD, function(x) {iNEXTbeta3D(data = x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_value",
#                                                             base = "coverage", nboot = 0, level=lcov$CoverageVals, FDdistM=dist_matS)}) 
#identical(beta.out.rawPA_FD[[1]], beta.out.rawPA_FD[[2]]) # TRUE
#save(beta.out.rawPA_FD, beta.out.rawPA_dFD, beta.out.rawPA_uFD, beta.out.rawPA_duFD, file = paste0(Results_BetaDivRData, "/beta.out.rawPATrait.RData"))


################################################################################
# Taxonomic Agg: ---------------------------------------------------------------
#beta.out.rawPA <- lapply(mAgg, function(x) {iNEXTbeta3D(data = x, diversity = "TD", datatype = "abundance",
#                                                      base = "coverage", nboot = 0, level=lcov$CoverageVals[(lcov$Type=="Agg" & lcov$Group=="Northern Range")])})              
# NOT RUNNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################################



# Format Beta PA ===============================================================
load(paste0(Results_BetaDivRData, "/beta.out.rawPATaxonomic.RData"))
load(paste0(Results_BetaDivRData, "/beta.out.rawPATrait.RData"))

beta.out.rawFormat <- function(dt){
  dtF <- lapply(dt, function(x) {lapply(x, function(y) {lapply(y, function(z){z$s1 <- as.numeric(str_split_fixed(z$Dataset, "\\.", 2)[,1]);z
                                                                              z$s2 <- as.numeric(str_split_fixed(z$Dataset, "\\.", 2)[,2]);z})})})
  dtF <- lapply(dtF, function(x) {lapply(x, function(y) {Map(cbind, y, "Metric" = names(y))})})
  dtF <- lapply(dtF, function(x) {lapply(x, function(y) {lapply(y, function(z){names(z)[5] <- "Estimate";z})})})
  dtF <- lapply(dtF, function(x) {lapply(x, function(y) {do.call(rbind, y)})})
  dtF <- lapply(dtF, function(x) {do.call(rbind, x)})
  dtF <- map_df(dtF, ~as.data.frame(.x), .id="Taxa")
  dtF <- dtF[!dtF$Metric %in% c("alpha", "gamma"),]

  return(dtF)
}

# Taxonomic: -------------------------------------------------------------------
beta.out.rawPA <- beta.out.rawFormat(beta.out.rawPA)
str(beta.out.rawPA)
beta.out.rawPA$Group <- "Northern Range"

beta.out.rawPA_d <- beta.out.rawFormat(beta.out.rawPA_d)
str(beta.out.rawPA_d)
beta.out.rawPA_d$Group <- "Disturbed"

beta.out.rawPA_u <- beta.out.rawFormat(beta.out.rawPA_u)
str(beta.out.rawPA_u)
beta.out.rawPA_u$Group <- "Undisturbed"

beta.out.rawPA_du <- beta.out.rawFormat(beta.out.rawPA_du)
str(beta.out.rawPA_du)
beta.out.rawPA_du$Group <- "DU"

beta.out.rawPATD <- as.data.frame(rbind(beta.out.rawPA, beta.out.rawPA_d, beta.out.rawPA_u, beta.out.rawPA_du))
range(beta.out.rawPATD$Estimate[beta.out.rawPATD$Metric=="beta"]) # 0.4191667 2.1909722 vals of less than 1

# Trait: -----------------------------------------------------------------------

formatFDBeta <- function(dt){
  dtFBeta <- dt[[1]]
  dtFBeta <- lapply(dtFBeta, function(y) {lapply(y, function(z){z$s1 <- as.numeric(str_split_fixed(z$Dataset, "\\.", 2)[,1]);z
  z$s2 <- as.numeric(str_split_fixed(z$Dataset, "\\.", 2)[,2]);z})})
  dtFBeta <- lapply(dtFBeta, function(y) {Map(cbind, y, "Metric" = names(y))})
  dtFBeta <- lapply(dtFBeta, function(y) {lapply(y, function(z){names(z)[5] <- "Estimate";z})})
  dtFBeta <- lapply(dtFBeta, function(y) {do.call(rbind, y)})
  dtFBeta <- do.call(rbind, dtFBeta)
  dtFBeta$Taxa <- rep("Fish (Trait)", nrow(dtFBeta))
  dtFBeta <- dtFBeta[!dtFBeta$Metric %in% c("alpha", "gamma"),]
  return(dtFBeta)
}

beta.out.rawPA_FD <- formatFDBeta(beta.out.rawPA_FD)
beta.out.rawPA_FD$Group <- "Northern Range"
beta.out.rawPA_dFD <- formatFDBeta(beta.out.rawPA_dFD)
beta.out.rawPA_dFD$Group <- "Disturbed"
beta.out.rawPA_uFD <- formatFDBeta(beta.out.rawPA_uFD)
beta.out.rawPA_uFD$Group <- "Undisturbed"
beta.out.rawPA_duFD <- formatFDBeta(beta.out.rawPA_duFD)
beta.out.rawPA_duFD$Group <- "DU"

names(beta.out.rawPA)
names(beta.out.rawPA_FD)

beta.out.rawPAFD <- as.data.frame(rbind(beta.out.rawPA_FD, beta.out.rawPA_dFD, beta.out.rawPA_uFD, beta.out.rawPA_duFD))
beta.out.rawPAFD$Taxa <- "Fish (Trait)"

range(beta.out.rawPAFD$Estimate[beta.out.rawPAFD$Metric=="beta"]) # 0.82 1.45 vals of less than 1

beta.out.rawPA <- as.data.frame(bind_rows(beta.out.rawPATD, beta.out.rawPAFD))


# Adding SC2: ------------------------------------------------------------------
beta.out.rawPA$SC2 <- rep(NA, nrow(beta.out.rawPA))

lcov

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9387 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Northern Range"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9547 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Northern Range"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9940 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Northern Range"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9980 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Northern Range"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9785 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Northern Range"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9893 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Northern Range"]<- "Cmax"


beta.out.rawPA$SC2[beta.out.rawPA$SC==0.8995 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Disturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9307 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Disturbed"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9745 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Disturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9909 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Disturbed"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9599 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Disturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9907 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Disturbed"]<- "Cmax"


beta.out.rawPA$SC2[beta.out.rawPA$SC==0.8515 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.8977 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9731 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9870 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9581 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9732 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="Undisturbed"]<- "Cmax"


beta.out.rawPA$SC2[beta.out.rawPA$SC==0.8515 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="DU"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.8977 & beta.out.rawPA$Taxa %in% c("Fish", "Fish (Trait)") & beta.out.rawPA$Group =="DU"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9731 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="DU"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9870 & beta.out.rawPA$Taxa %in% c("Invertebrates") & beta.out.rawPA$Group =="DU"]<- "Cmax"

beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9581 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="DU"]<- "Cmin"
beta.out.rawPA$SC2[beta.out.rawPA$SC==0.9732 & beta.out.rawPA$Taxa %in% c("Diatoms") & beta.out.rawPA$Group =="DU"]<- "Cmax"


beta.out.rawPA <- beta.out.rawPA[!is.na(beta.out.rawPA$SC2),]
beta.out.rawPA$Type= "PA"


# Sessions:=====================================================================
beta.out.rawPA$s1 <- as.character(beta.out.rawPA$s1)
beta.out.rawPA$s1 <- plyr::revalue(beta.out.rawPA$s1, c("20"="23",
                                                        "21"="46",
                                                        "22"="50"))
beta.out.rawPA$s1 <- as.numeric(beta.out.rawPA$s1)
beta.out.rawPA$s2 <- as.character(beta.out.rawPA$s2)
beta.out.rawPA$s2 <- plyr::revalue(beta.out.rawPA$s2, c("20"="23",
                                                        "21"="46",
                                                        "22"="50"))
beta.out.rawPA$s2 <- as.numeric(beta.out.rawPA$s2)


# Additional Formatting steps ==================================================
beta.out.rawPA <- beta.out.rawPA[beta.out.rawPA$Metric=="1-U",] # keep Jaccard
beta.out.rawPA$concat_sessions <- rep(NA, nrow(beta.out.rawPA))
beta.out.rawPA$concat_sessions <- apply(beta.out.rawPA[c("s1", "s2")], 1, function(x) paste(sort(x), collapse = "_"))
names(beta.out.rawPA)
beta.out.rawPA <- beta.out.rawPA %>%
  group_by(concat_sessions, Taxa, Order.q, Diversity, Metric, Group, SC2, Type) %>%
  slice(1)                                                      # remove dups of s1 & s2


beta.out.rawPA$Lag <- abs(beta.out.rawPA$s1-beta.out.rawPA$s2)  # add Lag
range(beta.out.rawPA$Lag)                                       # 0 49
beta.out.rawPA$LagSqrt <- sqrt(beta.out.rawPA$Lag)              # for analysis


################################################################################
# Data for analyses: ===========================================================
lBetaLag <- beta.out.rawPA[!beta.out.rawPA$Group =="DU",]
lBetaLag <- lBetaLag[!(lBetaLag$s1 == lBetaLag$s2),]
range(lBetaLag$Lag) # 1-49, OK

lBetaDU <- beta.out.rawPA[beta.out.rawPA$Group =="DU",]
lBetaDU <- lBetaDU[lBetaDU$s1 == lBetaDU$s2,]
range(lBetaDU$Lag)  # 0 , OK

save(lBetaLag, file=paste0(Results_BetaDivRData, "/lBetaLag.RData"))
save(lBetaDU, file=paste0(Results_BetaDivRData, "/lBetaDU.RData"))

# End of script ################################################################
################################################################################