################################################################################
# S0 Data Formatting for various analyses
# AFE
# Jan 2024 & update May 2024
################################################################################

################################################################################
# libraries: ===================================================================
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)

rm(list=ls())


formatted_data_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods/FormattedData_Trinidad"


################################################################################
# read data: ===================================================================

load(paste0(formatted_data_path, "/aggFishMacroInv2024.RData"))
load(paste0(formatted_data_path, "/BenthicInv1.RData"))
load(paste0(formatted_data_path, "/Diatoms.RData"))


f <- aggFishMacroInv2024                         # Fish & Macro-invertebrates 
bi <- BenthicInv1                                # Benthic Invertebrates 
dia <- Diatoms                                   # Diatom morphospecies


################################################################################
# format for analyses: =========================================================

c1 <- f %>% group_by(Site, Month, Year) %>% summarise(nDay= n_distinct(Day))
range(c1$nDay) # 1, OK 
range(f$Year) # 2011-2024

sort(unique(f$Genus))
sort(unique(f$Species))
sort(unique(bi$Classification))
sort(unique(dia$taxa))

names(f)
names(bi)
names(dia)

f$Taxa <- paste0(f$Genus, "_", f$Species)
fa <- f[,names(f) %in% c("Site", "Session", "Taxa", "Abundance")]

names(bi)<-c("Session", "Site", "Taxa", "Abundance")
names(dia)<-c("Session", "Site", "Taxa", "Abundance")

bi$Taxa <- gsub(" ", "_", bi$Taxa)
sort(unique(bi$Taxa))

# subsets sensitivity: ---------------------------------------------------------
unique(f$Taxa)
fsen <- subset(fa, !fa$Taxa %in% c("Atya_spp", "Macrobrachium_spp", "Eudaniela_garmani"))
length(unique(fsen$Taxa)) # 23, OK

fcru <- subset(fa, fa$Taxa %in% c("Atya_spp", "Macrobrachium_spp", "Eudaniela_garmani"))
fcru <- subset(fcru, fcru$Session < 20)

unique(fcru$Taxa)
fcru$Taxa[fcru$Taxa=="Eudaniela_garmani"] <- "crustacea_decapoda_pseudothelphusidae"
fcru$Taxa[fcru$Taxa=="Macrobrachium_spp"] <- "crustacea_spp"
fcru$Taxa[fcru$Taxa=="Atya_spp"] <- "crustacea_spp"

bisen <- as.data.frame(rbind(bi, fcru))
bisen <- bisen %>% group_by(Site, Session, Taxa) %>% summarise(Abundance=sum(Abundance)) # grouped

# Biomass: ---------------------------------------------------------------------
fb <- f[,names(f) %in% c("Site", "Session", "Taxa", "Biomass")]
fbsen <- subset(fb, !fb$Taxa %in% c("Atya_spp", "Macrobrachium_spp", "Eudaniela_garmani"))


################################################################################
# incidence, abundance & agg. abundance: =======================================

format_incidence <- function(x, type=NULL, v="NULL") {
  dt <- x
  if (v=="abundance"){
    dt_spread <- spread(dt, key="Taxa", value="Abundance", fill = 0)  
  }else{
    dt_spread <- spread(dt, key="Taxa", value="Biomass", fill = 0)
  }
  dt_spread$Site <- gsub(" ", "_", dt_spread$Site)
  vecCols <- names(dt_spread)[c(3:ncol(dt_spread))]
  if (type=="PA") {
    dt_spread[vecCols] <- apply(dt_spread[vecCols], 2, function(x) {x[x>0] <-1;x})
  }else{
    dt_spread <- dt_spread
  }
  dt_list <- split(dt_spread, f=dt_spread$Session)
  dt_list <- lapply(dt_list, function(x) {within(x, rm(Session))})
  dt_list <- lapply(dt_list, function(x) {as.data.frame(t(x))})
  dt_list <- lapply(dt_list, function(x) {names(x) <- x[1,]; x
  x <- x[-1,]; x})
  sps <- lapply(dt_list, function(x) {x <- rownames(x)})
  dt_list <- lapply(dt_list, function(x) {apply(x, 2, as.numeric)})
  dt_listII <- mapply(function(x,y) {rownames(x) <- y; x}, dt_list, sps, SIMPLIFY=F)
  dt_listII <- lapply(dt_listII, as.data.frame)
  return(dt_listII)
}

lAll <- list("fishmac"=as.data.frame(fa), "fish"=as.data.frame(fsen),
             "bi"=as.data.frame(bi), "dia"=as.data.frame(dia), "bimac"=as.data.frame(bisen)) # main list
lBio <- list("fishmacb"=as.data.frame(fb), "fishb"=as.data.frame(fbsen))

lAllPA <- lapply(lAll, function(x) {format_incidence(x, type="PA", v="abundance")})      # Presence-absence (Metacommunity)
lAllAbu <- lapply(lAll, function(x) {format_incidence(x, type="NULL",  v="abundance")})  # Site abundances  (Community)
lAllBio <- lapply(lBio, function(x) {format_incidence(x, type="NULL",  v="biomass")})    # Site biomasses (Fish Community)

lAllAbuAgg <- lapply(lAllAbu, function(x) {lapply(x, function(y) {rowSums(y)})})
lAllBioAgg <- lapply(lAllBio, function(x) {lapply(x, function(y) {rowSums(y)})})


################################################################################
# create lists D & U: ==========================================================
lPA_d <- lapply(lAllPA, function(x) {lapply(x, function(y) {y[names(y) %like% "_d"]})})
lPA_u <- lapply(lAllPA, function(x) {lapply(x, function(y) {y[names(y) %like% "_u"]})})

lAbu_d <- lapply(lAllAbu, function(x) {lapply(x, function(y) {y[names(y) %like% "_d"]})})
lAbuAgg_d <- lapply(lAbu_d, function(x) {lapply(x, function(y) {rowSums(y)})})

lAbu_u <- lapply(lAllAbu, function(x) {lapply(x, function(y) {y[names(y) %like% "_u"]})})
lAbuAgg_u <- lapply(lAbu_u, function(x) {lapply(x, function(y) {rowSums(y)})})

lBio_d <- lapply(lAllBio, function(x) {lapply(x, function(y) {y[names(y) %like% "_d"]})})
lBioAgg_d <- lapply(lBio_d, function(x) {lapply(x, function(y) {rowSums(y)})})

lBio_u <- lapply(lAllBio, function(x) {lapply(x, function(y) {y[names(y) %like% "_u"]})})
lBioAgg_u <- lapply(lBio_u, function(x) {lapply(x, function(y) {rowSums(y)})})


################################################################################
# rm macrocrustaceans: =========================================================
rm_fishmac_bimac <- function(x){
  x[["fishmac"]] <- NULL
  x[["bimac"]] <- NULL
  return(x)
}

lAllPA <- rm_fishmac_bimac(lAllPA)
lAllAbu <- rm_fishmac_bimac(lAllAbu)
lAllBio <- rm_fishmac_bimac(lAllBio)
lAllAbuAgg <- rm_fishmac_bimac(lAllAbuAgg)
lAllBioAgg <- rm_fishmac_bimac(lAllBioAgg)

lPA_d <- rm_fishmac_bimac(lPA_d)
lAbu_d <- rm_fishmac_bimac(lAbu_d)
lBio_d <- rm_fishmac_bimac(lBio_d)
lAbuAgg_d <- rm_fishmac_bimac(lAbuAgg_d)
lBioAgg_d <- rm_fishmac_bimac(lBioAgg_d)

lPA_u <- rm_fishmac_bimac(lPA_u)
lAbu_u <- rm_fishmac_bimac(lAbu_u)
lBio_u <- rm_fishmac_bimac(lBio_u)
lAbuAgg_u <- rm_fishmac_bimac(lAbuAgg_u)
lBioAgg_u <- rm_fishmac_bimac(lBioAgg_u)


################################################################################
# save objects: ================================================================
mainanalysisRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/mainanalysisRData"
setwd(mainanalysisRData)

# NR subsets: ------------------------------------------------------------------
save(lAllPA, file="lAllPA.RData")            # pa data (main)
save(lAllAbu, file="lAllAbu.RData")          # site abundance 
save(lAllBio, file="lAllBio.RData")          # site biomass
save(lAllAbuAgg, file="lAllAbuAgg.RData")    # agg abundance
save(lAllBioAgg, file="lAllBioAgg.RData")    # agg biomass


# Disturbed subsets: -----------------------------------------------------------
save(lPA_d, file="lPA_d.RData")            # pa data (main D)
save(lAbu_d, file="lAbu_d.RData")          # site abundance 
save(lBio_d, file="lBio_d.RData")          # site biomass
save(lAbuAgg_d, file="lAbuAgg_d.RData")    # agg abundance
save(lBioAgg_d, file="lBioAgg_d.RData")    # agg biomass

# Undisturbed subsets: ---------------------------------------------------------
save(lPA_u, file="lPA_u.RData")            # pa data (main U)
save(lAbu_u, file="lAbu_u.RData")          # site abundance 
save(lBio_u, file="lBio_u.RData")          # site biomass
save(lAbuAgg_u, file="lAbuAgg_u.RData")    # agg abundance
save(lBioAgg_u, file="lBioAgg_u.RData")    # agg biomass

setwd("C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs")


# END OF SCRIPT ################################################################
################################################################################