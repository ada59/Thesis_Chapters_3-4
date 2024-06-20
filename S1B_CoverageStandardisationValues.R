################################################################################
# Time lag analysis (pre-Chao)
# Date: January 2024
# AFE
################################################################################


# Libraries: ===================================================================
library(tseries)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)


getwd()

# Data =========================================================================
rm(list=ls())
mainanalysisRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/mainanalysisRData"

all_files <- list.files(mainanalysisRData, pattern = "\\.RData$", full.names = TRUE)

# Loop through the list and load each RData file
for (file in all_files) {
  load(file)
}

matI <- lapply(lAllPA, function(x) {do.call(cbind, x)})     
matI_d <- lapply(lPA_d, function(x) {do.call(cbind, x)})     
matI_u <- lapply(lPA_u, function(x) {do.call(cbind, x)})     

matAgg <- lapply(lAllAbuAgg, function(x) {do.call(cbind, x)})     
matAgg_d <- lapply(lAbuAgg_d, function(x) {do.call(cbind, x)})     
matAgg_u <- lapply(lAbuAgg_u, function(x) {do.call(cbind, x)})     


# Cmin & Cmax:==================================================================

# Incidence --------------------------------------------------------------------
nTfun <- function(x,y){
  nT <- data.frame(matrix(ncol=x, nrow = 1))
  names(nT) <- c(1:x)
  nT[1,] <- rep(y, ncol(nT))
  rownames(nT) <- "nT"
  return(nT)
}
nT_23_NR <-  nTfun(23,16)
nT_19_NR <-  nTfun(19,16)
nT_23_DU <-  nTfun(23,8)
nT_19_DU <-  nTfun(19,8)


CmaxI <- function(x, nT) {
  nT <- nT
  dt <- x
  cmaxI <- sapply(1:length(nT), function(i) rowSums(dt[, (sum(nT[1:i]) - sum(nT[i]) + 1) : sum(nT[1:i])] )) %>% rbind(as.integer(nT),.) %>% 
    apply(., 2, function(x) iNEXT.3D:::Coverage(x, 'incidence_freq', 2*x[1])) %>% min %>% round(., 4) 
  return(cmaxI)
}   
CminI <- function(x, nT=n) {
  nT <- nT
  dt <- x
  cminI <- sapply(1:length(nT), function(i) rowSums( dt[, (sum(nT[1:i]) - sum(nT[i]) + 1) : sum(nT[1:i])] )) %>% rbind(as.integer(nT),.) %>% 
    apply(., 2, function(x) iNEXT.3D:::Coverage(x, 'incidence_freq', x[1])) %>% min %>% round(., 4) 
  return(cminI)
}


lnt <- list(nT_23_NR, nT_19_NR, nT_19_NR)
lmaxI <- mapply(function(x, y) {c(CmaxI(x,y))}, matI, lnt, SIMPLIFY=T)
lminI <- mapply(function(x, y) {c(CminI(x,y))}, matI, lnt, SIMPLIFY=T)
lcovI <- c(lmaxI, lminI) # all coverage values for incidence data
lcovI <- data.frame("CoverageVals"=lcovI, 
                    "Coverage"=c("Cmax", "Cmax", "Cmax", "Cmin", "Cmin", "Cmin"),
                    "Group"=rep("Northern Range", times=6),
                    "Type"=rep("PA", times=6),
                    "Taxa"=c("Fish", "Invertebrates", "Diatoms", "Fish", "Invertebrates", "Diatoms"))
range(lcovI$CoverageVals)    # 0.9387 0.9980


lntDU <- list(nT_23_DU, nT_19_DU, nT_19_DU)
lmaxI_d <- mapply(function(x, y) {c(CmaxI(x,y))}, matI_d, lntDU, SIMPLIFY=T)
lminI_d <- mapply(function(x, y) {c(CminI(x,y))}, matI_d, lntDU, SIMPLIFY=T)
lcovI_d <- c(lmaxI_d, lminI_d) # all coverage values for incidence data
lcovI_d <- data.frame("CoverageVals"=lcovI_d, 
                    "Coverage"=c("Cmax", "Cmax", "Cmax", "Cmin", "Cmin", "Cmin"),
                    "Group"=rep("Disturbed", times=6),
                    "Type"=rep("PA", times=6),
                    "Taxa"=c("Fish", "Invertebrates", "Diatoms", "Fish", "Invertebrates", "Diatoms"))
range(lcovI_d$CoverageVals)  # 0.8995 0.9909

lmaxI_u <- mapply(function(x, y) {c(CmaxI(x,y))}, matI_u, lntDU, SIMPLIFY=T)
lminI_u <- mapply(function(x, y) {c(CminI(x,y))}, matI_u, lntDU, SIMPLIFY=T)
lcovI_u <- c(lmaxI_u, lminI_u) # all coverage values for incidence data
lcovI_u <- data.frame("CoverageVals"=lcovI_u, 
                      "Coverage"=c("Cmax", "Cmax", "Cmax", "Cmin", "Cmin", "Cmin"),
                      "Group"=rep("Undisturbed", times=6),
                      "Type"=rep("PA", times=6),
                      "Taxa"=c("Fish", "Invertebrates", "Diatoms", "Fish", "Invertebrates", "Diatoms"))
range(lcovI_u$CoverageVals)  # 0.8515 0.9870

as.data.frame(rbind(lcovI, lcovI_d, lcovI_u))


# Agg data: --------------------------------------------------------------------
CmaxA <- function(x) {
  Abun <- x
  cmaxA <- apply(Abun, 2, function(x) iNEXT.3D:::Coverage(x, 'abundance', 2*sum(x))) %>% min %>% round(., 4)
  return(cmaxA)
}

CminA <- function(x) {
  Abun <- x
  cminA <- apply(Abun, 2, function(x) iNEXT.3D:::Coverage(x, 'abundance', sum(x))) %>% min %>% round(., 4)
  return(cminA)
}
lAllAbuAgg <- lapply(lAllAbuAgg, function(x) {do.call(cbind, x)})
lAbuAgg_d <- lapply(lAbuAgg_d, function(x) {do.call(cbind, x)})
lAbuAgg_u <- lapply(lAbuAgg_u, function(x) {do.call(cbind, x)})

cmax_aggAll <- lapply(lAllAbuAgg, CmaxA) # OK
cmax_aggD <- lapply(lAbuAgg_d, CmaxA)
cmax_aggU <- lapply(lAbuAgg_u, CmaxA)

cmin_aggAll <- lapply(lAllAbuAgg, CminA) # OK
cmin_aggD <- lapply(lAbuAgg_d, CminA)
cmin_aggU <- lapply(lAbuAgg_u, CminA)

agg_coverages <- unlist(c(cmax_aggAll, cmax_aggD, cmax_aggU,
                          cmin_aggAll, cmin_aggD, cmin_aggU))

lcovAgg <- data.frame("CoverageVals"=agg_coverages, 
                    "Coverage"=c("Cmax", "Cmax", "Cmax","Cmax", "Cmax", "Cmax","Cmax", "Cmax", "Cmax",
                                 "Cmin", "Cmin", "Cmin","Cmin", "Cmin", "Cmin","Cmin", "Cmin", "Cmin"),
                    "Group"=c("Northern Range", "Northern Range", "Northern Range",
                              "Disturbed", "Disturbed", "Disturbed",
                              "Undisturbed", "Undisturbed", "Undisturbed",
                              "Northern Range", "Northern Range", "Northern Range",
                              "Disturbed", "Disturbed", "Disturbed",
                              "Undisturbed", "Undisturbed", "Undisturbed"),
                    "Type"=rep("Agg", times=18),
                    "Taxa"=rep(c("Fish", "Invertebrates", "Diatoms"), times=6))
lcovAgg
lcov <- as.data.frame(rbind(lcovI, lcovI_d, lcovI_u ,lcovAgg))

# Save: ------------------------------------------------------------------------
save(lcov, file = paste0(mainanalysisRData, "/lcov.RData"))


# End of script ################################################################
################################################################################