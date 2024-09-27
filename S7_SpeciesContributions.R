################################################################################
# S5 Species contributions
# AFE
# October 2023 & February 2014
################################################################################

# LIBRARIES ====================================================================
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(data.table)
library(stringr)
library(stringi)
library(pheatmap)
library(patchwork)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggplotify)


rm(list=ls())


################################################################################
# FUNCTION CHAO-RICOTTA ========================================================

#' dis1 computes the contribution of each species/node to the two types of dissimilarity measures.
#' @param x is the species-by-assemblages abundance matrix with species names as rownames.
#' @param q is value for the diversity order.
#' @param type is tax (taxonomic) or phy (phylogenetic).
#' @param type2 is "species" or "k"."species" means the contribution of each species/node to the two types of dissimilarity measures 
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity).
#' "k" means the contribution of each assemblage/location/site to the two types of dissimilarity measures 
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity). In the worked example, the contribution of each assemblage/stage is not computed.
#' @tree is the pylog object of the phylogenetic tree of all assemblages.
#' @return the contribution of each species/node to the two types of dissimilarity measures: Jaccard-type (1-U_qN) and Sorensen-type (1-C_qN)
dis1 <- function(x, q, type = "tax", type2 = "species", tree = NULL){
  if(type2 == "species"){
    FUN <- rowSums
  }else{
    FUN <- colSums
  }
  if(type == "tax"){
    x <- as.matrix(x)
    x <- x[rowSums(x)>0, ]
    N <- ncol(x)
    zbar <- rowSums(x)/N
    x1 <- x[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    if(q==0){
      UqN <- FUN(x==0)/((N-1)*(sum(rowSums(x)>0)))
      CqN <- FUN(x==0)/((N-1)*(sum(apply(x, 2, function(i){sum(i>0)}))))
    }else if(q==2){
      UqN <- FUN((x1-zbar1)^2)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1-zbar1)^2)/((1-N^(1-q))*sum(x1^q))
    }else if(q!=1){
      UqN <- FUN((x1)^q-(zbar1)^q)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1)^q-(zbar1)^q)/((1-N^(1-q))*sum(x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(x1*log(x2), na.rm = T)/((sum(x)*log(N)))
      CqN <- UqN
    }
  }else{
    Li <- c(tree$leaves, tree$nodes)
    cumtree = function(a, tree){
      a <- a[names(tree$leaves)]
      for(i in 1:length(tree$parts)){
        a[1+length(a)] <- sum(a[tree$parts[[i]]])
        names(a)[length(a)] <- names(tree$parts)[i]
      }
      a
    }
    ai <- apply(x, 2, cumtree, tree)
    wt <- apply(ai, 1, function(x1)(sum(x1))^q/sum(Li*rowSums(ai, na.rm = T)^q))
    N <- ncol(ai)
    zbar <- rowSums(ai)/N
    x1 <- ai[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    Li <- Li[zbar>0]
    T1 <- sum(rowSums(x1)*Li)
    if(q==0){
      if(type2 == "species"){
        rn <- nrow(x1)
        UqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))})/((N-1)*sum(Li)) 
        CqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))/((N-1)*sum(Li*rowSums(x1!=0)))})
      }else{
        UqN <- apply(x1, 2, function(x){sum(Li[x==0])})/((N-1)*sum(Li)) 
        CqN <- apply(x1, 2, function(x){sum(Li[x==0])/((N-1)*sum(Li*colSums(x1!=0)))})
      }
      
    }else if(q==2){
      UqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else if(q!=1){
      UqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(Li*x1*log(x2), na.rm = T)/(T1*log(N))
      CqN <- UqN
    }
  }
  
  # c(sum(UqN), sum(CqN))
  rbind(UqN, CqN)
}


################################################################################
# READ & FORMAT DATA ===========================================================
mainanalysisRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/mainanalysisRData"
all_files <- list.files(mainanalysisRData, pattern = "\\.RData$", full.names = TRUE)
for (file in all_files) {
  load(file)
} # load files

process <- function(dt, nsamples=NULL){
  dt <- lapply(dt, function(x) {do.call(rbind, x)})
  dt <- lapply(dt, function(x){x$Session <- str_split_fixed(rownames(x), "\\.",2)[,1];x
  x$Taxa <- str_split_fixed(rownames(x), "\\.",2)[,2];x})
  
  dt <- lapply(dt, function(x) {x$Sums <- rowSums(x[,c(1:nsamples)]);x})
  dt <- lapply(dt, function(x) {x[,names(x) %in% c("Session", "Taxa", "Sums")]})
  dt <- lapply(dt, function(x) {x$Session <- as.numeric(x$Session);x})
  dt <- lapply(dt, function(x) {data.frame(spread(x, key="Session", value="Sums", fill=0))})
  dt <- lapply(dt, function(x) {rownames(x) <- x$Taxa;x
  x <-x[,-1];x
  names(x) <- gsub("X", "", names(x));x})
}

lAllPA <- process(lAllPA, nsamples = 16)
lPA_d <- process(lPA_d, nsamples = 8)
lPA_u <- process(lPA_u, nsamples = 8)

lAllPA$fishsen <- lAllPA$fish
lAllPA$fishsen <- lAllPA$fishsen[,c(1:19)] # sensitivity analysis (Fish period 2011-2015)


names(lAllPA$fish)[names(lAllPA$fish)==21] <- 46
names(lAllPA$fish)[names(lAllPA$fish)==22] <- 50
names(lAllPA$fish)[names(lAllPA$fish)==23] <- 54
names(lAllPA$fish)[names(lAllPA$fish)==20] <- 23

names(lPA_d$fish)[names(lPA_d$fish)==21] <- 46
names(lPA_d$fish)[names(lPA_d$fish)==22] <- 50
names(lPA_d$fish)[names(lPA_d$fish)==23] <- 54
names(lPA_d$fish)[names(lPA_d$fish)==20] <- 23

names(lPA_u$fish)[names(lPA_u$fish)==21] <- 46
names(lPA_u$fish)[names(lPA_u$fish)==22] <- 50
names(lPA_u$fish)[names(lPA_u$fish)==23] <- 54
names(lPA_u$fish)[names(lPA_u$fish)==20] <- 23                     # appropriate format for function dis

extract_info <- function(input_string) {
  words <- strsplit(input_string, "_")[[1]]
  return(paste0(substr(words[1], 1, 1), "_", words[2]))
}
extract_third_word <- function(input_string) {
  words <- strsplit(input_string, "_")[[1]]
  
  if (length(words) == 3) {
    return(words[3])
  } else {
    return(input_string)
  }
}

rownames(lAllPA$fish) <- sapply(rownames(lAllPA$fish), extract_info)
rownames(lPA_d$fish) <- sapply(rownames(lPA_d$fish), extract_info)
rownames(lPA_u$fish) <- sapply(rownames(lPA_u$fish), extract_info)                # name abbv in heatmaps

rownames(lAllPA$bi) <- sapply(rownames(lAllPA$bi), extract_third_word)
rownames(lAllPA$bi) <- str_to_title(rownames(lAllPA$bi))
rownames(lAllPA$bi)[rownames(lAllPA$bi)=="Unknown"] <- "Lepidoptera"              # corrected after extract_third_word
rownames(lAllPA$bi)[rownames(lAllPA$bi)=="Entognatha_collembola"] <- "Collembola" # corrected after extract_third_word

rownames(lPA_d$bi) <- sapply(rownames(lPA_d$bi), extract_third_word)
rownames(lPA_d$bi) <- str_to_title(rownames(lPA_d$bi))
rownames(lPA_d$bi)[rownames(lPA_d$bi)=="Unknown"] <- "Lepidoptera"
rownames(lPA_d$bi)[rownames(lPA_d$bi)=="Entognatha_collembola"] <- "Collembola"

rownames(lPA_u$bi) <- sapply(rownames(lPA_u$bi), extract_third_word)
rownames(lPA_u$bi) <- str_to_title(rownames(lPA_u$bi))
rownames(lPA_u$bi)[rownames(lPA_u$bi)=="Unknown"] <- "Lepidoptera"
rownames(lPA_u$bi)[rownames(lPA_u$bi)=="Entognatha_collembola"] <- "Collembola"


rownames(lAllPA$fishsen) <- sapply(rownames(lAllPA$fishsen), extract_info)
    


lAllPA_FO <- lAllPA
lPA_d_FO <- lPA_d
lPA_u_FO <- lPA_u


# Format for pair-wise dissimilarities: ----------------------------------------
combosdf <- function(x, y){
  list_of_dfs <- unlist(lapply(x, function(a) lapply(y, function (b) as.data.frame(do.call(cbind,list(a, b))))), recursive=FALSE) # needed structure to compute lag beta diversity
  return(list_of_dfs)
}
structureforBeta <- function(dt){
  rows <- lapply(dt, function(x) {rownames(x)})
  l2 <- mapply(function(x,y){combosdf(x,y)}, dt, dt, SIMPLIFY = F) 
  l2 <- mapply(function(x,y){lapply(x, function(z) {rownames(z) <-y; z})}, l2, rows, SIMPLIFY = F) 
  l2 <- lapply(l2, function(x) {lapply(x, function(y) {as.data.frame(y)})})
  l2 <- lapply(l2, function(x) {Map(cbind, x, "Combo" = names(x))})
  l2 <- lapply(l2, function(x) {lapply(x, function(y) {y$s1 <- as.numeric(str_split_fixed(y$Combo, "\\.", 2)[,1]);y
  y$s2 <- as.numeric(str_split_fixed(y$Combo, "\\.", 2)[,2]);y})})
  l3 <- lapply(l2, function(x) {as.data.frame(do.call(rbind, x))})
  l3 <- lapply(l3, function(x) {x$ConcatCombo <- apply(x[c("s1", "s2")], 1, function(z) paste(sort(z), collapse = "_"));x})
  l3 <- lapply(l3, function(x) {x$Taxa <- str_split_fixed(rownames(x),"\\.",3 )[,3];x})
  l4 <- lapply(l3, function(x) {x <- x %>%
    group_by(ConcatCombo, Taxa) %>%
    slice(1)})
  l4 <- lapply(l4, function(x) {x[!(x$s1==x$s2),]})
  l5 <- lapply(l4, function(x) {split(x, f=x$ConcatCombo)})
  l5 <- lapply(l5, function(x) {lapply(x, function(y){as.data.frame(y)})})
  l5 <- lapply(l5, function(x) {lapply(x, function(y){rownames(y)=y$Taxa;y})})
  l5 <- lapply(l5, function(x) {lapply(x, function(y){within(y, rm(Combo, s1, s2, ConcatCombo, Taxa))})})
  return(l5)
}

lAllPA <- structureforBeta(lAllPA)
lPA_d <- structureforBeta(lPA_d)
lPA_u <- structureforBeta(lPA_u)


################################################################################
# FREQ OCCURRENCES =============================================================
class(lAllPA_FO$fish)

lAllPA_FO <- lapply(lAllPA_FO, function(x) {x$Taxa <- rownames(x);x})
lAllPA_FO <- lapply(lAllPA_FO, function(x) {gather(x, key="Session", value="Freq", -ncol(x))})
lAllPA_FO <- lapply(lAllPA_FO, function(x) {x$Session <- as.integer(x$Session);x
x$Freq <- as.integer(x$Freq);x})

lPA_d_FO <- lapply(lPA_d_FO, function(x) {x$Taxa <- rownames(x);x})
lPA_d_FO <- lapply(lPA_d_FO, function(x) {gather(x, key="Session", value="Freq", -ncol(x))})
lPA_d_FO <- lapply(lPA_d_FO, function(x) {x$Session <- as.integer(x$Session);x
x$Freq <- as.integer(x$Freq);x})

lPA_u_FO <- lapply(lPA_u_FO, function(x) {x$Taxa <- rownames(x);x})
lPA_u_FO <- lapply(lPA_u_FO, function(x) {gather(x, key="Session", value="Freq", -ncol(x))})
lPA_u_FO <- lapply(lPA_u_FO, function(x) {x$Session <- as.integer(x$Session);x
x$Freq <- as.integer(x$Freq);x})

c1check <- lapply(lAllPA_FO, function(x) {x %>%
  group_by(Taxa) %>%
  summarise(Count = n())}) # OK, a number for each session even if it is 0

Lms <- function(x){
  store <- list()
  for (i in unique(x$Taxa)){
    subsetTaxa <- subset(x, x$Taxa==i)
    lmFreq <- glm(Freq ~ Session, family = poisson(link = "log"), data = subsetTaxa)
    slope <- coef(lmFreq)["Session"]
    pval <- summary(lmFreq)$coefficients["Session", "Pr(>|z|)"]
    taxa <- i
    store[[i]] <- t(as.data.frame(rbind(taxa, slope, pval)))
  }
  return(store)
}

lAllPAFreq <- lapply(lAllPA_FO, function(x) {Lms(x)})
lAllPAFreq <- lapply(lAllPAFreq, function(x) {as.data.frame(do.call(rbind, x))})
lAllPAFreq <- lapply(lAllPAFreq, function(x) {x$color <- ifelse(x$pval <=0.05 & x$slope <0,"red", "#FF9966");x
x$color <- ifelse(x$pval <=0.05 & x$slope >0,"blue", x$color);x
x$color <- ifelse(x$pval >0.05 & x$slope >0,"#56B4E9", x$color);x}) 

lPA_d_Freq <- lapply(lPA_d_FO, function(x) {Lms(x)})
lPA_d_Freq <- lapply(lPA_d_Freq, function(x) {as.data.frame(do.call(rbind, x))})
lPA_d_Freq <- lapply(lPA_d_Freq, function(x) {x$color <- ifelse(x$pval <=0.05 & x$slope <0,"red", "#FF9966");x
x$color <- ifelse(x$pval <=0.05 & x$slope >0,"blue", x$color);x
x$color <- ifelse(x$pval >0.05 & x$slope >0,"#56B4E9", x$color);x}) 

lPA_u_Freq <- lapply(lPA_u_FO, function(x) {Lms(x)})
lPA_u_Freq <- lapply(lPA_u_Freq, function(x) {as.data.frame(do.call(rbind, x))})
lPA_u_Freq <- lapply(lPA_u_Freq, function(x) {x$color <- ifelse(x$pval <=0.05 & x$slope <0,"red", "#FF9966");x
x$color <- ifelse(x$pval <=0.05 & x$slope >0,"blue", x$color);x
x$color <- ifelse(x$pval >0.05 & x$slope >0,"#56B4E9", x$color);x}) 


################################################################################
# CONTRIBUTIONS ================================================================
contributions <- function(dt, orderq=NULL, typeSL=NULL){
  cont <- lapply(dt, function(x) {lapply(x, function(y) {t(dis1(y, q=orderq, type="tax", type2 = typeSL))})})
  cont <- lapply(cont, function(x) {Map(cbind, x, "Combo" = names(x))})
  cont <- lapply(cont, function(x) {as.data.frame(do.call(rbind,x))})
  cont <- lapply(cont, function(x) {x$Order.q <- paste0("q", orderq);x
  x$s1 <- str_split_fixed(x$Combo, "_", 2)[,1];x
  x$s2 <- str_split_fixed(x$Combo, "_", 2)[,2];x
  x$Lag <- abs(as.numeric(x$s1)-as.numeric(x$s2));x})
}


################################################################################
# FORMAT HEATMAPS ==============================================================
PlotsPath <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
diss <- function(x, orderq=NULL){
SpsContPA <- contributions(x, orderq = orderq, typeSL="species")
SpsContPA <- lapply(SpsContPA, function(x) {x$Taxa <- str_split_fixed(rownames(x), "\\.", 2)[,1];x})
SpsContPA <- lapply(SpsContPA, function(x) {x$UqN <- as.numeric(x$UqN);x}) # select Jaccard
SpsContPA <- lapply(SpsContPA, function(x) {x %>% group_by(Taxa, Lag) %>% summarise(Cmean = mean(UqN))}) # mean cont per taxa & lag
SpsContPA <- lapply(SpsContPA, function(x) {x[order(as.numeric(x$Lag), decreasing = F),]})
SpsContPA <- lapply(SpsContPA, function(x) {spread(x, key="Lag", value="Cmean")})
SpsContPA <- lapply(SpsContPA, function(x) {as.data.frame(x)})
SpsContPA <- lapply(SpsContPA, function(x) {rownames(x) <- x$Taxa; x})
SpsContPA <- lapply(SpsContPA, function(x) {within(x, rm(Taxa))})
SpsContPA <- lapply(SpsContPA, function(x) {x[,order(as.numeric(names(x)))]})
SpsContPA <- lapply(SpsContPA, function(x) {x$Rows <- rowSums(x, na.rm=T);x}) # positive only if at least 1 cont.
SpsContPA <- lapply(SpsContPA, function(x) {subset(x, x$Rows > 0)})           # rm sps with no contribution
SpsContPA <- lapply(SpsContPA, function(x) {within(x, rm(Rows))})
return(SpsContPA)
}

SpsCont0PA <- diss(lAllPA, orderq = 0)
SpsCont2PA <- diss(lAllPA, orderq = 2)
SpsCont0d <- diss(lPA_d, orderq = 0)
SpsCont2d <- diss(lPA_d, orderq = 2)
SpsCont0u <- diss(lPA_u, orderq = 0)
SpsCont2u <- diss(lPA_u, orderq = 2)


colfuncF <- colorRampPalette(c("white", "#79b9e7"))
colfuncF2 <- colorRampPalette(c("white", "#2b4b9c"))

colfuncI <- colorRampPalette(c("white", "#e7bd26"))
colfuncI2 <- colorRampPalette(c("white", "#ad8e1d"))

colfuncD <- colorRampPalette(c("white", "#96b24f"))
colfuncD2 <- colorRampPalette(c("white", "#2b5108"))


# C5 NR: -----------------------------------------------------------------------

# Fish: 
Cont0F <- SpsCont0PA[["fish"]]
Cont0F$color <- lAllPAFreq$fish$color[match(rownames(Cont0F), lAllPAFreq$fish$taxa)]
rownames(Cont0F)[rownames(Cont0F)=="S_marmoratus"] <- "S_marmoratus    "       # for a slightly better alignment
Cont2F <- SpsCont2PA[["fish"]]
Cont2F$color <- lAllPAFreq$fish$color[match(rownames(Cont2F), lAllPAFreq$fish$taxa)]
rownames(Cont2F)[rownames(Cont2F)=="A_bimaculatus"] <- "A_bimaculatus   "    # for a slightly better alignment

# Invertebrates:
Cont0I <- SpsCont0PA[["bi"]]
Cont0I$color <- lAllPAFreq$bi$color[match(rownames(Cont0I), lAllPAFreq$bi$taxa)]
rownames(Cont0I)[rownames(Cont0I)=="Crustacea_spp"] <- "Crustacea_spp     "  # for a slightly better alignment
Cont2I <- SpsCont2PA[["bi"]]
Cont2I$color <- lAllPAFreq$bi$color[match(rownames(Cont2I), lAllPAFreq$bi$taxa)]

# Diatoms:
Cont0D <- SpsCont0PA[["dia"]]
Cont0D$color <- lAllPAFreq$dia$color[match(rownames(Cont0D), lAllPAFreq$dia$taxa)]
Cont2D <- SpsCont2PA[["dia"]]
Cont2D$color <- lAllPAFreq$dia$color[match(rownames(Cont2D), lAllPAFreq$dia$taxa)]

# Fish (SEN): 
Cont0FSEN <- SpsCont0PA[["fishsen"]]
Cont0FSEN$color <- lAllPAFreq$fishsen$color[match(rownames(Cont0FSEN), lAllPAFreq$fishsen$taxa)]
rownames(Cont0FSEN)[rownames(Cont0FSEN)=="S_marmoratus"] <- "S_marmoratus    "       # for a slightly better alignment
Cont2FSEN <- SpsCont2PA[["fishsen"]]
Cont2FSEN$color <- lAllPAFreq$fishsen$color[match(rownames(Cont2FSEN), lAllPAFreq$fishsen$taxa)]
rownames(Cont2FSEN)[rownames(Cont2FSEN)=="A_bimaculatus"] <- "A_bimaculatus   "    # for a slightly better alignment


# C6 Disturbed: ----------------------------------------------------------------

# Fish:
Cont0Fd <- SpsCont0d[["fish"]]
Cont0Fd$color <- lPA_d_Freq$fish$color[match(rownames(Cont0Fd), lPA_d_Freq$fish$taxa)]
rownames(Cont0Fd)[rownames(Cont0Fd)=="S_marmoratus"] <- "S_marmoratus    "        # for a slightly better alignment
Cont2Fd <- SpsCont2d[["fish"]]
Cont2Fd$color <- lPA_d_Freq$fish$color[match(rownames(Cont2Fd), lPA_d_Freq$fish$taxa)]
rownames(Cont2Fd)[rownames(Cont2Fd)=="A_bimaculatus"] <- "A_bimaculatus   "     # for a slightly better alignment

# Invertebrates:
Cont0Id <- SpsCont0d[["bi"]]
Cont0Id$color <- lPA_d_Freq$bi$color[match(rownames(Cont0Id), lPA_d_Freq$bi$taxa)]
rownames(Cont0Id)[rownames(Cont0Id)=="Crustacea_spp"] <- "Crustacea_spp        "# for a slightly better alignment
Cont2Id <- SpsCont2d[["bi"]]
Cont2Id$color <- lPA_d_Freq$bi$color[match(rownames(Cont2Id), lPA_d_Freq$bi$taxa)]

# Diatoms:
Cont0Dd <- SpsCont0d[["dia"]]
Cont0Dd$color <- lPA_d_Freq$dia$color[match(rownames(Cont0Dd), lPA_d_Freq$dia$taxa)]
Cont2Dd <- SpsCont2d[["dia"]]
Cont2Dd$color <- lPA_d_Freq$dia$color[match(rownames(Cont2Dd), lPA_d_Freq$dia$taxa)]


# C6 Undisturbed ---------------------------------------------------------------

# Fish:
Cont0Fu <- SpsCont0u[["fish"]]
Cont0Fu$color <- lPA_u_Freq$fish$color[match(rownames(Cont0Fu), lPA_u_Freq$fish$taxa)]
rownames(Cont0Fu)[rownames(Cont0Fu)=="S_marmoratus"] <- "S_marmoratus    "      # for a slightly better alignment
Cont2Fu <- SpsCont2u[["fish"]]
Cont2Fu$color <- lPA_u_Freq$fish$color[match(rownames(Cont2Fu), lPA_u_Freq$fish$taxa)]
rownames(Cont2Fu)[rownames(Cont2Fu)=="A_bimaculatus"] <- "A_bimaculatus   "     # for a slightly better alignment

# Invertebrates:
Cont0Iu <- SpsCont0u[["bi"]]
Cont0Iu$color <- lPA_u_Freq$bi$color[match(rownames(Cont0Iu), lPA_u_Freq$bi$taxa)]
rownames(Cont0Iu)[rownames(Cont0Iu)=="Crustacea_spp"] <- "Crustacea_spp        "# for a slightly better alignment
Cont2Iu <- SpsCont2u[["bi"]]
Cont2Iu$color <- lPA_u_Freq$bi$color[match(rownames(Cont2Iu), lPA_u_Freq$bi$taxa)]

# Diatoms:
Cont0Du <- SpsCont0u[["dia"]]
Cont0Du$color <- lPA_u_Freq$dia$color[match(rownames(Cont0Du), lPA_u_Freq$dia$taxa)]
Cont2Du <- SpsCont2u[["dia"]]
Cont2Du$color <- lPA_u_Freq$dia$color[match(rownames(Cont2Du), lPA_u_Freq$dia$taxa)]


################################################################################
# HEATMAPS: ====================================================================
# In q=0 diss, common sps may not have an effect, bc only richness counts
# yet all taxa have a contribution to q=2 diss, even if for rare sps this 
# cont is almost negligible.


# Chapter 5 NR: ================================================================
# 53-3 cols that don't exist
dim(Cont0F) # 50 + 1 for color
names(Cont0F)

pNR_F_0 <- pheatmap(Cont0F[,c(1:50)], cluster_cols=FALSE, color = colfuncF(16), cluster_rows = F, show_colnames = F)
cols <- Cont0F[order(match(rownames(Cont0F), pNR_F_0$gtable$grobs[[2]]$label)), ]$color
pNR_F_0$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_F_0
pNR_F_0 <- as.ggplot(pNR_F_0) # fish order q0

pNR_F_2 <- pheatmap(Cont2F[,c(1:50)], cluster_cols=FALSE, color = colfuncF2(16), cluster_rows = F)
cols <- Cont2F[order(match(rownames(Cont2F), pNR_F_2$gtable$grobs[[3]]$label)), ]$color
pNR_F_2$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_F_2
pNR_F_2 <- as.ggplot(pNR_F_2) # fish order q2

pNR_I_0 <- pheatmap(Cont0I[,c(1:18)], cluster_cols=FALSE, color = colfuncI(16), cluster_rows = F, show_colnames = F)
cols <- Cont0I[order(match(rownames(Cont0I), pNR_F_0$gtable$grobs[[2]]$label)), ]$color
pNR_I_0$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_I_0
pNR_I_0 <- as.ggplot(pNR_I_0) # inv order q0

pNR_I_2 <- pheatmap(Cont2I[,c(1:18)], cluster_cols=FALSE, color = colfuncI2(16), cluster_rows = F)
cols <- Cont2I[order(match(rownames(Cont2I), pNR_I_2$gtable$grobs[[3]]$label)), ]$color
pNR_I_2$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_I_2
pNR_I_2 <- as.ggplot(pNR_I_2) # inv order q0

pNR_D_0 <- pheatmap(Cont0D[,c(1:18)], cluster_cols=FALSE, color = colfuncD(16), cluster_rows = F, show_colnames = F)
cols <- Cont0D[order(match(rownames(Cont0D), pNR_D_0$gtable$grobs[[2]]$label)), ]$color
pNR_D_0$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_D_0
pNR_D_0 <- as.ggplot(pNR_D_0) # dia order q0

pNR_D_2 <- pheatmap(Cont2D[,c(1:18)], cluster_cols=FALSE, color = colfuncD2(16), cluster_rows = F) 
cols <- Cont2D[order(match(rownames(Cont2D), pNR_D_2$gtable$grobs[[3]]$label)), ]$color
pNR_D_2$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_D_2
pNR_D_2 <- as.ggplot(pNR_D_2) # dia order q2

pNR_FSEN_0 <- pheatmap(Cont0FSEN[,c(1:18)], cluster_cols=FALSE, color = colfuncF(16), cluster_rows = F, show_colnames = F)
cols <- Cont0FSEN[order(match(rownames(Cont0FSEN), pNR_FSEN_0$gtable$grobs[[2]]$label)), ]$color
pNR_FSEN_0$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_FSEN_0
pNR_FSEN_0 <- as.ggplot(pNR_FSEN_0) # fish (sen) order q0

pNR_FSEN_2 <- pheatmap(Cont2FSEN[,c(1:18)], cluster_cols=FALSE, color = colfuncF2(16), cluster_rows = F)
cols <- Cont2FSEN[order(match(rownames(Cont2FSEN), pNR_FSEN_2$gtable$grobs[[2]]$label)), ]$color
pNR_FSEN_2$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_FSEN_2
pNR_FSEN_2 <- as.ggplot(pNR_FSEN_2) # fish (sen) order q2


fish_cont <- pNR_F_0 / pNR_F_2 + plot_layout(heights = c(1, 2))   
bi_cont <- pNR_I_0 / pNR_I_2 + plot_layout(heights = c(1, 2))
dia_cont <- pNR_D_0 / pNR_D_2 + plot_layout(heights = c(1, 2)) # other edits on Power Point
fish_cont_sen <- pNR_FSEN_0 / pNR_FSEN_2 + plot_layout(heights = c(1, 2))   

ggsave(paste0(PlotsPath, "/Chapter5/S7_C5_fish_cont.png"), fish_cont, width = 10, height = 8)
ggsave(paste0(PlotsPath, "/Chapter5/S7_C5_bi_cont.png"), bi_cont, width = 10, height = 8)
ggsave(paste0(PlotsPath, "/Chapter5/S7_C5_dia_cont.png"), dia_cont, width = 10, height = 8)
ggsave(paste0(PlotsPath, "/Chapter5/S7_C5_fishsen_cont.png"), fish_cont_sen, width = 10, height = 8)


# Chapter 6 Disturbed: ---------------------------------------------------------
pNR_F_0d <- pheatmap(Cont0Fd[,c(1:50)], cluster_cols=FALSE, color = colfuncF(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Fd[order(match(rownames(Cont0Fd), pNR_F_0d$gtable$grobs[[2]]$label)), ]$color
pNR_F_0d$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_F_0d
pNR_F_0d <- as.ggplot(pNR_F_0d)

pNR_F_2d <- pheatmap(Cont2Fd[,c(1:50)], cluster_cols=FALSE, color = colfuncF2(16), cluster_rows = F)
cols <- Cont2Fd[order(match(rownames(Cont2Fd), pNR_F_2d$gtable$grobs[[3]]$label)), ]$color
pNR_F_2d$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_F_2d
pNR_F_2d <- as.ggplot(pNR_F_2d)

pNR_I_0d <- pheatmap(Cont0Id[,c(1:18)], cluster_cols=FALSE, color = colfuncI(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Id[order(match(rownames(Cont0Id), pNR_I_0d$gtable$grobs[[2]]$label)), ]$color
pNR_I_0d$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_I_0d
pNR_I_0d <- as.ggplot(pNR_I_0d)

pNR_I_2d <- pheatmap(Cont2Id[,c(1:18)], cluster_cols=FALSE, color = colfuncI2(16), cluster_rows = F)
cols <- Cont2Id[order(match(rownames(Cont2Id), pNR_I_2d$gtable$grobs[[3]]$label)), ]$color
pNR_I_2d$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_I_2d
pNR_I_2d <- as.ggplot(pNR_I_2d)

pNR_D_0d <- pheatmap(Cont0Dd[,c(1:18)], cluster_cols=FALSE, color = colfuncD(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Dd[order(match(rownames(Cont0Dd), pNR_D_0d$gtable$grobs[[2]]$label)), ]$color
pNR_D_0d$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_D_0d
pNR_D_0d <- as.ggplot(pNR_D_0d)

pNR_D_2d <- pheatmap(Cont2Dd[,c(1:18)], cluster_cols=FALSE, color = colfuncD2(16), cluster_rows = F)
cols <- Cont2Dd[order(match(rownames(Cont2Dd), pNR_D_2d$gtable$grobs[[3]]$label)), ]$color
pNR_D_2d$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_D_2d
pNR_D_2d <- as.ggplot(pNR_D_2d)


ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q0_fish_cont.png"), pNR_F_0d, width = 10, height = 4) # disturbed order q0
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q0_bi_cont.png"), pNR_I_0d, width = 10, height = 4)
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q0_dia_cont.png"), pNR_D_0d, width = 10, height = 4)

ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q2_fish_cont.png"), pNR_F_2d, width = 10, height = 8) # disturbed order q2
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q2_bi_cont.png"), pNR_I_2d, width = 10, height = 8)
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_DIS_q2_dia_cont.png"), pNR_D_2d, width = 10, height = 8)


# Chapter 6 Undisturbed: -------------------------------------------------------
pNR_F_0u <- pheatmap(Cont0Fu[,c(1:50)], cluster_cols=FALSE, color = colfuncF(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Fu[order(match(rownames(Cont0Fu), pNR_F_0u$gtable$grobs[[2]]$label)), ]$color
pNR_F_0u$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_F_0u
pNR_F_0u <- as.ggplot(pNR_F_0u)

pNR_F_2u <- pheatmap(Cont2Fu[,c(1:50)], cluster_cols=FALSE, color = colfuncF2(16), cluster_rows = F)
cols <- Cont2Fu[order(match(rownames(Cont2Fu), pNR_F_2u$gtable$grobs[[3]]$label)), ]$color
pNR_F_2u$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_F_2u
pNR_F_2u <- as.ggplot(pNR_F_2u)

pNR_I_0u <- pheatmap(Cont0Iu[,c(1:18)], cluster_cols=FALSE, color = colfuncI(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Iu[order(match(rownames(Cont0Iu), pNR_I_0u$gtable$grobs[[2]]$label)), ]$color
pNR_I_0u$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_I_0u
pNR_I_0u <- as.ggplot(pNR_I_0u)

pNR_I_2u <- pheatmap(Cont2Iu[,c(1:18)], cluster_cols=FALSE, color = colfuncI2(16), cluster_rows = F)
cols <- Cont2Iu[order(match(rownames(Cont2Iu), pNR_I_2u$gtable$grobs[[3]]$label)), ]$color
pNR_I_2u$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_I_2u
pNR_I_2u <- as.ggplot(pNR_I_2u)

pNR_D_0u <- pheatmap(Cont0Du[,c(1:18)], cluster_cols=FALSE, color = colfuncD(16), cluster_rows = F, show_colnames = F)
cols <- Cont0Du[order(match(rownames(Cont0Du), pNR_D_0u$gtable$grobs[[2]]$label)), ]$color
pNR_D_0u$gtable$grobs[[2]]$gp <- gpar(col=cols)
pNR_D_0u
pNR_D_0u <- as.ggplot(pNR_D_0u)

pNR_D_2u <- pheatmap(Cont2Du[,c(1:18)], cluster_cols=FALSE, color = colfuncD2(16), cluster_rows = F)
cols <- Cont2Du[order(match(rownames(Cont2Du), pNR_D_2u$gtable$grobs[[3]]$label)), ]$color
pNR_D_2u$gtable$grobs[[3]]$gp <- gpar(col=cols)
pNR_D_2u
pNR_D_2u <- as.ggplot(pNR_D_2u)


ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q0_fish_cont.png"), pNR_F_0u, width = 10, height = 4) # undisturbed order q0
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q0_bi_cont.png"), pNR_I_0u, width = 10, height = 4)
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q0_dia_cont.png"), pNR_D_0u, width = 10, height = 4)

ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q2_fish_cont.png"), pNR_F_2u, width = 10, height = 8) # undisturbed order q2
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q2_bi_cont.png"), pNR_I_2u, width = 10, height = 8)
ggsave(paste0(PlotsPath, "/Chapter6/S7_C6_UNDIS_q2_dia_cont.png"), pNR_D_2u, width = 10, height = 8)


# End of script ################################################################
################################################################################