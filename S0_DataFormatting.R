################################################################################
# S0 Data Formatting (Initial for Amy's database downloads)
# AFE
# Jan 2024
################################################################################


# libraries: ===================================================================
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)

rm(list=ls())

raw_data_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/raw_datasets/RawDat_TrinidadNR"
formatted_data_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/raw_datasets/FormattedDataOutputs"

################################################################################
# read data: ===================================================================

# fish: ------------------------------------------------------------------------
fish <- read.csv2(paste0(raw_data_path,"/fishAFE2023_typos_corrected.csv"), h=T)
str(fish)
fishlengths <- read.csv2(paste0(raw_data_path,"/fish2023AFE_lengths.csv"), h=T) # sent a posteriori by Amy, typos corrected too
str(fishlengths)
identical(sort(fishlengths$weight), sort(fish$weight)) # TRUE (so only the order might be different)

# NOTE: continue with fish lengths file.

# inv: -------------------------------------------------------------------------
inv <- read.csv(paste0(raw_data_path,"/BioTIME_Trinidad_Invertebrates.csv"), h=T)

# dia: -------------------------------------------------------------------------
dia <- read.csv(paste0(raw_data_path,"/BioTIME_Trinidad_Diatoms.csv"), h=T)

################################################################################
# format data ==================================================================

# fish: ------------------------------------------------------------------------

################################################################################
# inv: -------------------------------------------------------------------------
str(inv)
sum(is.na(inv)) # 0

range(inv$timestep) # already named, Jan 2011 to August 2015 (yet the replicabilities chosen for 2011 are unknown)
range(inv$abundance) # 1 1616

# inv, taxonomy: ---------------------------------------------------------------
length(unique(inv$Classification)) # 47
sort(unique(inv$Classification))


# LARGER GROUPS:
# Hydracarina (unranked under order Trombidiformes)
# Crustacea (subphylum)
# Collembola (subclass)
# Coleoptera (order)
# Diptera  (order)
# Ephemeroptera (order)
# Hemiptera (order)
# Lepidoptera  (order)
# Odonata (order)
# Plecoptera (order)
# Trichoptera (order)
# Insecta (class)
# Oligochaeta
# Turbellaria


inv$Classification <- tolower(inv$Classification) # fix same names with differences in capital letters
inv$Classification[inv$Classification=="turbellaria turbellaria"] <-"turbellaria"

sum(inv$Classification=="crustacea decapoda unknown") # 2
sum(inv$Classification=="crustacea unknown")          # 18

inv$Classification[inv$Classification %in% c("crustacea decapoda unknown", "crustacea unknown")] <- "crustacea spp"

inv$Classification[inv$Classification=="oligochaeta unknown"] <- "oligochaeta"
inv$Classification[inv$Classification=="entognatha collembola unknown"] <- "entognatha collembola"
# "entognatha collembola unknown" to "entognatha collembola"
# "oligochaeta unknown" to ""oligochaeta"
sort(unique(inv$Classification))


sum(inv$Classification=="insecta coleoptera unknown") # 75 
sum(inv$Classification=="insecta diptera unknown")    # 18 
sum(inv$Classification=="insecta unknown")            # 15 
sum(inv$Classification=="insecta trichoptera unknown")# 3
3+75+18+15 # 111 occurrences
111/6104*100  # 1.8 %

# 1) Remove all insecta unknown:
unknowninsecta <- c("insecta coleoptera unknown", "insecta diptera unknown",
                    "insecta unknown", "insecta trichoptera unknown")
inv1 <- inv[!inv$Classification %in% unknowninsecta,]
length(unique(inv1$Classification))  # 30
sort(unique(inv1$Classification))
inv1 <- inv1 %>% group_by(timestep, Site, Classification) %>% summarise(abundance=sum(abundance))

# 2) keep as insecta spp.:
inv2 <- inv
inv2$Classification <- ifelse(inv2$Classification %in% unknowninsecta, "insecta spp", inv2$Classification)
length(unique(inv2$Classification))  # 31
sort(unique(inv2$Classification))
inv2 <- inv2 %>% group_by(timestep, Site, Classification) %>% summarise(abundance=sum(abundance))

BenthicInv1 <- inv1
BenthicInv2 <- inv2
save(BenthicInv1, file=paste0(formatted_data_path, "/BenthicInv1.RData"))  # no insecta unknwons
save(BenthicInv2, file=paste0(formatted_data_path, "/BenthicInv2.RData"))  # with insecta spp

################################################################################
# NOTES INVERTEBRATES:##########################################################
# 1) some unknown records pooled or removed after consulting with Amy.


################################################################################
# dia: -------------------------------------------------------------------------
str(dia)
sum(is.na(dia)) # 0

range(dia$timestep) # already named, Jan 2011 to August 2015 (yet the replicabilities chosen for 2011 are unknown)
range(dia$abundance) # 1 3408

# dia, taxonomy: ---------------------------------------------------------------
length(unique(dia$morphospecies)) # 26 
unique(dia$morphospecies)

#"A" (F)  "AE"(F) "AG" (F) "AH"(F) "AL"(F) "AM" (F) "AP"(F) "AQ" (F) "AV"(F) 
#"BB"(F) "BD" (F) "BJ" (F) "BK"(F) "BT"(F) "BV" (F)
#"C" (F)  "CB"(F) "CC" (F) "CI" (F)
#"D"(F)  "E" (F) "G"(F)  "H"(F)  "L"(F)  "N"(F)  U

sum(dia$morphospecies=="U")
unique(dia$site[dia$morphospecies=="U"])
unique(dia$timestep[dia$morphospecies=="U"])

names(dia) <- c("timestep", "site", "taxa", "abundance")

Diatoms <- dia
save(Diatoms, file=paste0(formatted_data_path, "/Diatoms.RData"))

################################################################################
# NOTES DIATOMS:################################################################
# 1) There are 2 lists: a more conservative and less conservative list.
# The more conservative list grouped together what were initially given different letters, 
# as they were fairly similar and we could not be certain that they were not just 
# variations of the same species. 

# 2) The PNAS dataset is the conservative list

# 3) Morphospecies U is the only one not in the key, however after consulting with
# Amy, this must be a key error and "U" is a legit category. She said is NOT likely "Unkown".

# 4) Morphospecies key in raw_data_path (conservative is ConListvCL)

# 5) "F" means found in key

# 6) PENDING: ask Amy for the downloadable database to see if I can disentangle 
# what U is (perhaps another letter in diatomConList)


# End of script ################################################################
################################################################################