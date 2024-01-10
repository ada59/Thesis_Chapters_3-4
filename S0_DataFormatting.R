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
head(fishlengths)

# step 1:-----------------------------------------------------------------------
sort(unique(fishlengths$Field))   # comments on weight estimations
sort(unique(fishlengths$Field.2)) # 0
sort(unique(fishlengths$Field.3)) # 0
sort(unique(fishlengths$Field.4)) # 0

fish <- fishlengths
fish <- within(fish, rm(Field, Field.2, Field.3, Field.4, Fish.Family, Genus)) # rm uneeded cols


# step 2: ----------------------------------------------------------------------
class(fish$Date)
fish$Date <- gsub("/", "_", fish$Date)
fish$Day <- str_split_fixed(fish$Date, "_", 3)[,1]
fish$Month <- str_split_fixed(fish$Date, "_", 3)[,2]
fish$Year <- str_split_fixed(fish$Date, "_", 3)[,3] # split info on date


# step 3: ----------------------------------------------------------------------
fish <- fish %>% relocate(c(Day, Month, Year), .before=Species)
fish <- fish %>% relocate(Disturbed, .after=Site)
fish <- fish %>% relocate(c(number.caught, number.seen), .before=weight)
fish <- fish %>% relocate(length, .after =weight)
fish <- fish %>% relocate(c(females, males, juveniles), .after =length) 
fish <- fish %>% relocate(HCvsEF, .before = weight)
fish <- within(fish, rm(Date)) # better reorganise data


# step 4: ----------------------------------------------------------------------
sort(unique(fish$Year))
sum(fish$Year=="") # 1
# View(fish[fish$Year=="",])  # Seems it is just the first obs, OK
fish <- fish[!fish$Year=="",] # make sure no blank records kept
fish %>% group_by(Year) %>% summarise(n_distinct(Site)) # OK


#step 5: -----------------------------------------------------------------------

# 2010:
fish$Session <- ifelse(fish$Year=="2010", 0, NA)  # will be removed prior to analyses

# 2011
fish$Session <- ifelse(fish$Year=="2011" & fish$Month %in% c("01", "02"), 1, fish$Session)
fish$Session <- ifelse(fish$Year=="2011" & fish$Month %in% c("05", "06"), 2, fish$Session) # Repeatability in Maracas D & U in June 2011
fish$Session <- ifelse(fish$Year=="2011" & fish$Month %in% c("08"), 3, fish$Session)       # Repeatability in Acono D & U in August 2011
fish$Session <- ifelse(fish$Year=="2011" & fish$Month %in% c("10", "11"), 4, fish$Session) 

# 2012
fish$Session <- ifelse(fish$Year=="2012" & fish$Month %in% c("01"), 5, fish$Session)
fish$Session <- ifelse(fish$Year=="2012" & fish$Month %in% c("05"), 6, fish$Session) 
fish$Session <- ifelse(fish$Year=="2012" & fish$Month %in% c("07"), 7, fish$Session)       
fish$Session <- ifelse(fish$Year=="2012" & fish$Month %in% c("10", "11"), 8, fish$Session) 

# 2013
fish$Session <- ifelse(fish$Year=="2013" & fish$Month %in% c("01", "02"), 9, fish$Session)
fish$Session <- ifelse(fish$Year=="2013" & fish$Month %in% c("04", "05"), 10, fish$Session) 
fish$Session <- ifelse(fish$Year=="2013" & fish$Month %in% c("07"), 11, fish$Session)       
fish$Session <- ifelse(fish$Year=="2013" & fish$Month %in% c("10", "11"), 12, fish$Session) 

# 2014
fish$Session <- ifelse(fish$Year=="2014" & fish$Month %in% c("01", "02"), 13, fish$Session)
fish$Session <- ifelse(fish$Year=="2014" & fish$Month %in% c("04", "05"), 14, fish$Session) 
fish$Session <- ifelse(fish$Year=="2014" & fish$Month %in% c("08"), 15, fish$Session)       
fish$Session <- ifelse(fish$Year=="2014" & fish$Month %in% c("10", "11", "12"), 16, fish$Session) 

# 2015
fish$Session <- ifelse(fish$Year=="2015" & fish$Month %in% c("01", "02"), 17, fish$Session)
fish$Session <- ifelse(fish$Year=="2015" & fish$Month %in% c("04", "05"), 18, fish$Session) 
fish$Session <- ifelse(fish$Year=="2015" & fish$Month %in% c("07", "08"), 19, fish$Session)       

# 2016
fish$Session <- ifelse(fish$Year=="2016" & fish$Month %in% c("08"), 20, fish$Session)

# 2022
fish$Session <- ifelse(fish$Year=="2022" & fish$Month %in% c("05"), 21, fish$Session) # Different dates in Maracas bc electrofishing broke (OK)

# 2023
fish$Session <- ifelse(fish$Year=="2023" & fish$Month %in% c("05"), 22, fish$Session)


# Seasons:
fish$Season <- ifelse(fish$Month %in% c("01", "02"), 2, NA)
fish$Season <- ifelse(fish$Month %in% c("04", "05", "06"), 3, fish$Season)
fish$Season <- ifelse(fish$Month %in% c("07", "08"), 4, fish$Season)
fish$Season <- ifelse(fish$Month %in% c("10", "11", "12"), 1 , fish$Season) 

# add sessions and seasons


# step 6 -----------------------------------------------------------------------
reps <- fish %>% group_by(Year, Month, Site) %>% summarise(nDay=n_distinct(Day), nID=n_distinct(sampleID)) # Repeatabilities
reps2 <- fish %>% group_by(Year, Month, Site) %>% distinct(Day)

fish <- fish[!fish$Year=="2010",]                     # rm session 0

length(unique(fish$Species[fish$Site=="Maracas u" & fish$Year=="2011" & fish$Season=="3" & fish$Day=="01"])) # 11, Pick
length(unique(fish$Species[fish$Site=="Maracas u" & fish$Year=="2011" & fish$Season=="3" & fish$Day=="03"])) # 9

length(unique(fish$Species[fish$Site=="Maracas d" & fish$Year=="2011" & fish$Season=="3" & fish$Day=="01"])) # 11, Pick
length(unique(fish$Species[fish$Site=="Maracas d" & fish$Year=="2011" & fish$Season=="3" & fish$Day=="03"])) # 8

length(unique(fish$Species[fish$Site=="Acono u" & fish$Year=="2011" & fish$Season=="4" & fish$Day=="12"])) # 6
length(unique(fish$Species[fish$Site=="Acono u" & fish$Year=="2011" & fish$Season=="4" & fish$Day=="15"])) # 7, Pick

length(unique(fish$Species[fish$Site=="Acono d" & fish$Year=="2011" & fish$Season=="4" & fish$Day=="12"])) # 7
length(unique(fish$Species[fish$Site=="Acono d" & fish$Year=="2011" & fish$Season=="4" & fish$Day=="15"])) # 9, Pick

fish <- fish[!(fish$Site %like% "Maracas" & fish$Year=="2011" & fish$Season=="3" & fish$Day=="03"),]# rm repeatability
fish <- fish[!(fish$Site %like% "Acono" & fish$Year=="2011" & fish$Season=="4" & fish$Day=="12"),]  # rm repeatability


# rm unecessary timesteps and repeatabilities


# step 7: ----------------------------------------------------------------------
keysamplingdates <- fish[,names(fish) %in% c("Session", "Day", "Month", "Year", "Site", "Disturbed")]
keysamplingdates <- distinct(keysamplingdates) # 352
save(keysamplingdates, file="keysamplingdates.RData")


# step 8: ----------------------------------------------------------------------
check <- fish %>% group_by(Site, Year, Month, Day) %>% summarise(n=n_distinct(sampleID)) # OK, always 1

# NOTE for below: Rest guppies sampled when 1000 in seen: 
# View(fish[fish$Species=="Poecilia reticulata" & fish$number.seen>0,])
# The below data were updated when entering the data in May-June 2023.

fish[fish$Site=="Quare u" & fish$sampleID==172 & fish$Species=="Poecilia reticulata",]        # Leave as is
fish[fish$Site=="Upper Aripo d" & fish$sampleID==361 & fish$Species=="Poecilia reticulata",]  # OK
fish[fish$Site=="Acono d" & fish$sampleID==377 & fish$Species=="Poecilia reticulata",]        # OK
fish[fish$Site=="Quare u" & fish$sampleID==390 & fish$Species=="Poecilia reticulata",]        # OK

# Other checks: 
# Rest sampled fish when big numbers in seen (2022-2023). 
# View(fish[fish$number.seen > 10,]) #Keep as is for now

# Check D & U well assigned:
sort(unique(fish$Disturbed))
sort(unique(fish$Site[fish$Disturbed=="disturbed"])) # OK
sort(unique(fish$Site[fish$Disturbed=="undisturbed"])) # OK

# Typos in latin names:
sort(unique(fish$Species))
fish$Species <- tolower(fish$Species) # OK
sort(unique(fish$Species))

# Remove unidentified record:
sum(fish$Species %in% c("a (possibly a. hartii)"))  # 1
fish <- fish[!fish$Species %in% c("a (possibly a. hartii)"),]

# check whether seen & observed make sense:
range(fish$number.caught[fish$number.seen>0]) # 0, OK
range(fish$number.seen[fish$number.caught>0]) # 0, OK

# Check weights: 
fish$weight <- as.numeric(fish$weight)
sum(is.na(fish$weight)) # 0, OK
sum(fish$weight==0)     # 0, OK
sum(fish$weight=="")    # 0, OK


# other mixed checks


# step 9 (taxonomy):------------------------------------------------------------
# Change P. sphenops to P. sp & O mossambicus to O sp:
sum(fish$Species=="poecilia sphenops") # 1
sum(fish$number.caught[fish$Species=="poecilia sphenops"]) # 1
fish$Species[fish$Species=="poecilia sphenops"] <- "poecilia sp"
fish$Species[fish$Species=="oreochromis mossambicus"] <- "oreochromis sp"
sort(unique(fish$Species))


sum(fish$Species %in% c("atya spp.", "eudaniela garmani",
                        "macrobrachium acanthurus", "macrobrachium carcinus",
                        "macrobrachium crenulatum", "macrobrachium faustinum",
                        "macrobrachium jelskii", "macrobrachium olfersii",
                        "macrobrachium spp."))  # 580 records
names(fish)



fish$Genus <- str_split_fixed(fish$Species, " ", 2)[,1]
fish$Species <- str_split_fixed(fish$Species, " ", 2)[,2]
fish$Genus <- str_to_title(fish$Genus)

sort(unique(fish$Genus))
sort(unique(fish$Species))

fish$Species[fish$Genus %in% "Macrobrachium"] <- "spp" # pool all records of macrobranchium
fish$Genus[fish$Genus=="Agonostomus"] <- "Dajaus"
fish$Genus[fish$Genus=="Rivulus"] <- "Anablepsoides"
fish$Species[fish$Species=="riisei "] <- "riisei"
unique(fish$Species[fish$Genus=="Odontostilbe"])
fish$Species[fish$Genus=="Odontostilbe"] <- "pulchra"
fish$Species[fish$Species=="spp."] <- "spp"

sort(unique(fish$Species))
sort(unique(fish$Genus))

fish$Taxa <- paste0(fish$Genus, " ", fish$Species)    # add Taxa column
fish <- fish[!fish$Taxa=="Poecilia sp",]              # rm unknown record

# taxonomic updates


################################################################################
# Create Datasets: =============================================================

fish$TotalB <- ifelse(fish$number.seen > 0, (fish$weight * fish$number.seen), (fish$number.caught * fish$weight))
fish$TotalA <- ifelse(fish$number.seen == 0, fish$number.caught, fish$number.seen) # very large number of guppies, appears accurate to rawdata!

unaggFishMacroInv2023 <- fish                   # unaggregated dataset
range(unaggFishMacroInv2023$TotalA)
range(unaggFishMacroInv2023$TotalB)

sum(fish$number.caught[fish$number.seen>0]>0)    # 0
sum(fish$number.seen[fish$number.caught>0]>0)    # 0 
sum(fish$number.seen==0 & fish$number.caught==0) # 0

aggFishMacroInv2023 <- fish %>% group_by(Site, Day, Month, Year, Genus, Species, Session, Season) %>%
  summarise(Biomass=sum(TotalB), Abundance=sum(TotalA)) # aggregated up to 2023
aggFishMacroInv2023 <- as.data.frame(aggFishMacroInv2023)
aggFishMacroInv2010_15 <- as.data.frame(subset(aggFishMacroInv2023, aggFishMacroInv2023$Year < 2016)) # aggregated up to 2015

range(aggFishMacroInv2023$Abundance)
range(aggFishMacroInv2010_15$Abundance)
range(aggFishMacroInv2023$Biomass)
range(aggFishMacroInv2010_15$Biomass)                   

save(unaggFishMacroInv2023, file=paste0(formatted_data_path, "/unaggFishMacroInv2023.RData"))  # doesn't contain pooled data for macrobranchium
save(aggFishMacroInv2023, file=paste0(formatted_data_path, "/aggFishMacroInv2023.RData"))
save(aggFishMacroInv2010_15, file=paste0(formatted_data_path, "/aggFishMacroInv2010_15.RData"))

str(unaggFishMacroInv2023)
str(aggFishMacroInv2023)
str(aggFishMacroInv2010_15)


################################################################################
# NOTES FISH:###################################################################
# 1 ) There are two repeatabilities in this series: 12th and 15th August 2011 
# Acono and 1st and 3rd June 2011 Maracas (step 6)
# 2) Picked repeatabilities: most parsimonious option --> session when more species were reported (step 6)
# 3) multiple taxonomic updates made (step 9)


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