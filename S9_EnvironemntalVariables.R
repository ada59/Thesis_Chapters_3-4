################################################################################
# S9 Temperature & rainfall as predictors of α diversity
# AFE
# update July 2024
################################################################################

# Libraries: ===================================================================
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(AICcmodavg)
library(broom)
library(patchwork)
library(readxl)
library(data.table)
library(lubridate)
library(corrplot)
library(psych)
library(car)

rm(list=ls())

table_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Tables"
plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"


################################################################################
# READ DATA ====================================================================
piarco <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods/RawData_Trinidad/Environmental/piarco_metoffice.xlsx"
multiplesheets <- function(fname) {
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  names(data_frame) <- sheets
  return(data_frame)
}


lt <- multiplesheets(piarco)
lt <- lapply(lt, function(x) {names(x) <- gsub(" ", "", names(x));x})
lt <- lapply(lt, function(x) {x$TempMean <- (x$TempMax + x$TempMin)/2;x}) # av. temp daily
lt <- lapply(lt, function(x) {x$Date <- paste0(x$Day, "-", x$Month, "-", x$Year);x})
lt <- lapply(lt, function(x) {x$Date <- lubridate::dmy(x$Date);x})


dt <- as.data.frame(do.call(rbind, lt))
sum(dt$Precip=="TR") # 4 in Dec 2017
dt$Precip[dt$Precip=="TR"] <- NA
sum(is.na(dt$Precip))
rownames(dt) <-NULL

check <- dt %>% group_by(Year) %>% summarise(nmonth=n_distinct(Month)) # OK


################################################################################
# DAILY & MONTHLY ==============================================================
str(dt)
dt$Precip <- as.numeric(dt$Precip)
pairs.panels(dt[,c(4:7)])
dtav <- dt %>% group_by(Year, Month) %>% summarise(av_Precip=mean(Precip, na.rm=T),
                                                   t_Precip=sum(Precip, na.rm=T),
                                                   av_TempMax=mean(TempMax),
                                                   av_TempMean=mean(TempMean),
                                                   av_TempMin=mean(TempMin)) # OK
pairs.panels(dtav[,c(3:7)])


# Plot trends monthly ----------------------------------------------------------
dt_long <- gather(dt, key="Variable", value="Value", -c(1:3,8))
dt_long <- dt_long %>% mutate(YearMonth = ym(paste(Year, Month, sep = "-")))
dt_long$YearMonth <- sub("-\\d{2}$", "", dt_long$YearMonth)
dt_long <- dt_long[!dt_long$Variable == "TempMean",]
dt_long$Variable <- plyr::revalue(dt_long$Variable, c("Precip"="Precipitation (mm)",
                                                      "TempMax"="Max Temperature (°C)",
                                                      "TempMin"="Min Temperature (°C)"))


dt_long$Season <- ifelse(dt_long$Month %in% c(1,2,3,4,5), "Dry", "Rainy")
season_boundaries <- dt_long %>%
  group_by(Year, Season) %>%
  summarize(Start_Date = as.character(min(YearMonth)),
            End_Date = as.character(max(YearMonth))) %>%
  ungroup()    # shading seasons
#season_boundaries$End_Date <- gsub("05", "06", season_boundaries$End_Date) # to avoid gaps in plot
#season_boundaries$Start_Date <- gsub("01", "12", season_boundaries$End_Date) # to avoid gaps in plot
str(season_boundaries)

shaded_regions_fishdata <- data.frame(
  xmin = c("2010-12","2016-07", "2022-04"),
  xmax = c("2015-09", "2016-09", "2022-06"),
  ymin = -Inf,
  ymax = Inf
)              # shading fish sampling periods
str(shaded_regions_fishdata)

shaded_regions_fishdata_inv <- data.frame(
  xmin = c("2010-01","2015-09", "2016-09", "2022-06"),
  xmax = c("2010-12", "2016-07", "2022-04", "2022-12"),
  ymin = -Inf,
  ymax = Inf
)              # shading fish sampling periods
str(shaded_regions_fishdata_inv)

unique_yearmonths <- sort(unique(dt_long$YearMonth))
dt_long$YearMonth <- factor(dt_long$YearMonth, levels = unique_yearmonths)
(trends_monthly <- ggplot(dt_long, aes(x = YearMonth, y = Value,  group = YearMonth)) +
  geom_boxplot() +
  #scale_fill_manual(breaks =c("Max Temperature (°C)", "Min Temperature (°C)", "Precipitation (mm)"), values=c("#D55E00", "#56B4E9", "#009E73")) +
  geom_rect(data = shaded_regions_fishdata_inv, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill="darkgrey",alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 4) +
  labs(title = "",
       x = "Year-Month",
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
        axis.title.x = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 15),
        legend.position = "none")) 
(trends_monthly <- trends_monthly +
    geom_rect(data = season_boundaries, aes(xmin = Start_Date, xmax = End_Date, ymin = -Inf, ymax = Inf, fill=Season), alpha = 0.1, inherit.aes = FALSE))
ggsave(paste0(plot_path, "/Chapter5/S9_C5_MonthyClimatic.png"), trends_monthly, width = 14, height = 6.5)
# warnings OK, simply removes a 4 daily observations in 2017 for which there's no data


################################################################################
# QUARTERLY ====================================================================
dt_long$Quaterly <- ifelse(dt_long$Month %in% c(1,2,3), 1, NA)
dt_long$Quaterly <- ifelse(dt_long$Month %in% c(4,5,6), 2, dt_long$Quaterly)
dt_long$Quaterly <- ifelse(dt_long$Month %in% c(7,8,9), 3, dt_long$Quaterly)
dt_long$Quaterly <- ifelse(dt_long$Month %in% c(10,11,12), 4, dt_long$Quaterly)
dt_long$YearQuarter <- paste0(dt_long$Year, "-", dt_long$Quaterly)
unique_yearquarter <- sort(unique(dt_long$YearQuarter))
dt_long$YearQuarter <- factor(dt_long$YearQuarter, levels = unique_yearquarter)


(trends_quarter <- ggplot(dt_long, aes(x = YearQuarter, y = Value, fill=Variable, group = YearQuarter)) +
    geom_boxplot() +
    facet_wrap(~ Variable, scales = "free_y", nrow = 4) +
    labs(title = "Quarterly Averages of Climatic Variables",
         x = "Year-Month",
         y = "") +
    theme_classic() +
    scale_fill_manual(breaks =c("Max Temperature (°C)", "Min Temperature (°C)", "Precipitation (mm)"), values=c("#D55E00", "#56B4E9", "#009E73")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")) 
# not saved, just used to visualise quarterly patterns.


################################################################################
# SAMPLING SEASON DURATION =====================================================
load("C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Thesis_GeneralMethods/FormattedData_Trinidad/aggFish2024.RData")
aggFish2024 <- aggFish2024 %>%
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))

duration <- aggFish2024 %>%
  group_by(Session) %>%
  summarize(DurationDays = as.numeric(max(Date) - min(Date)) + 1)  # add 1 to include both start and end days
mean(duration$DurationDays) # 22 days (av) / 40 days (max) [DURATION OF A SEASON]
print(duration, n=nrow(duration))

first_days <- aggFish2024 %>%
  group_by(Session) %>%
  summarize(first_day = min(Date))

duration_between_samp <- first_days %>%
  mutate(days_to_next = lead(first_day) - first_day) %>%
  na.omit() 
min(duration_between_samp$days_to_next)  # 55 (May 7 July sampling in 2012)
mean(duration_between_samp$days_to_next[duration_between_samp$Session<19]) # 91


################################################################################
# FORMAT FOR ANALYSIS ==========================================================

# NOTE: ========================================================================
# sensitivity: av. of the 1 months from 1st sampling date, and
# av. 1 months from 1st sampling date & plus up to last sampling date.


# Create temp & rainfall data points: ------------------------------------------

sampling <- aggFish2024 %>%
  group_by(Session) %>%
  summarize(earliest = min(Date),
            latest = max(Date))
sampling$Session <- as.character(sampling$Session)
sampling$Session <- plyr::revalue(sampling$Session, c("21"="46",
                                                      "22"="50",
                                                      "23"="54"))
sampling$Session[sampling$Session=="20"] <- "23"
sampling$Session <- as.numeric(sampling$Session)
sampling <- sampling[sampling$Session < 50,] # no temp & rain data for 2023 & 2024


window <- 30


periods <- lapply(sampling$earliest, function(date) {start_date <- date - window # or sampling 2
end_date <- date - 1
filtered_data <- dt %>%
      filter(Date >= start_date & Date <= end_date)
    return(filtered_data)}) 
names(periods) <- sampling$Session
periods <- lapply(periods, function(x) {x %>% summarise(AvPrecip=mean(Precip),
                                                        AvTempMax=mean(TempMax), 
                                                        AvTempMin=mean(TempMin))}) 
periods <- Map(cbind, periods, "Session"=names(periods))
periods <- as.data.frame(do.call(rbind, periods)) # dataset for main analysis


periods_sen <- lapply(1:length(sampling$earliest), function(i) { # or sampling 2
    start_date <- sampling$earliest[i] 
    end_date <- sampling$latest[i]
    window_start <- start_date - window
    window_end <- end_date # or -1, but for this just end_date OK
    filtered_data <- dt %>%
      filter(Date >= window_start & Date <= window_end)
    return(filtered_data)  
  })
names(periods_sen) <- sampling$Session
periods_sen <- lapply(periods_sen, function(x) {x %>% summarise(AvPrecip=mean(Precip),
                                                            AvTempMax=mean(TempMax), 
                                                            AvTempMin=mean(TempMin))}) 
periods_sen <- Map(cbind, periods_sen, "Session"=names(periods_sen))
periods_sen <- as.data.frame(do.call(rbind, periods_sen)) # dataset for sensitivity analysis

periods$Session <- as.numeric(periods$Session)
periods_sen$Session <- as.numeric(periods_sen$Session)

#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnv.jpg"), width = 500, height = 500)
#pairs.panels(periods[,c(1:3)])      
#dev.off()
#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnvSensitivity.jpg"), width = 500, height = 500)
#pairs.panels(periods_sen[,c(1:3)])       
#dev.off()

#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnv2011_15.jpg"), width = 500, height = 500)
#pairs.panels(periods[periods$Session < 20, c(1:3)])      
#dev.off()
#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnvSensitivity2011_15.jpg"), width = 500, height = 500)
#pairs.panels(periods_sen[periods_sen$Session < 20, c(1:3)])       
#dev.off()

test_quadratic <- lm(AvPrecip ~ AvTempMax + I(AvTempMax^2), data=periods)
summary(test_quadratic)
test_quadratic2 <- lm(AvPrecip ~ AvTempMax + I(AvTempMax^2), data=periods_sen)
summary(test_quadratic2) # OK, both ways.


# Alpha diversity subset: ------------------------------------------------------
Results_AlphaDivRData <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Results_AlphaDivRData"
load(paste0(Results_AlphaDivRData, "/S3_lPA_Alpha.RData"))
str(lPA_Alpha)
lPA_Alpha <- lPA_Alpha[!lPA_Alpha$Order.q == 1,]               # rm order q1
lPA_Alpha <- lPA_Alpha[!lPA_Alpha$SC2 == "Cmax",]              # rm order Cmax
lPA_Alpha <- lPA_Alpha[lPA_Alpha$Group == "Northern Range",]   # keep only NR data
lPA_Alpha <- lPA_Alpha[lPA_Alpha$Session < 50,]                # no climate data 2023 & 2024
lPA_Alpha$Taxa[lPA_Alpha$Taxa=="Fish (Trait)"] <- "FTrait"
range(lPA_Alpha$Estimate)

diversity <- split(lPA_Alpha, list(lPA_Alpha$Order.q, lPA_Alpha$Taxa)) 
org <- function(x){
  x <- c(
    x[grep("Diatoms", names(x))],
    x[grep("Invertebrates", names(x))],
    x[grep("Fish", names(x))],
    x[grep("FTrait", names(x))]
  )
  return(x)
} 

diversity <- org(diversity)


# Combine: ---------------------------------------------------------------------

comb <- lapply(diversity, function(x) {as.data.frame(merge(x, periods, by="Session"))})
comb <- lapply(comb, function(x) {as.data.frame(x[,names(x) %in% c("Estimate", "AvPrecip", "AvTempMax")])})


comb_sen <- lapply(diversity, function(x) {as.data.frame(merge(x, periods_sen, by="Session"))})
comb_sen <- lapply(comb_sen, function(x) {as.data.frame(x[,names(x) %in% c("Estimate", "AvPrecip", "AvTempMax")])})


################################################################################
# MODELS =======================================================================

make_horitz <- function(x){
  data <- as.data.frame(do.call(rbind, x))
  data$Code <- paste0("q", str_split_fixed(rownames(data), "\\.", 3)[,1], " ",
                        str_split_fixed(rownames(data), "\\.", 3)[,2])
  data_h <- data %>%
    pivot_wider(
      names_from = Code,  # Column whose values will become new column names
      values_from = Estimate
    )
  return(data_h)
}
comb_horitz <- make_horitz(comb)
comb_sen_horitz <- make_horitz(comb_sen)

#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnvDiv.jpg"), width = 900, height = 900)
#pairs.panels(comb_horitz, cex.cor=1)
#dev.off()
#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_PairsPanelsEnvDivSen.jpg"), width = 900, height = 900)
#pairs.panels(comb_sen_horitz, cex.cor=1)
#dev.off()

run_models <- function(dt){
  m0 <- lm(Estimate ~ AvPrecip*AvTempMax, data = dt)
  m1 <- lm(Estimate ~ AvPrecip + AvTempMax + I(AvTempMax^2), data = dt) 
  m2 <- lm(Estimate ~ AvTempMax + AvPrecip + I(AvPrecip^2), data = dt) 
  m3 <- lm(Estimate ~ AvPrecip + AvTempMax, data = dt)
  m4 <- lm(Estimate ~ AvTempMax + I(AvTempMax^2), data = dt)
  m5 <- lm(Estimate ~ AvPrecip + I(AvPrecip^2), data = dt)
  m6 <- lm(Estimate ~ AvPrecip, data = dt)
  m7 <- lm(Estimate ~ AvTempMax, data = dt)
  models <- list("Interaction"=m0, 
                 "Main&TempQ"=m1, 
                 "Main&PrecipQ"=m2,
                 "MainLinear"=m3, 
                 "TempQ"=m4, 
                 "PrecipQ"=m5, 
                 "Precip"=m6,
                 "Temp"=m7)
  return(models)
} 

comb_mod <- lapply(comb, function(x) {run_models(dt=x)})          
comb_sen_mod <- lapply(comb_sen, function(x) {run_models(dt=x)})


################################################################################
# MODEL SELECTION ==============================================================

R2_AIC_Tab <- function(dt){
  R2 <- lapply(dt, function(x) {lapply(x, function(y) {round(summary(y)$adj.r.squared, 2)})}) # R2
  R2 <- lapply(R2, function(x) {as.data.frame(do.call(rbind, x))})
  compare <- lapply(dt, function(x) {aictab(x)})                                              # rank by AICc
  compare <- mapply(function(x,y) {x$adjR2 <- y$V1[match(x$Modnames, rownames(y))];x}, compare, R2, SIMPLIFY = F)
  raw_tabs <- do.call(rbind, compare)
  
  raw_tabs$Code <- rownames(raw_tabs) # organise & layout for SM
  rownames(raw_tabs) <- 1:nrow(raw_tabs)
  raw_tabs <- separate(raw_tabs, Code, into = paste0("Column", 1:3), sep = "\\.")          # re-assign names
  raw_tabs <- raw_tabs[, ! names(raw_tabs) %in% c("Column3")]                              # not needed
  names(raw_tabs)[names(raw_tabs)=="Column1"] <- "Orderq"
  names(raw_tabs)[names(raw_tabs)=="Column2"] <- "Taxa"
  raw_tabs <- raw_tabs %>% relocate(c(Taxa, Orderq), .before = Modnames)
  raw_tabs <- raw_tabs[, ! names(raw_tabs)  %in% c("K", "ModelLik", "AICcWt", "Cum.Wt")]   # unnecessary vals
  raw_tabs[,c(4:7)] <- sapply(raw_tabs[,c(4:7)], function(x) {round(x, 2)})                # all rounded
  
  raw_tabs$Taxa <- factor(raw_tabs$Taxa, levels = c("Diatoms", "Invertebrates", 
                                                    "Fish", "FTrait"))                     # correct order
  raw_tabs <- raw_tabs[order(raw_tabs$Taxa), ]                                             # correct order
  raw_tabs <- raw_tabs %>% arrange(Taxa, Orderq)                           # correct order
  formatted_tabs <- raw_tabs
  
  return(formatted_tabs)
} ### consider model averaging if there's time

comb_aic <- R2_AIC_Tab(comb_mod)
comb_sen_aic <- R2_AIC_Tab(comb_sen_mod)

#write.csv(comb_aic, file=paste0(table_path, "/Chapter5/S9_AlphaEnv_AIC_Tables.csv"), row.names = F)
#write.csv(comb_sen_aic, file=paste0(table_path, "/Chapter5/S9_AlphaEnvSen_AIC_Tables.csv"), row.names = F)


comb_best <- subset(comb_aic, comb_aic$Delta_AICc==0)
comb_best <- split(comb_best, f=list(comb_best$Orderq, comb_best$Taxa)) 
comb_best <- lapply(comb_best, function(x) {x$Modnames})                                      
comb_best <- mapply(function(x,y) {x[y]}, comb_mod, comb_best, SIMPLIFY=T)                   # best according to AICc
comb_sum <- lapply(comb_best, function(x) {summary(x)})                                      # summary best

comb_sen_best <- subset(comb_sen_aic, comb_sen_aic$Delta_AICc==0)
comb_sen_best <- split(comb_sen_best, f=list(comb_sen_best$Orderq, comb_sen_best$Taxa)) 
comb_sen_best <- lapply(comb_sen_best, function(x) {x$Modnames})                                      
comb_sen_best <- mapply(function(x,y) {x[y]}, comb_sen_mod, comb_sen_best, SIMPLIFY=T)       # best according to AICc
comb_sen_sum <- lapply(comb_sen_best, function(x) {summary(x)})                              # summary best


################################################################################
# MODEL VALIDATION: ============================================================

# normality & homoscedasticity -------------------------------------------------
comb_sresid <- lapply(comb_best, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)})           # standardised residuals
comb_sen_sresid <- lapply(comb_sen_best, function(x) {(x$residuals - mean(x$residuals))/sd(x$residuals)})   # standardised residuals

namesComb <- names(comb_sresid)
namesComb <- gsub("\\.", " ", namesComb)

namesCombSen <- names(comb_sen_sresid)
namesCombSen <- gsub("\\.", " ", namesCombSen)


#plot_path <- "C:/Users/afe1/OneDrive - University of St Andrews/PHD/0_GLOBAL_THESIS_CHAPTER_GITHUB_REPOSITORIES/Shared-SharedSMs/Plots"
#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_AlphaEnv_NormalityAssumption.jpg"), width = 900, height = 600)
#par(mfrow=c(2,4))
#mapply(function(x,y){
#  qqnorm(x, cex=2,  pch=20, main = y, cex.main=1.7)
#  qqline(x, lty=2, lwd=2)
#  }, comb_sresid, namesComb, SIMPLIFY = T) 
#dev.off()  # assess normality of residuals (or comb_sen)

#jpeg(file=paste0(plot_path,"/Chapter5/S9_C5_AlphaEnv_HomoscedasticityAssumption.jpg"), width = 900, height = 600)
#par(mfrow=c(2,4))
#mapply(function(x,y,z){
# plot(x ~ y$fitted.values, pch=20, cex=2, cex.lab=1.5, main=z, cex.main=1.7)
#}, comb_sresid, comb_best, namesComb, SIMPLIFY = T)
#dev.off()  # assess homoscedasticity  (or comb_sen)
#par(mfrow=c(2,4))
#lapply(comb_best, function(x) {plot(x, which=1)})
#dev.off()


# NOTES: =======================================================================
# Meets both assumptions roughly (better fit for q2 than q0)
# Sen similar but worse fit



################################################################################
# AUTOCORRELATION: =============================================================
lapply(comb_mod, function(x) {lapply(x, function(y) {durbinWatsonTest(y)})})
lapply(comb_sen_mod, function(x) {lapply(x, function(y) {durbinWatsonTest(y)})})

auto_durbin_tab <- function(x){
  auto_p <- x
  auto_p <- as.data.frame(do.call(rbind, auto_p))
  auto_p$Code <- rownames(auto_p) # organise & layout for SM
  rownames(auto_p) <- 1:nrow(auto_p)
  auto_p <- separate(auto_p, Code, into = paste0("Column", 1:3), sep = "\\.")     
  auto_p <- auto_p[, ! names(auto_p) %in% c("Column3", "alternative")]                       
  names(auto_p)[names(auto_p)=="Column1"] <- "Orderq"
  names(auto_p)[names(auto_p)=="Column2"] <- "Taxa"
  auto_p <- auto_p %>% relocate(c(Taxa, Orderq), .before = r)
  auto_p[,c(3:5)] <- sapply(auto_p[,c(3:5)], function(x) {round(as.numeric(x), 3)})           
  return(auto_p)
}

auto <- lapply(comb_best, function(y) {durbinWatsonTest(y)})  # best model
auto <- auto_durbin_tab(auto)
auto_sen <- lapply(comb_sen_best, function(y) {durbinWatsonTest(y)})  # best model (sensitivity)
auto_sen <- auto_durbin_tab(auto_sen)

#write.csv(auto, file=paste0(table_path, "/Chapter5/S9_AlphaEnv_Auto_Tables.csv"), row.names = F)
#write.csv(auto_sen, file=paste0(table_path, "/Chapter5/S9_AlphaEnv_AutoSen_Tables.csv"), row.names = F)


################################################################################
# COEFS BEST: ==================================================================
retrieve_coefs <- function(best){
  coefs <- lapply(best, function(x) {as.data.frame(summary(x)$coefficients)}) 
  coefs <- lapply(coefs, function(x) {x$Term <- rownames(x);x})  
  coefs <- do.call(rbind, coefs)
  coefs$Code <- rownames(coefs)
  rownames(coefs) <- 1:nrow(coefs)
  
  coefs <- separate(coefs, Code, into = paste0("Column", 1:4), sep = "\\.") # warning OK
  names(coefs)[names(coefs)=="Column1"] <- "Orderq"
  names(coefs)[names(coefs)=="Column2"] <- "Taxa"
  names(coefs)[names(coefs)=="Column3"] <- "Model"
  names(coefs)[names(coefs)=="Column4"] <- "Variable"
  coefs <- coefs[,!names(coefs) %in% c("Model", "Term")]
  coefs <- coefs %>% relocate(c(Taxa, Orderq, Variable), .before = Estimate)
  return(coefs)
}


comb_coefs <- retrieve_coefs(comb_best) 
comb_coefs$Variable <- plyr::revalue(comb_coefs$Variable, c("(Intercept)"="Intercept")) 

comb_sen_coefs <- retrieve_coefs(comb_sen_best) 
comb_sen_coefs$Variable <- plyr::revalue(comb_sen_coefs$Variable, c("(Intercept)"="Intercept")) 

write.csv(comb_coefs, file=paste0(table_path, "/Chapter5/S9_C5_AlphaEnv_Best_Coef_Tables.csv"), row.names = F)
write.csv(comb_sen_coefs, file=paste0(table_path, "/Chapter5/S9_C5_AlphaEnv_Best_Coef_Tables_Sen.csv"), row.names = F)

################################################################################
# PREDICT: =====================================================================
# https://biologyforfun.wordpress.com/2014/04/08/interpreting-interaction-coefficient-in-r-part1-lm/

sum_best <- function(x){
  R2adj <- rep(summary(x)$adj.r.squared, times=1000)
  sum <- summary(x)$coefficients
  p_val <- rep(sum[2,4], times=1000)
  summa <- as.data.frame(cbind("R2"=R2adj, "pval"=p_val))
  return(summa)
}

pred_best <- function(bestmodel, rawdata){
  best_x <- data.frame("AvPrecip"=seq(from = min(rawdata$AvPrecip), to = max(rawdata$AvPrecip), length= 1000),
                       "AvTempMax"=seq(from = min(rawdata$AvTempMax), to = max(rawdata$AvTempMax), length= 1000))
  pred <- data.frame(best_x, predict(bestmodel, best_x, interval="confidence"), sum_best(bestmodel))
  return(pred)
}


comb_pred <- mapply(function(x,y) {pred_best(x,y)}, comb_best, comb, SIMPLIFY = F)
comb_sen_pred <- mapply(function(x,y) {pred_best(x,y)}, comb_sen_best, comb_sen, SIMPLIFY = F)


################################################################################
# PLOTS: =======================================================================
comb_pred <- as.data.frame(do.call(rbind, comb_pred))
comb_sen_pred <- as.data.frame(do.call(rbind, comb_sen_pred))

process_df <- function(x){
  x$Code <- rownames(x)
  rownames(x) <- 1:nrow(x)
  
  x <- separate(x, Code, into = paste0("Column", 1:4), sep = "\\.") # warning OK
  names(x)[names(x)=="Column1"] <- "Orderq"
  names(x)[names(x)=="Column2"] <- "Taxa"
  names(x)[names(x)=="Column3"] <- "Model"
  x <- x[,!names(x) %in% c("Column4")]
  x$TaxaOrderq <- paste0(x$Taxa, " " ,x$Orderq)
  x <- x %>% relocate(c(TaxaOrderq, Taxa, Orderq, Model), .before = AvPrecip)
  return(x)
}


comb_pred <- process_df(comb_pred)
comb_pred <- comb_pred[comb_pred$pval < 0.1,]              # linear regressions
sort(unique(comb_pred$TaxaOrderq))                         # "Diatoms 2" & "Invertebrates 2"

comb_sen_pred <- process_df(comb_sen_pred)
comb_sen_pred <- comb_sen_pred[comb_sen_pred$pval < 0.1,]  # linear regressions
sort(unique(comb_sen_pred$TaxaOrderq))                     # "Invertebrates 2" (and worse fits overall)


comb_raw_plots <- as.data.frame(do.call(rbind, comb))
comb_raw_plots$Orderq <- str_split_fixed(rownames(comb_raw_plots), "\\.", 3)[,1]
comb_raw_plots$Taxa <- str_split_fixed(rownames(comb_raw_plots), "\\.", 3)[,2]
comb_raw_plots <- comb_raw_plots[comb_raw_plots$Taxa %in% c("Diatoms", "Invertebrates") & comb_raw_plots$Orderq %in% c(2),]              

comb_raw_sen_plots <- as.data.frame(do.call(rbind, comb_sen))
comb_raw_sen_plots$Orderq <- str_split_fixed(rownames(comb_raw_sen_plots), "\\.", 3)[,1]
comb_raw_sen_plots$Taxa <- str_split_fixed(rownames(comb_raw_sen_plots), "\\.", 3)[,2]
comb_raw_sen_plots <- comb_raw_sen_plots[comb_raw_sen_plots$Taxa == "Invertebrates" & comb_raw_sen_plots$Orderq == 2,]              


(p <- ggplot(data=comb_raw_plots, aes(x=AvPrecip, y=Estimate))+
  labs(x="Average Precipitation (mm)", y="evenness (α q2) cov. Cmin")+
  geom_point(data=comb_raw_plots, aes(x=AvPrecip, y=Estimate), size=1.5, shape=16, alpha=0.2)+
  geom_line(data=comb_pred, aes(x=AvPrecip, y=fit), size=0.8, linetype="solid", color="#1b3765")+
  geom_ribbon(data=comb_pred, aes(x=AvPrecip, y=fit, ymin = lwr, ymax = upr), fill="white", alpha = 0.1, color="#1b3765")+
  theme_classic()+
  ylim(15,25.5)+
  facet_wrap(~ Taxa, nrow=1))


(p_sen <- ggplot(data=comb_raw_sen_plots, aes(x=AvPrecip, y=Estimate))+
    labs(x="Average Precipitation (mm)", y="evenness (α q2) cov. Cmin")+
    geom_point(data=comb_raw_sen_plots, aes(x=AvPrecip, y=Estimate), size=1.5, shape=16, alpha=0.2)+
    geom_line(data=comb_sen_pred, aes(x=AvPrecip, y=fit), size=0.8, linetype="solid", color="#1b3765")+
    ylim(15,25.5)+
    geom_ribbon(data=comb_sen_pred, aes(x=AvPrecip, y=fit, ymin = lwr, ymax = upr), fill="white", alpha = 0.1, color="#1b3765")+
    theme_classic())


ggsave(paste0(plot_path, "/Chapter5/S9_C5_PlotEffect_AvPrecip.png"), p, width = 6, height = 4) 
ggsave(paste0(plot_path, "/Chapter5/S9_C5_PlotEffect_AvPrecip_Sen.png"), p_sen, width = 3, height = 4)   


# End of script ################################################################
################################################################################