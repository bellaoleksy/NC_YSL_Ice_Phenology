library(tidyverse)
library(ggpubr)
library(mgcv)
library(gratia)
library(plyr)
library(vegan)
library(Matrix)
library(corrplot)
library(itsadug)
library(broom)
library(huxtable)
library(lme4)
library(lmerTest)
library(ggeffects)
library(officer)
library(dataRetrieval)
library(lemon)
library(Hmisc)

#Look into QDO, PDO, ENSO
source("0_functions.R")
summarize <- dplyr::summarize
rename <- dplyr::rename

#Ice On

Timing<-read.csv(file="Data/R/YSL_Ice.csv",header=T,sep=",")
Weather<-read.csv(file="Data/R/Yellowstone_Snow_Rain.csv",header=T,sep=",")

yellowstone_on <- Timing %>%
  select(Year, IceOnDate, IceOnJulian) %>%
  mutate(IceOnJulian = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                 TRUE ~ IceOnJulian),
         IceOnDate = ymd(parse_date_time(paste(Year, IceOnJulian), orders = "yj")),
         water_year = calcWaterYear(IceOnDate)) %>%
  rename(Year_on = Year) %>%
  group_by(water_year) %>%
  slice(which.max(IceOnDate))

#Ice off
yellowstone_off <- Timing %>%
  select(Year, IceOffDate, IceOffJulian) %>%
  mutate(IceOffDate = ymd(parse_date_time(paste(Year, IceOffJulian), orders = "yj")),
         water_year = calcWaterYear(IceOffDate)) %>%
  dplyr::rename(Year_off = Year)

#Combine for full phenology
yellowstone_phenology <- full_join(yellowstone_off, 
                                   yellowstone_on) %>%
  mutate(j_on_wy = hydro.day(IceOnDate),
         j_off_wy = hydro.day(IceOffDate),
         ice_days = j_off_wy - j_on_wy)

min(yellowstone_phenology$ice_days, na.rm=TRUE)
max(yellowstone_phenology$ice_days, na.rm=TRUE)
median(yellowstone_phenology$ice_days, na.rm=TRUE)

min(yellowstone_phenology$j_on_wy, na.rm=TRUE) # december 2
max(yellowstone_phenology$j_on_wy, na.rm=TRUE) # january 28
median(yellowstone_phenology$j_on_wy, na.rm=TRUE) #december 26

min(yellowstone_phenology$IceOffDate, na.rm=TRUE) 
max(yellowstone_phenology$IceOffDate, na.rm=TRUE) 
median(yellowstone_phenology$IceOffJulian, na.rm=TRUE) #may 2024




#Annual weather by water-year

#Need to know how many water-years are have no data for certain variables
# wx_QAQC <- Weather %>%
#   mutate(
#     Date = mdy(Date),
#     year = year(Date),
#     month = month(Date),
#     month_name = month(Date, label = TRUE),
#     water_year = dataRetrieval::calcWaterYear(Date)) %>%
#   select(-SnowDepth.mm,
#          -contains(".F"),
#          -contains(".in")) %>%
#   pivot_longer(max.C:snow.mm) %>%
#   group_by(water_year, name, month) %>% 
#   summarise_all(~sum(is.na(.))) %>%
#   filter(value > 15) %>%
#   select(water_year, name, month, value) %>%
#   pivot_wider(names_from = name,
#               values_from = value)
# # Filter OUT all the year-months where over half the days are missing Wx data
# 
# rain.mm_filtered <- Weather %>%
#   mutate(
#     Date = mdy(Date),
#     year = year(Date),
#     month = month(Date),
#     month_name = month(Date, label = TRUE),
#     water_year = dataRetrieval::calcWaterYear(Date)) %>%
#   select(Date, month:water_year, rain.mm) %>%
#   filter(water_year %in% wx_QAQC$water_year &
#            month %in% wx_QAQC$month)
# 
# snow.mm_filtered <- Weather %>%
#   mutate(
#     Date = mdy(Date),
#     year = year(Date),
#     month = month(Date),
#     month_name = month(Date, label = TRUE),
#     water_year = dataRetrieval::calcWaterYear(Date)) %>%
#   select(Date, month:water_year, snow.mm) %>%
#   filter(water_year %in% wx_QAQC$water_year &
#            month %in% wx_QAQC$month)
# 
#            # rain.mm %in% wx_QAQC$rain.mm &
#            # max.C %in% wx_QAQC$max.C &
#            # mean.C %in% wx_QAQC$mean.C &
#            # min.C %in% wx_QAQC$min.C)
# 
# 
yellowstone_wx_wy <- Weather %>%
  mutate(
    Date = mdy(Date),
    year = year(Date),
    month = month(Date),
    month_name = month(Date, label = TRUE),
    water_year = dataRetrieval::calcWaterYear(Date)) %>%
  group_by(water_year) %>%
  summarize(AnnualMax=sum(max.C,na.rm=T),
            AnnualMin=sum(min.C,na.rm=T),
            AnnualRain=sum(rain.mm,na.rm=T),
            AnnualSnow=sum(snow.mm,na.rm=T))

#Get rid of Inf values and replace with NAs
na_strings <- c(-Inf,Inf)
yellowstone_wx_wy <- yellowstone_wx_wy %>% naniar::replace_with_na_all(condition = ~.x %in% na_strings)
# 
# Weather_2007 <- Weather %>%
#   mutate(Date=mdy(Date), 
#          water_year=calcWaterYear(Date)) %>%
#   filter(water_year==2007)
# #There's barely any data for this water_year... some November dates but nothing into the spring. 
# #How much more data could be missing?
# 
# naniar::vis_miss(Weather)
# #There's a LOT of missing data for SnowDepth.mm particularly toward the end of the dataset
# #It makes me not want to rely on this too heavily, especially since we see an *apparent*
# #decline in max snow depth later in the record. It's probably just because the data are spotty.
# 
# #Visualize the missing data over time
# ggplot(Weather %>%
#          mutate(Date=mdy(Date), 
#                 water_year=calcWaterYear(Date)), 
#        aes(x = Date, 
#            y = SnowDepth.mm)) + 
#   naniar::geom_miss_point(alpha=0.05)
# 
# #Do we have that issue with say snow.mm?
# ggplot(Weather %>%
#          mutate(Date=mdy(Date), 
#                 water_year=calcWaterYear(Date)), 
#        aes(x = Date, 
#            y = snow.mm)) + 
#   naniar::geom_miss_point(alpha=0.05)
# #Not quite as bad... 
# 
# 
# #Season weather by water-year
yellowstone_wx_seasons <- Weather %>%
  mutate(
    Date = mdy(Date),
    year = year(Date),
    month = month(Date),
    month_name = month(Date, label = TRUE),
    water_year = dataRetrieval::calcWaterYear(Date),
    #This is a little janky but we want to pretend that Sept is part of the same
    #water year as the next month (Oct) for the purposes of summarizing the data below
    season = case_when(Month %in% c("Mar","Apr","May") ~ "Spring",
                       Month %in% c("Jun","Jul","Aug","Sep") ~ "Summer",
                       Month %in% c("Oct","Nov") ~ "Fall",
                       Month %in% c("Dec","Jan","Feb") ~ "Winter")) %>%
  group_by(water_year, season) %>%
  dplyr::summarize(
    Max = sum(max.C, na.rm = T),
    Min = sum(min.C, na.rm = T),
    Rain = sum(rain.mm, na.rm = T),
    Snow = sum(snow.mm, na.rm = T),
    TempSum = sum(mean.C, na.rm = T),
    SnowDepth = max(SnowDepth.mm, na.rm = T)
  ) %>%
  #Convert dataframe from long to wide format
  pivot_wider(
    names_from = "season",
    names_sep = "",
    names_glue = "{season}{.value}",
    values_from = c(Max:SnowDepth)
  )

#Monthly weather by water-year
yellowstone_wx_monthly <- Weather %>%
  mutate(
    Date = mdy(Date),
    year = year(Date),
    month = month(Date),
    month_name = month(Date, label = TRUE),
    water_year = dataRetrieval::calcWaterYear(Date)) %>%
  group_by(water_year, month_name) %>%
  dplyr::summarize(  
    Max = sum(max.C, na.rm = T),
    Min = sum(min.C, na.rm = T),
    Rain = sum(rain.mm, na.rm = T),
    Snow = sum(snow.mm, na.rm = T),
    TempSum = sum(mean.C, na.rm = T),
    SnowDepth = max(SnowDepth.mm, na.rm = T)
  ) %>%
  #Convert dataframe from long to wide format
  pivot_wider(
    names_from = "month_name",
    names_sep = "",
    names_glue = "{month_name}{.value}",
    values_from = c(Max:SnowDepth)
  ) 


#Join them all together
yellowstone_full <- yellowstone_phenology %>%
  left_join(., yellowstone_wx_wy) %>%
  left_join(., yellowstone_wx_monthly) %>%
  left_join(., yellowstone_wx_seasons)
#Get rid of Inf values again
yellowstone_full <- yellowstone_full %>% naniar::replace_with_na_all(condition = ~.x %in% na_strings)


# CORR - ICE ON -----------------------------------------------------------

# Have to narrow down our candidate variables before modeling.
# Anything spurious? Colinear?

## Here I am looking at the massive list of potential predictors of ice-in or
## ice-out and seeing which ones are highly correlated (>0.7)

res2 <- rcorr(as.matrix(yellowstone_full[,11:ncol(yellowstone_full)]))
res2


#Create dataframe of pearson r and p-values
YSL_correlations<-flattenCorrMatrix(res2$r, res2$P)

#In an effort to limit the number of covariates, here I am looking at all correlations > 0.7
#and then below making a choice on which of the pair I should keep.
YSL_collinear_trim <- YSL_correlations %>%
  filter(abs(cor) > 0.7) %>%
  arrange(desc(cor)) 

head(YSL_collinear_trim)

#Look for potentially spurious correlations
res3 <- rcorr(as.matrix(yellowstone_full[,8:ncol(yellowstone_full)]))
YSL_ice_correlations<-flattenCorrMatrix(res3$r, res3$P)
YSL_ice_correlations_trim <-YSL_ice_correlations%>%
  filter(row %in% c("j_on_wy","j_off_wy","ice_days")) %>%
  filter(p < 0.05) 
#Quite a lot of spurious associates with ice-on

#Extract variable names
IceOnVars_full <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p < 0.05) %>% 
  arrange(row) %>%
  filter(row=="j_on_wy") %>%
  filter(!grepl('Feb|Spring|Jun|j_off_wy|ice_days', column))  %>% #spurious likely relationships
  arrange(p) %>%
  pull(column)

#We need to narrow down the IceOffVars. How many are collinear? 
YSL_collinear_trim %>% 
  filter(row %in% IceOnVars_full) %>%
  arrange(row)

IceOnVars <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p < 0.05) %>% 
  arrange(row) %>%
  filter(row=="j_on_wy") %>%
  filter(!grepl('Feb|Spring|Jun|j_off_wy|ice_days|Depth', column)) %>%
  filter(!column %in% c("FallSnow", "WinterMax", "WinterMin")) %>%#corr with NovSnow
  arrange(p) %>%
  pull(column)


IceOffVars_full <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p <= 0.05) %>% 
  arrange(row) %>%
  filter(row=="j_off_wy") %>%
  filter(!grepl('Jun|Nov|ice_days', column)) %>%  #spurious likely relationships
  arrange(desc(abs(cor)))

#We need to narrow down the IceOffVars. How many are collinear? 
YSL_collinear_trim %>% 
  filter(row %in% IceOffVars) %>%
  arrange(row)
IceOffVars_full

#Trim IceOffVars
IceOffVars <- IceOffVars_full %>%
  filter(!column %in% c("AprMin", #corr w SpringMin, AprTempSum
                        "SpringTempSum", #corr w AprTempSum
                        "SpringMax",#corr w AprMax
                        "AnnualMax")) %>% #seems spurious-- correlates strong with summer
  pull(column)

#Ice days variables
IceDaysVars <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p <= 0.01) %>% 
  arrange(row) %>%
  filter(row=="ice_days") %>%
  filter(!grepl('June', column)) %>%  #spurious likely relationships
  pull(column)


# ~ MODELS - ICE ON ----------------------------------------------------------------

### I added Family Gamma here since observations are always > 0
mod0_iceOn_new <- gam(j_on_wy ~ s(water_year),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod0_iceOn_new)
# summary(mod0_iceOn)
draw(mod0_iceOn_new)
# draw(mod0_iceOn)
appraise(mod0_iceOn_new)
# appraise(mod0_iceOn)
#These results are different as with Lusha's code, but still no trend

### Start with all fall variables that were highly correlated
IceOnVars

mod1_iceOn_new <- gam(j_on_wy ~ s(WinterTempSum) + s(FallMax) +
                        s(FallRain) + s(FallTempSum) + s(WinterRain)  + s(NovSnow),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod1_iceOn_new)
# summary(mod1_iceOn)
draw(mod1_iceOn_new)
# draw(mod1_iceOn)
appraise(mod1_iceOn_new)
# appraise(mod1_iceOn_new)

mod2_iceOn_new <- gam(j_on_wy ~ s(WinterTempSum) + s(FallMax) +
                        s(FallRain)  + s(WinterRain)  + s(NovSnow),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod2_iceOn_new)
draw(mod2_iceOn_new)
appraise(mod2_iceOn_new)


mod3_iceOn_new <- gam(j_on_wy ~  s(WinterTempSum) + s(FallMax) +
                        s(FallRain)   + s(NovSnow),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod3_iceOn_new)
draw(mod3_iceOn_new)
appraise(mod3_iceOn_new)

mod4_iceOn_new <- gam(j_on_wy ~  s(FallMax) +
                        s(FallRain)   + s(NovSnow),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod4_iceOn_new)
draw(mod4_iceOn_new)
appraise(mod4_iceOn_new)


compareML(mod2_iceOn_new, mod3_iceOn_new) 
compareML(mod3_iceOn_new, mod4_iceOn_new) 




# ~ MODELS - ICE OFF----------------------------------------------------------------

mod0_iceOff_new <- gam(j_off_wy ~ s(water_year),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceOff_new)
draw(mod0_iceOff_new)
appraise(mod0_iceOff_new)
#These are the same as with Lusha's code


IceOffVars

mod1_iceOff_new <- gam(j_off_wy ~ s(AprTempSum) + s(AprMax) + s(MayTempSum) +
                         s(MayMin) + s(SpringRain) + s(SpringMin) + s(SpringSnow),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod1_iceOff_new)
draw(mod1_iceOff_new)
appraise(mod1_iceOff_new)
#Identical. Great.

#Mod 2
mod2_iceOff_new <- gam(j_off_wy ~  s(AprTempSum) + s(AprMax) + s(MayTempSum) +
                     s(MayMin) + s(SpringMin) + s(SpringSnow),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod2_iceOff_new)
draw(mod2_iceOff_new)
appraise(mod2_iceOff_new)

#Mod 3
mod3_iceOff_new <- gam(j_off_wy ~  s(AprTempSum) + s(AprMax)  +
                         s(MayMin) + s(SpringMin) + s(SpringSnow),
                       family=Gamma(link="log"),
                       data = yellowstone_full,
                       # correlatiOff = corCAR1(form = ~ Year),
                       method = "REML")
summary(mod3_iceOff_new)
draw(mod3_iceOff_new)
appraise(mod3_iceOff_new)

#Mod 4
mod4_iceOff_new <- gam(j_off_wy ~  s(AprMax)  +
                         s(MayMin) + s(SpringMin) + s(SpringSnow),
                       family=Gamma(link="log"),
                       data = yellowstone_full,
                       # correlatiOff = corCAR1(form = ~ Year),
                       method = "REML")
summary(mod4_iceOff_new)
draw(mod4_iceOff_new)
appraise(mod4_iceOff_new)


compareML(mod3_iceOff_new, mod4_iceOff_new)
compareML(mod3_iceOff_new, mod2_iceOff_new)


# ~ ANNUAL Max snow depth ----------------------------------------------------------------


mod0_maxSnowDepth_new <- gam(SnowDepth ~ s(water_year),
                         family=Gamma(link="log"),
                         data = yellowstone_wx_wy,
                         # correlation = corCAR1(form = ~ Year),
                         method = "REML")
summary(mod0_maxSnowDepth_new)
summary(mod0_maxSnowDepth)
draw(mod0_maxSnowDepth_new)
draw(mod0_maxSnowDepth)
appraise(mod0_maxSnowDepth_new)
appraise(mod0_maxSnowDepth)
