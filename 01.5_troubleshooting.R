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

#Custom function to make sure sum() doesn't return 0 for column of NAs
sum2 <- function(x) {
  if(all(is.na(x))){
    c(x[0],NA)} else {
      sum(x,na.rm = TRUE)}
}



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
  summarize(AnnualMax=sum2(max.C),
            AnnualMin=sum2(min.C),
            AnnualRain=sum2(rain.mm),
            AnnualSnow=sum2(snow.mm))

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
    Max = sum2(max.C),
    Min = sum2(min.C),
    Rain = sum2(rain.mm),
    Snow = sum2(snow.mm),
    TempSum = sum2(mean.C),
    SnowDepth = max(SnowDepth.mm)
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
    Max = sum2(max.C),
    Min = sum2(min.C),
    Rain = sum2(rain.mm),
    Snow = sum2(snow.mm),
    TempSum = sum2(mean.C),
    SnowDepth = max(SnowDepth.mm)
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
  filter(!column %in% c("FallSnow", #corr with Nov Snow
                        "WinterMin")) %>% #corr w/ winterMin
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
  filter(!column %in% c("AprTempSum", #corr w AprMax
                        "SpringMax",#corr w MayTempSum
                        "SpringMin")) %>% #corr w SpringTempSum
  filter(!grepl('Annual', column)) %>%
  pull(column)

#Ice days variables
IceDaysVars_full <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p <= 0.01) %>% 
  arrange(row) %>%
  filter(row=="ice_days") %>%
  filter(!grepl('June', column)) %>%  #spurious likely relationships
  arrange(p) %>%
  pull(column)

#We need to narrow down the IceDaysVars How many are collinear? 
YSL_collinear_trim %>% 
  filter(row %in% IceDaysVars_full) %>%
  arrange(row)

# ~ MODELS - ICE ON ----------------------------------------------------------------

### I added Family Gamma here since observations are always > 0
mod0_iceOn <- gam(j_on_wy ~ s(water_year),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod0_iceOn)
# summary(mod0_iceOn)
draw(mod0_iceOn)
# draw(mod0_iceOn)
appraise(mod0_iceOn)
# appraise(mod0_iceOn)
#These results are different as with Lusha's code, but still no trend

### Start with all fall variables that were highly correlated
IceOnVars

mod1_iceOn <- gam(j_on_wy ~ s(NovSnow) + s(DecMin) +
                        s(WinterTempSum) + s(DecTempSum) + s(WinterMax)  + s(JanMin),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod1_iceOn)
# summary(mod1_iceOn)
draw(mod1_iceOn)
# draw(mod1_iceOn)
appraise(mod1_iceOn)
# appraise(mod1_iceOn)

mod2_iceOn <- gam(j_on_wy ~ s(NovSnow) + s(DecMin) + 
                    s(DecTempSum) + s(WinterMax)  + s(JanMin),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod2_iceOn)
draw(mod2_iceOn)
appraise(mod2_iceOn)


mod3_iceOn <- gam(j_on_wy ~  s(NovSnow) + s(DecMin) + 
                     s(WinterMax)  + s(JanMin),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod3_iceOn)
draw(mod3_iceOn)
appraise(mod3_iceOn)

mod4_iceOn <- gam(j_on_wy ~   s(NovSnow) + s(DecMin) + 
                     s(JanMin),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod4_iceOn)
draw(mod4_iceOn)
appraise(mod4_iceOn)


mod5_iceOn <- gam(j_on_wy ~   s(NovSnow) + s(DecMin),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  method = "REML")
summary(mod5_iceOn)
draw(mod5_iceOn)
appraise(mod5_iceOn)

compareML(mod2_iceOn, mod3_iceOn) 
compareML(mod3_iceOn, mod4_iceOn) 
compareML(mod5_iceOn, mod4_iceOn) 




#.....  Ice On Figure -----------------------------------------------------------
#annotate panel letters inside plot
panelLetter.normal <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  hjustvar = c(-0.5) ,
  vjustvar = c(1.5))

summary(mod4_iceOn)

# >>>> Panel A -- Ice On vs. Cumul. Nov. Snow --------------------------------------

new_data <-
  with(yellowstone_full,
       expand.grid(
         NovSnow = seq(
           min(NovSnow, na.rm = TRUE),
           max(NovSnow, na.rm =
                 TRUE),
           length = 200
         ),
         DecMin = median(DecMin, na.rm =
                             TRUE),
         JanMin = median(JanMin, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOn)$linkinv
pred_NovSnow <- predict(mod4_iceOn, new_data, type = "link", se.fit = TRUE)
pred_NovSnow <- cbind(pred_NovSnow, new_data)
pred_NovSnow <- transform(pred_NovSnow, lwr_ci = ilink(fit - (2 * se.fit)),
                          upr_ci = ilink(fit + (2 * se.fit)),
                          fitted = ilink(fit))

pred_NovSnow <- pred_NovSnow %>%
  select(NovSnow, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_NovSnow = lwr_ci,
                upr_ci_NovSnow = upr_ci,
                fitted_NovSnow = fitted)

#Modify axis labels: fed DOY -> dates
breaks_IceOn<-c(60, 80, 100, 120)
# c("30-Nov","20-Dec","10-Jan","30-Jan") #corresponding cal dates


Fig3A <-
  ggplot(pred_NovSnow, aes(x = NovSnow, y = fitted_NovSnow)) +
  geom_ribbon(aes(ymin = lwr_ci_NovSnow, ymax = upr_ci_NovSnow), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=NovSnow,
                             y=j_on_wy,
                             fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative Nov. Snow (mm)",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,
                     labels=c("30-Nov","20-Dec","10-Jan","30-Jan"))+
  # scale_x_continuous(breaks=c(-35,-30,-25,-20,-15),limits=c(-35,-15))+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  theme_pubr(border=TRUE, base_size=8)+
  theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2)) +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="a",
                fontface="bold"))


# # >>>>  Panel B -- Ice On vs. Cumulative min Dec  -------------------------------------

summary(mod4_iceOn)

new_data <-
  with(yellowstone_full,
       expand.grid(
         DecMin = seq(
           min(DecMin, na.rm = TRUE),
           max(DecMin, na.rm =
                 TRUE),
           length = 200
         ),
         NovSnow = median(NovSnow, na.rm =
                            TRUE),
         JanMin = median(JanMin, na.rm=TRUE)
       ))


ilink <- family(mod4_iceOn)$linkinv
pred_DecMin <- predict(mod4_iceOn, new_data, type = "link", se.fit = TRUE)
pred_DecMin <- cbind(pred_DecMin, new_data)
pred_DecMin <- transform(pred_DecMin, lwr_ci = ilink(fit - (2 * se.fit)),
                          upr_ci = ilink(fit + (2 * se.fit)),
                          fitted = ilink(fit))
pred_DecMin <- pred_DecMin %>%
  select(DecMin, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_DecMin = lwr_ci,
                upr_ci_DecMin = upr_ci,
                fitted_DecMin = fitted)

#Modify axis labels: fed DOY -> dates
breaks_IceOn<-c(60, 80, 100, 120)
# c("30-Nov","20-Dec","10-Jan","30-Jan") #corresponding cal dates


Fig3B <-
  ggplot(pred_DecMin, aes(x = DecMin, y = fitted_DecMin)) +
  geom_ribbon(aes(ymin = lwr_ci_DecMin, ymax = upr_ci_DecMin), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=DecMin,
                                        y=j_on_wy,
                                        fill=water_year), 
             shape=21,alpha=0.9)+
  labs(x="Cumulative min. Dec. temperatures (°C)",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,
                     labels=c("30-Nov","20-Dec","10-Jan","30-Jan"))+
  # scale_x_continuous(breaks=labels_IsoMax_fed,labels=c("12-Dec","01-Jan","21-Jan"),limits=c(70,120))+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  theme_pubr(border=TRUE, base_size=8)+
  theme(
    axis.title.y=element_blank(),
    plot.margin=unit(c(0.5,0,0.5,0), "lines"),
    axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="b",
                fontface="bold"))





# # >>>>  Panel C -- Ice On vs. Cumul min Jan -------------------------------------
summary(mod4_iceOn)

new_data <-
  with(yellowstone_full,
       expand.grid(
         JanMin = seq(
           min(JanMin, na.rm = TRUE),
           max(JanMin, na.rm =
                 TRUE),
           length = 200
         ),
         NovSnow = median(NovSnow, na.rm =
                             TRUE),
         DecMin = median(DecMin, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOn)$linkinv
pred_JanMin <- predict(mod4_iceOn, new_data, type = "link", se.fit = TRUE)
pred_JanMin <- cbind(pred_JanMin, new_data)
pred_JanMin <- transform(pred_JanMin, lwr_ci = ilink(fit - (2 * se.fit)),
                           upr_ci = ilink(fit + (2 * se.fit)),
                           fitted = ilink(fit))
pred_JanMin <- pred_JanMin %>%
  select(JanMin, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_JanMin = lwr_ci,
                upr_ci_JanMin = upr_ci,
                fitted_JanMin = fitted)

breaks_IceOn<-c(60, 80, 100, 120)
# c("30-Nov","20-Dec","10-Jan","30-Jan") #corresponding cal dates


Fig3C <-
  ggplot(pred_JanMin, aes(x = JanMin, y = fitted_JanMin)) +
  geom_ribbon(aes(ymin = lwr_ci_JanMin, ymax = upr_ci_JanMin), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=JanMin,
                             y=j_on_wy,
                             fill=water_year), 
             shape=21,alpha=0.9)+
  labs(x="Cumulative min. Jan. temperatures (°C)",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,
                     labels=c("30-Nov","20-Dec","10-Jan","30-Jan"))+
  # scale_x_continuous(breaks=c(0,75,150,225,300,375),limits=c(0,400))+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  theme_pubr(border=TRUE, base_size=8)+
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank(),
    plot.margin=unit(c(0.5,0,0.5,0), "lines"),
    axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="c",
                fontface="bold"))


# ~ MODELS - ICE OFF----------------------------------------------------------------

mod0_iceOff <- gam(j_off_wy ~ s(water_year),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceOff)
draw(mod0_iceOff)
appraise(mod0_iceOff)
#These are the same as with Lusha's code


IceOffVars

mod1_iceOff <- gam(j_off_wy ~ s(AprMax) + s(SpringTempSum) + s(MayTempSum) +
                         s(MayMin) + s(SpringSnow) + s(AprMin) + s(MarMax),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod1_iceOff)
draw(mod1_iceOff)
appraise(mod1_iceOff)

#Mod 2
mod2_iceOff <- gam(j_off_wy ~  s(AprMax) + s(SpringTempSum) + s(MayTempSum) +
                     s(MayMin) + s(SpringSnow) + s(AprMin),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod2_iceOff)
draw(mod2_iceOff)
appraise(mod2_iceOff)

#Mod 3
mod3_iceOff <- gam(j_off_wy ~  s(AprMax) + s(SpringTempSum) +
                     s(MayMin) + s(SpringSnow) + s(AprMin),
                       family=Gamma(link="log"),
                       data = yellowstone_full,
                       # correlatiOff = corCAR1(form = ~ Year),
                       method = "REML")
summary(mod3_iceOff)
draw(mod3_iceOff)
appraise(mod3_iceOff)

#Mod 4
mod4_iceOff <- gam(j_off_wy ~  s(AprMax) + s(SpringTempSum) +
                     s(MayMin) + s(SpringSnow) ,
                       family=Gamma(link="log"),
                       data = yellowstone_full,
                       # correlatiOff = corCAR1(form = ~ Year),
                       method = "REML")
summary(mod4_iceOff)
draw(mod4_iceOff)
appraise(mod4_iceOff)

#Mod 5
mod5_iceOff <- gam(j_off_wy ~  s(AprMax) + s(SpringTempSum) +
                     s(MayMin) ,
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod5_iceOff)
draw(mod5_iceOff)
appraise(mod5_iceOff)


compareML(mod3_iceOff, mod4_iceOff)
compareML(mod3_iceOff, mod2_iceOff)
compareML(mod5_iceOff, mod4_iceOff)


#.....  Ice Off Figure -----------------------------------------------------------

summary(mod4_iceOff)

# >>>> Panel A -- Ice Off vs. AprMax ------------------------------------


new_data <-
  with(yellowstone_full,
       expand.grid(
         AprMax = seq(
           min(AprMax, na.rm = TRUE),
           max(AprMax, na.rm =
                 TRUE),
           length = 200
         ),
         SpringTempSum = median(SpringTempSum, na.rm =
                               TRUE),
         MayMin = median(MayMin, na.rm=TRUE),
         SpringSnow = median(SpringSnow, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOff)$linkinv
pred_AprMax <- predict(mod4_iceOff, new_data, type = "link", se.fit = TRUE)
pred_AprMax <- cbind(pred_AprMax, new_data)
pred_AprMax <- transform(pred_AprMax, lwr_ci = ilink(fit - (2 * se.fit)),
                             upr_ci = ilink(fit + (2 * se.fit)),
                             fitted = ilink(fit))
pred_AprMax <- pred_AprMax %>%
  select(AprMax, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_AprMax = lwr_ci,
                upr_ci_AprMax = upr_ci,
                fitted_AprMax = fitted)

#Modify axis labels: 
breaks_IceOff<-c(210, 225, 240, 255)
# figure out what DOY corresponds to as a date
# c("29-Apr","14-May","29-May","13-Jun")


Fig3D <-
  ggplot(pred_AprMax, aes(x = AprMax, y = fitted_AprMax)) +
  geom_ribbon(aes(ymin = lwr_ci_AprMax, ymax = upr_ci_AprMax), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=AprMax,
                              y=j_off_wy,
                              fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative max. April temperatures (°C)",
  y="Ice-off date")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,
                     labels=c("29-Apr","14-May","29-May","13-Jun"))+
  theme_pubr(border=TRUE, base_size=8)+
  theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2)) +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="d",
                fontface="bold"))


# >>>> Panel B -- Ice Off vs. SpringTempSum -----------------------------------
summary(mod4_iceOff)


new_data <-
  with(yellowstone_full,
       expand.grid(
         SpringTempSum = seq(
           min(SpringTempSum, na.rm = TRUE),
           max(SpringTempSum, na.rm =
                 TRUE),
           length = 200
         ),
         AprMax = median(AprMax, na.rm =
                               TRUE),
         MayMin = median(MayMin, na.rm=TRUE),
         SpringSnow = median(SpringSnow, na.rm=TRUE)
       ))


ilink <- family(mod4_iceOff)$linkinv
pred_SpringTempSum <- predict(mod4_iceOff, new_data, type = "link", se.fit = TRUE)
pred_SpringTempSum <- cbind(pred_SpringTempSum, new_data)
pred_SpringTempSum <- transform(pred_SpringTempSum, lwr_ci = ilink(fit - (2 * se.fit)),
                             upr_ci = ilink(fit + (2 * se.fit)),
                             fitted = ilink(fit))
pred_SpringTempSum <- pred_SpringTempSum %>%
  select(SpringTempSum, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_SpringTempSum = lwr_ci,
                upr_ci_SpringTempSum = upr_ci,
                fitted_SpringTempSum = fitted)

breaks_IceOff<-c(210, 225, 240, 255)
# figure out what DOY corresponds to as a date
# c("29-Apr","14-May","29-May","13-Jun")


Fig3E <-
  ggplot(pred_SpringTempSum, aes(x = SpringTempSum, y = fitted_SpringTempSum)) +
  geom_ribbon(aes(ymin = lwr_ci_SpringTempSum, ymax = upr_ci_SpringTempSum), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=SpringTempSum,
                              y=j_off_wy,
                              fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative spring temperatures (°C)",
       y="Ice off")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,
                     labels=c("29-Apr","14-May","29-May","13-Jun"))+
  theme_pubr(border=TRUE, base_size=8)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "lines"),
        axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="e",
                fontface="bold"))



# >>>> Panel C - Ice off versus MayMin ----------------------------------
summary(mod4_iceOff)


new_data <-
  with(yellowstone_full,
       expand.grid(
         MayMin = seq(
           min(MayMin, na.rm = TRUE),
           max(MayMin, na.rm =
                 TRUE),
           length = 200
         ),
         AprMax = median(AprMax, na.rm =
                               TRUE),
         SpringTempSum = median(SpringTempSum, na.rm=TRUE),
         SpringSnow = median(SpringSnow, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOff)$linkinv
pred_MayMin <- predict(mod4_iceOff, new_data, type = "link", se.fit = TRUE)
pred_MayMin <- cbind(pred_MayMin, new_data)
pred_MayMin <- transform(pred_MayMin, lwr_ci = ilink(fit - (2 * se.fit)),
                            upr_ci = ilink(fit + (2 * se.fit)),
                            fitted = ilink(fit))
pred_MayMin <- pred_MayMin %>%
  select(MayMin, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_MayMin = lwr_ci,
                upr_ci_MayMin = upr_ci,
                fitted_MayMin = fitted)


breaks_IceOff<-c(210, 225, 240, 255)
# figure out what DOY corresponds to as a date
# c("29-Apr","14-May","29-May","13-Jun")


Fig3F <-
  ggplot(pred_MayMin, aes(x = MayMin, y = fitted_MayMin)) +
  geom_ribbon(aes(ymin = lwr_ci_MayMin, ymax = upr_ci_MayMin), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=MayMin,
                              y=j_off_wy,
                              fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative min. May temperatures (C)",
       y="Ice off")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,
                     labels=c("29-Apr","14-May","29-May","13-Jun"))+
  theme_pubr(border=TRUE, base_size=8)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "lines"),
        axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="f",
                fontface="bold"))

# >>>> Panel D - Ice off versus SpringSnow ----------------------------------
summary(mod4_iceOff)


new_data <-
  with(yellowstone_full,
       expand.grid(
         SpringSnow = seq(
           min(SpringSnow, na.rm = TRUE),
           max(SpringSnow, na.rm =
                 TRUE),
           length = 200
         ),
         AprMax = median(AprMax, na.rm =
                           TRUE),
         SpringTempSum = median(SpringTempSum, na.rm=TRUE),
         MayMin = median(MayMin, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOff)$linkinv
pred_SpringSnow <- predict(mod4_iceOff, new_data, type = "link", se.fit = TRUE)
pred_SpringSnow <- cbind(pred_SpringSnow, new_data)
pred_SpringSnow <- transform(pred_SpringSnow, lwr_ci = ilink(fit - (2 * se.fit)),
                         upr_ci = ilink(fit + (2 * se.fit)),
                         fitted = ilink(fit))
pred_SpringSnow <- pred_SpringSnow %>%
  select(SpringSnow, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_SpringSnow = lwr_ci,
                upr_ci_SpringSnow = upr_ci,
                fitted_SpringSnow = fitted)


breaks_IceOff<-c(210, 225, 240, 255)
# figure out what DOY corresponds to as a date
# c("29-Apr","14-May","29-May","13-Jun")


Fig3G <-
  ggplot(pred_SpringSnow, aes(x = SpringSnow, y = fitted_SpringSnow)) +
  geom_ribbon(aes(ymin = lwr_ci_SpringSnow, ymax = upr_ci_SpringSnow), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=SpringSnow,
                                        y=j_off_wy,
                                        fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative spring snow (mm)",
       y="Ice off")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,
                     labels=c("29-Apr","14-May","29-May","13-Jun"))+
  theme_pubr(border=TRUE, base_size=8)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "lines"),
        axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="g",
                fontface="bold"))


# FIGURE 3  - ice on and off drivers ----------------------------------------------------------------
##Vertically aligned
Fig3A_vert <- Fig3A +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.title.y=element_blank(),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))+
  labs(title="Ice-on")
Fig3B_vert <- Fig3B +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))
Fig3C_vert <- Fig3C +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))

Fig3D_vert <- Fig3D +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.title.y=element_blank(),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))+
  labs(title="Ice-off")
Fig3E_vert <- Fig3E +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))
Fig3F_vert <- Fig3F +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))
Fig3G_vert <- Fig3G +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))

combined <- (Fig3A_vert+ #a
               Fig3D_vert+ #d
               Fig3B_vert+ #b
               Fig3E_vert+ #e
               Fig3C_vert+ #c
               Fig3F_vert+
               patchwork::plot_spacer()+
               Fig3G_vert) & #f
  scale_fill_gradient( low = "white", high = "black",
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              frame.colour = "black",
                                              ticks = TRUE, 
                                              nbin = 10,
                                              label.position = "bottom",
                                              barwidth = 10,
                                              barheight = 1.3, 
                                              title.vjust = 0.75,
                                              direction = 'horizontal'),
                       "Year") &
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=8))

combined <- combined  +
  patchwork::plot_layout(ncol = 2, guides="collect")

combined
# ggsave("Figures/Figure3_GAMS_IceOn_IceOff.pdf", width=6, height=8,units="in", dpi=600)
ggsave("Figures/Figure3_GAMS_IceOn_IceOff.png", width=6, height=8,units="in", dpi=600)




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
