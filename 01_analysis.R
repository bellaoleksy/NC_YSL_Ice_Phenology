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

Timing<-read.csv(file="data/YSL_ice_phenology_20240220.csv",header=T,sep=",")
Weather<-read.csv(file="data/Yellowstone_Snow_Rain.csv",header=T,sep=",")

yellowstone_phenology <- Timing %>%
  mutate(IceOnDate = ymd(CY_IceOn),
         IceOffDate = ymd(CY_IceOff),
         j_on_wy = hydro.day(IceOnDate),
         j_off_wy = hydro.day(IceOffDate),
         ice_days = j_off_wy - j_on_wy)



# Summary statistics -------------------------------------------------------


min(yellowstone_phenology$ice_days, na.rm=TRUE)
max(yellowstone_phenology$ice_days, na.rm=TRUE)
median(yellowstone_phenology$ice_days, na.rm=TRUE)
mean(yellowstone_phenology$ice_days, na.rm=TRUE)


min(yellowstone_phenology$j_on_wy, na.rm=TRUE) # 275 + 62 = december 2
max(yellowstone_phenology$j_on_wy, na.rm=TRUE) # january 28
median(yellowstone_phenology$j_on_wy, na.rm=TRUE) #275 + 87 december 27

min(yellowstone_phenology$j_off_wy, na.rm=TRUE) #DOY 120 April 29th
max(yellowstone_phenology$j_off_wy, na.rm=TRUE) #DOY 165 June 13th 
median(yellowstone_phenology$j_off_wy, na.rm=TRUE) #DOY 145 May 24th 

#How many observations of ice on?
colSums(!is.na(yellowstone_phenology))
#IceOn = 74
#IceOff = 84
# Perc IceOn missing is
(96-74)/96
# Perc IceOff missing is
(96-84)/96

#Annual weather by water-year

#Custom function to make sure sum() doesn't return 0 for column of NAs
sum2 <- function(x) {
  if(all(is.na(x))){
    c(x[0],NA)} else {
      sum(x,na.rm = TRUE)}
}

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



# Correlations with climate variables -----------------------------------------------------------

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
  filter(!column %in% c("FallMax", #corr with OctTempSum, FallTempSum
                        "NovTempSum")) %>% #corr w/ FallTempSum
  arrange(p) %>%
  pull(column)


IceOffVars_full <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p <= 0.05) %>% 
  arrange(row) %>%
  filter(row=="j_off_wy") %>%
  filter(!grepl('Jun|Nov|ice_days', column)) %>%  #spurious likely relationships
  arrange(p) %>%
  pull(column)

#We need to narrow down the IceOffVars. How many are collinear? 
YSL_collinear_trim %>% 
  filter(row %in% IceOffVars_full) %>%
  arrange(row)
IceOffVars_full

#Trim IceOffVars

IceOffVars <- YSL_ice_correlations_trim %>%
  group_by(row) %>%
  filter(p <= 0.05) %>% 
  arrange(row) %>%
  filter(row=="j_off_wy") %>%
  filter(!grepl('Jun|Nov|ice_days', column)) %>%  #spurious likely relationships
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
IceDaysVars_full



# ~ COMPARATIVE LAKES ANALYSIS  ----------------------------------------------------


## Non-Yellowstone lakes
non_ysl <- read_csv(here::here("data/other_phenology.txt"))

non_ysl <- non_ysl %>%
  # select necessary columns
  select(lake, lakecode, start_year, iceOn, iceOff, orig_duration) %>% 
  # filter out lakes: 
  filter(lakecode == "JK25"| # Lake Haukivesi
           lakecode == "JK02"| # Lake Kallavesi
           lakecode == "GW369"| # Lake Kallsjön
           lakecode == "JK03"| # Näsijärvi
           lakecode == "JK05"| # Päijänne
           lakecode == "JK40"| # Pielinen
           lakecode == "NG1" # Lake Baikal
  ) %>%
  # below is jp original
  # caluclate on and off julian dates for water year
  mutate(j_on_wy = hydro.day(iceOn),
         j_off_wy = hydro.day(iceOff)) %>%
  # filter out data only from 1927 to 2022
  filter(start_year >= 1927, start_year <=2022) %>%
  mutate(ice_days = j_off_wy - j_on_wy)

#modify YSL dataset for binding
ysl_bind <- yellowstone_phenology %>%
  mutate(start_year = water_year-1,
         lake = "yellowstone") %>%
  select(-contains(c("_MD","CY_")), -water_year) %>%
  rename(iceOn=IceOnDate,
         iceOff=IceOffDate)

str(ysl_bind)
str(non_ysl)

# combine data sets
full_data <- bind_rows(non_ysl, ysl_bind) %>%
  arrange(lake, start_year)


#=.. GAMs all lakes ----------------------------------------------------------

# for-loop fitting models
lake_names <- full_data %>%
  pull(lake) %>%
  unique()

# ice-on fits
ice_on_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, j_on_wy~s(start_year, k=3), lake_name)
  ice_on_list[[lake]] <- response_mod
}
names(ice_on_list) <- lake_names

# ice-off fits
ice_off_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, j_off_wy~s(start_year, k=3), lake_name)
  ice_off_list[[lake]] <- response_mod
}
names(ice_off_list) <- lake_names

# ice duration fits
ice_days_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, ice_days~s(start_year, k=3), lake_name)
  ice_days_list[[lake]] <- response_mod
}
names(ice_days_list) <- lake_names

# save list of response lists
ice_response_list = list(duration = ice_days_list,
                         off_date = ice_off_list,
                         on_date = ice_on_list)

# gam summaries
# ice-on summary
on_stats <- gam_summary(full_data, j_on_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_on") 
on_stats

# ice-off summary
off_stats <- gam_summary(full_data, j_off_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_off") 
off_stats
# ice duration summary
duration_stats <- gam_summary(full_data, ice_days~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_duration") 
duration_stats

# >>> write csv of gam summaries ####
bind_rows(on_stats, off_stats, duration_stats) %>%
  write_csv("Figures/MS/gam_stat_table.csv")




# Figure 2 - all lakes ----------------------------------------------------

hi_res <- full_data %>% 
  select(lake, start_year, ice_days, j_on_wy, j_off_wy) %>%
  pivot_longer(ice_days:j_off_wy) %>%
  mutate(name = factor(name,
                       levels=c("j_on_wy",
                                "j_off_wy",
                                "ice_days"))) %>%
  ggplot(aes(x = start_year,
             y = value)) +
  geom_point(shape=21, fill="grey50", alpha=0.9) +
  theme_bw() +
  facet_grid(name ~ lake, scales = "free_y") +
  geom_smooth(method = "gam", color="black", formula = y ~s(x, k=3)) + 
  labs(y = "Day of water year") +
  scale_x_continuous(breaks=c(1930,1960,1990,2020),
                     limits=c(1930,2020))+
  theme_pubr(border=TRUE, base_size=8)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
hi_res

ggsave(filename = here::here("Figures/MS/Figure2_full_fig_hi-res.pdf"),
       plot = hi_res,
       dpi = 600)


# ~ MODELS - ICE ON ----------------------------------------------------------------

### I added Family Gamma here since observations are always > 0
mod0_iceOn <- gam(j_on_wy ~ s(water_year),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod0_iceOn)
draw(mod0_iceOn)
appraise(mod0_iceOn)

### Start with all fall variables that were highly correlated
IceOnVars

mod1_iceOn <- gam(j_on_wy ~ s(FallTempSum) + s(OctTempSum) +
                        s(OctMin) + s(FallMin) + s(DecMin)  + s(DecTempSum),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod1_iceOn)
draw(mod1_iceOn)
appraise(mod1_iceOn)

mod2_iceOn <- gam(j_on_wy ~  s(FallTempSum) + 
                    s(OctMin) + s(FallMin) + s(DecMin)  + s(DecTempSum),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod2_iceOn)
draw(mod2_iceOn)
appraise(mod2_iceOn)


mod3_iceOn <- gam(j_on_wy ~   s(FallTempSum) + 
                    s(OctMin) + s(FallMin) + s(DecMin),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod3_iceOn)
draw(mod3_iceOn)
appraise(mod3_iceOn)

mod4_iceOn <- gam(j_on_wy ~  s(FallTempSum) + 
                     s(FallMin) + s(DecMin),
                      family=Gamma(link="log"),
                      data = yellowstone_full,
                      method = "REML")
summary(mod4_iceOn)
draw(mod4_iceOn)
appraise(mod4_iceOn)


mod5_iceOn <- gam(j_on_wy ~   s(FallTempSum) + 
                     s(DecMin),
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

# >>>> Panel A -- Ice On vs. Cumul. Fall Temp--------------------------------------

new_data <-
  with(yellowstone_full,
       expand.grid(
         FallTempSum = seq(
           min(FallTempSum, na.rm = TRUE),
           max(FallTempSum, na.rm =
                 TRUE),
           length = 200
         ),
         DecMin = median(DecMin, na.rm =
                             TRUE),
         FallMin = median(FallMin, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOn)$linkinv
pred_FallTempSum <- predict(mod4_iceOn, new_data, type = "link", se.fit = TRUE)
pred_FallTempSum <- cbind(pred_FallTempSum, new_data)
pred_FallTempSum <- transform(pred_FallTempSum, lwr_ci = ilink(fit - (2 * se.fit)),
                          upr_ci = ilink(fit + (2 * se.fit)),
                          fitted = ilink(fit))

pred_FallTempSum <- pred_FallTempSum %>%
  select(FallTempSum, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_FallTempSum = lwr_ci,
                upr_ci_FallTempSum = upr_ci,
                fitted_FallTempSum = fitted)

#Modify axis labels: fed DOY -> dates
breaks_IceOn<-c(60, 80, 100, 120)
# c("30-Nov","20-Dec","10-Jan","30-Jan") #corresponding cal dates


Fig3A <-
  ggplot(pred_FallTempSum, aes(x = FallTempSum, y = fitted_FallTempSum)) +
  geom_ribbon(aes(ymin = lwr_ci_FallTempSum, ymax = upr_ci_FallTempSum), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=FallTempSum,
                             y=j_on_wy,
                             fill=water_year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative Fall Temperature (°C)",
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
Fig3A

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
         FallTempSum = median(FallTempSum, na.rm =
                            TRUE),
         FallMin = median(FallMin, na.rm=TRUE)
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
Fig3B




# # >>>>  Panel C -- Ice On vs. Cumul Fall Min Temp -------------------------------------
summary(mod4_iceOn)

new_data <-
  with(yellowstone_full,
       expand.grid(
         FallMin = seq(
           min(FallMin, na.rm = TRUE),
           max(FallMin, na.rm =
                 TRUE),
           length = 200
         ),
         FallTempSum = median(FallTempSum, na.rm =
                             TRUE),
         DecMin = median(DecMin, na.rm=TRUE)
       ))

ilink <- family(mod4_iceOn)$linkinv
pred_FallMin <- predict(mod4_iceOn, new_data, type = "link", se.fit = TRUE)
pred_FallMin <- cbind(pred_FallMin, new_data)
pred_FallMin <- transform(pred_FallMin, lwr_ci = ilink(fit - (2 * se.fit)),
                           upr_ci = ilink(fit + (2 * se.fit)),
                           fitted = ilink(fit))
pred_FallMin <- pred_FallMin %>%
  select(FallMin, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_FallMin = lwr_ci,
                upr_ci_FallMin = upr_ci,
                fitted_FallMin = fitted)

breaks_IceOn<-c(60, 80, 100, 120)
# c("30-Nov","20-Dec","10-Jan","30-Jan") #corresponding cal dates


Fig3C <-
  ggplot(pred_FallMin, aes(x = FallMin, y = fitted_FallMin)) +
  geom_ribbon(aes(ymin = lwr_ci_FallMin, ymax = upr_ci_FallMin), alpha = 0.2) +
  geom_line() +
  geom_point(data=yellowstone_full, aes(x=FallMin,
                             y=j_on_wy,
                             fill=water_year), 
             shape=21,alpha=0.9)+
  labs(x="Cumulative min. fall temperatures (°C)",
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
Fig3C

# Export Table Ice On Models ----------------------------------------------

models <- list(
  mod1_iceOn,
  mod2_iceOn,
  mod3_iceOn,
  mod4_iceOn,
  mod5_iceOn
)

compile_gam_outputs <- function(models) {
  map_df(models, ~{
    tidied <- tidy(.x)
    dev_exp <- 1 - .x$deviance / .x$null.deviance
    log_likelihood <- as.numeric(logLik(.x))
    aic <- AIC(.x)
    bind_cols(tidied, dev_exp = dev_exp, log_likelihood = log_likelihood, aic = aic)
  }, .id = "model_name")
}

# Compile GAM model output for Ice On
compiled_outputs_iceon <- compile_gam_outputs(models) %>%
  write_csv("Figures/MS/gam_stat_table_YSL_ice_on.csv")




# ~ MODELS - ICE OFF----------------------------------------------------------------

mod0_iceOff <- gam(j_off_wy ~ s(water_year),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceOff)
draw(mod0_iceOff)
appraise(mod0_iceOff)

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

gam.check(mod1_iceOff)

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
mod4_iceOff <- gam(j_off_wy ~  s(AprMax) + 
                     s(SpringTempSum) +
                     s(MayMin) + s(SpringSnow) ,
                       family=Gamma(link="log"),
                       data = yellowstone_full,
                       method = "REML")
summary(mod4_iceOff)
draw(mod4_iceOff)
appraise(mod4_iceOff)

gam.check(mod4_iceOff)

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


# Figure 3 - ice on and off drivers -----------------------------------------------------------

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
Fig3E


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
Fig3F

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
Fig3G

# Export Table Ice Off Models ----------------------------------------------

models <- list(
  mod1_iceOff,
  mod2_iceOff,
  mod3_iceOff,
  mod4_iceOff,
  mod5_iceOff
)

compile_gam_outputs <- function(models) {
  map_df(models, ~{
    tidied <- tidy(.x)
    dev_exp <- 1 - .x$deviance / .x$null.deviance
    log_likelihood <- as.numeric(logLik(.x))
    aic <- AIC(.x)
    bind_cols(tidied, dev_exp = dev_exp, log_likelihood = log_likelihood, aic = aic)
  }, .id = "model_name")
}

# Compile GAM model output for Ice On
compiled_outputs_iceon <- compile_gam_outputs(models) %>%
  write_csv("Figures/MS/gam_stat_table_YSL_ice_off.csv")


# ........>> COMPILE FIGURE 3 <<........ ----------------------------------------------------------------
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
               patchwork::guide_area()+ #place legend in plot space
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
# ggsave("Figures/MS/Figure3_GAMS_IceOn_IceOff.pdf", width=6, height=8,units="in", dpi=600)
ggsave("Figures/MS/Figure3_GAMS_IceOn_IceOff.png", width=6, height=8,units="in", dpi=600)


# FIGURE 4  - trends in drivers ----------------------------------------------------------------


# >>>> Panel A - Cumul Fall Temp. -------------------------------------------------


mod0_FallTempSum <- gam(FallTempSum ~ s(water_year) ,
                      # family=Gamma(link="log"),
                      data = yellowstone_full,
                      # correlation = corCAR1(form = ~ water_year),
                      method = "REML")
summary(mod0_FallTempSum)
draw(mod0_FallTempSum)
report(mod0_FallTempSum)
acf(residuals(mod0_FallTempSum), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
FallTempSumPred <- cbind(years,
                       data.frame(predict(
                         mod0_FallTempSum, years,
                         type = "response",
                         se.fit = TRUE
                       )))

### Calculate upper and lower bounds
FallTempSumPred <- transform(FallTempSumPred,
                           upper = fit + (2 * se.fit),
                           lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_FallTempSum) #in theory gratia::derivatives should work here
m1.d.2 <- gratia::derivatives(mod0_FallTempSum)
m1.d.3 <- gratia::fderiv(mod0_FallTempSum)

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")


#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(FallTempSumPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Positive rate of change in cumulative snow both the early part of the record and late part of the record

#Add a column for periods of time when the trend is accelerating
FallTempSumPred <- cbind(FallTempSumPred, data.frame(incr=unlist(m1.dsig$incr)),
                       data.frame(decr=unlist(m1.dsig$decr)))


Fig4A <- FallTempSumPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="a",
                fontface="bold"))

Fig4A


# >>>> Panel B - Cumul Dec min --------------------------------------------


mod0_DecMin <- gam(DecMin ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_DecMin)
draw(mod0_DecMin)

# acf(residuals(mod0_DecMin), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
DecMinPred <- cbind(years,
                           data.frame(predict(
                             mod0_DecMin, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
DecMinPred <- transform(DecMinPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_DecMin) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(DecMinPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Winter temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
DecMinPred <- cbind(DecMinPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


Fig4B <- DecMinPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=DecMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", linewidth=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", linewidth=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative Dec. minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="b",
                fontface="bold"))

Fig4B


# >>>> Panel C - Cumul Fall Min Temp  --------------------------------------------


mod0_FallMin <- gam(FallMin ~ s(water_year),
                   # family=Gamma(link="log"),
                   data = yellowstone_full,
                   correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                   #specifies the correlation argument of gam
                   method = "REML")
summary(mod0_FallMin)
draw(mod0_FallMin)

# acf(residuals(mod0_FallMin), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
FallMinPred <- cbind(years,
                    data.frame(predict(
                      mod0_FallMin, years,
                      type = "response",
                      se.fit = TRUE
                    )))

### Calculate upper and lower bounds
FallMinPred <- transform(FallMinPred,
                        upper = fit + (2 * se.fit),
                        lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_FallMin) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or Janreasing trends
m1.dsig <- signifD(FallMinPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Winter temps have been getting warmer in the last few Janades!

#Add a column for periods of time when the trend is accelerating
FallMinPred <- cbind(FallMinPred, data.frame(incr=unlist(m1.dsig$incr)),
                    data.frame(decr=unlist(m1.dsig$decr)))


Fig4C <- FallMinPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="c",
                fontface="bold"))

Fig4C


# >>>> Panel D - Cumul Apr max --------------------------------------------


mod0_AprMax <- gam(AprMax ~ s(water_year),
                   # family=Gamma(link="log"),
                   data = yellowstone_full,
                   correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                   #specifies the correlation argument of gam
                   method = "REML")
summary(mod0_AprMax)
draw(mod0_AprMax)

# acf(residuals(mod0_AprMax), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
AprMaxPred <- cbind(years,
                    data.frame(predict(
                      mod0_AprMax, years,
                      type = "response",
                      se.fit = TRUE
                    )))

### Calculate upper and lower bounds
AprMaxPred <- transform(AprMaxPred,
                        upper = fit + (2 * se.fit),
                        lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_AprMax) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or Janreasing trends
m1.dsig <- signifD(AprMaxPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Winter temps have been getting warmer in the last few Janades!

#Add a column for periods of time when the trend is accelerating
AprMaxPred <- cbind(AprMaxPred, data.frame(incr=unlist(m1.dsig$incr)),
                    data.frame(decr=unlist(m1.dsig$decr)))


Fig4D <- AprMaxPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=AprMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative April maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="d",
                fontface="bold"))

Fig4D



# >>>> Panel E - Cumul Spring temp ----------------------------------------

mod0_SpringTempSum <- gam(SpringTempSum ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_SpringTempSum)
draw(mod0_SpringTempSum)

# acf(residuals(mod0_MaxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SumSpringTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_SpringTempSum, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
SumSpringTempPred <- transform(SumSpringTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_SpringTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative
plot.Deriv(m1.d)
## Maximum Spring temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumSpringTempPred <- cbind(SumSpringTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


Fig4E <- SumSpringTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring mean temperatures (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="e",
                fontface="bold"))

Fig4E



# >>>> Panel F - Cumul May min --------------------------------------------


mod0_MayMin <- gam(MayMin ~ s(water_year),
                   # family=Gamma(link="log"),
                   data = yellowstone_full,
                   correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                   #specifies the correlation argument of gam
                   method = "REML")
summary(mod0_MayMin)
draw(mod0_MayMin)

# acf(residuals(mod0_MayMin), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MayMinPred <- cbind(years,
                    data.frame(predict(
                      mod0_MayMin, years,
                      type = "response",
                      se.fit = TRUE
                    )))

### Calculate upper and lower bounds
MayMinPred <- transform(MayMinPred,
                        upper = fit + (2 * se.fit),
                        lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MayMin) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or Mayreasing trends
m1.dsig <- signifD(MayMinPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Winter temps have been getting warmer in the last few Mayades!

#Add a column for periods of time when the trend is accelerating
MayMinPred <- cbind(MayMinPred, data.frame(incr=unlist(m1.dsig$incr)),
                    data.frame(decr=unlist(m1.dsig$decr)))


Fig4F <- MayMinPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=MayMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative May minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="f",
                fontface="bold"))

Fig4F


# >>>> Panel G - Cumul Spring snow ----------------------------------------


mod0_cumulSpringSnow <- gam(SpringSnow ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = yellowstone_full,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulSpringSnow)
draw(mod0_cumulSpringSnow)

# acf(residuals(mod0_cumulSpringSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulSpringSnowPred <- cbind(years,
                             data.frame(predict(
                               mod0_cumulSpringSnow, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
cumulSpringSnowPred <- transform(cumulSpringSnowPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulSpringSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulSpringSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
cumulSpringSnowPred <- cbind(cumulSpringSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


Fig4G <- cumulSpringSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="g",
                fontface="bold"))

Fig4G

# ........>> COMPILE FIGURE 4 <<........ ----------------------------------------------------------------
##Vertically aligned
Fig4A_vert <- Fig4A +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title="Ice-on")
Fig4B_vert <- Fig4B +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
Fig4C_vert <- Fig4C +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"))

Fig4D_vert <- Fig4D +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(title="Ice-off")
Fig4E_vert <- Fig4E +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
Fig4F_vert <- Fig4F +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
Fig4G_vert <- Fig4G +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"))

combined <- (Fig4A_vert+ #a
               Fig4D_vert+ #d
               Fig4B_vert+ #b
               Fig4E_vert+ #e
               Fig4C_vert+ #c
               Fig4F_vert+
               patchwork::plot_spacer()+ #place legend in plot space
               Fig4G_vert) & #f
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=8))

combined <- combined  +
  patchwork::plot_layout(ncol = 2, guides="collect")

combined


left_panel <- (Fig4A_vert+ #a
  Fig4B_vert+ #b
  Fig4C_vert+ #c
  patchwork::plot_spacer()) +
  patchwork::plot_layout(nrow=4) & 
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=8))

right_panel <- (Fig4D_vert+ #a
                 Fig4E_vert+ #b
                 Fig4F_vert+ #c
                Fig4G_vert) +
  patchwork::plot_layout(nrow=4) & 
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=8))

combined <- left_panel | right_panel
combined

ggsave("Figures/MS/Figure4_TrendInCovariates.pdf", width=8, height=10,units="in", dpi=600)
ggsave("Figures/MS/Figure4_TrendInCovariates.png", width=8, height=10,units="in", dpi=600)


# ~ MODELS - ICE DURATION----------------------------------------------------------------

mod0_iceDur <- gam(ice_days ~ s(water_year),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceDur)
draw(mod0_iceDur)
appraise(mod0_iceDur)

IceDaysVars

mod1_iceDur <- gam(ice_days ~ s(FallTempSum) + s(DecTempSum) + s(DecMin) +
                     s(SpringTempSum) + s(WinterMax) + s(MayMin) ,
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod1_iceDur)
draw(mod1_iceDur)
appraise(mod1_iceDur)

mod2_iceDur <- gam(ice_days ~ s(FallTempSum) + s(DecMin) +
                     s(SpringTempSum) + s(WinterMax) + s(MayMin) ,
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod2_iceDur)
draw(mod2_iceDur)
appraise(mod2_iceDur)


mod3_iceDur <- gam(ice_days ~ s(FallTempSum) + s(DecMin) +
                     s(SpringTempSum) + s(WinterMax),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod3_iceDur)
draw(mod3_iceDur)
appraise(mod3_iceDur)


mod4_iceDur <- gam(ice_days ~ s(FallTempSum) + s(DecMin) +
                     s(SpringTempSum),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod4_iceDur)
draw(mod4_iceDur)
appraise(mod4_iceDur)



# Export Table Ice On Models ----------------------------------------------

models <- list(
  mod1_iceDur,
  mod2_iceDur,
  mod3_iceDur,
  mod4_iceDur
)

compile_gam_outputs <- function(models) {
  map_df(models, ~{
    tidied <- tidy(.x)
    dev_exp <- 1 - .x$deviance / .x$null.deviance
    log_likelihood <- as.numeric(logLik(.x))
    aic <- AIC(.x)
    bind_cols(tidied, dev_exp = dev_exp, log_likelihood = log_likelihood, aic = aic)
  }, .id = "model_name")
}

# Compile GAM model output for Ice On
compiled_outputs_iceon <- compile_gam_outputs(models) %>%
  write_csv("Figures/MS/gam_stat_table_YSL_ice_duration.csv")




# * -----------------------------------------------------------------------
# * -----------------------------------------------------------------------
# * -----------------------------------------------------------------------
# * -----------------------------------------------------------------------


# Discussion points -------------------------------------------------------

#Is spring rain correlated with spring snow?
plot(yellowstone_full$SpringRain,
     yellowstone_full$SpringSnow)
cor.test(yellowstone_full$SpringRain,
     yellowstone_full$SpringSnow)


# * -----------------------------------------------------------------------
# * -----------------------------------------------------------------------
# * -----------------------------------------------------------------------
# **SUPPLEMENTAL FIGURES** ------------------------------------------------


# ~ Annual snow ----------------------------------------------------------------


mod0_cumulSnow <- gam(AnnualSnow ~ s(water_year) ,
                      # family=Gamma(link="log"),
                      data = yellowstone_full,
                      # correlation = corCAR1(form = ~ water_year),
                      method = "REML")
summary(mod0_cumulSnow)
draw(mod0_cumulSnow)
# report(mod0_cumulSnow)
# acf(residuals(mod0_cumulSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulSnowPred <- cbind(years,
                       data.frame(predict(
                         mod0_cumulSnow, years,
                         type = "response",
                         se.fit = TRUE
                       )))

### Calculate upper and lower bounds
cumulSnowPred <- transform(cumulSnowPred,
                           upper = fit + (2 * se.fit),
                           lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulSnow) #in theory gratia::derivatives should work here
m1.d.2 <- gratia::derivatives(mod0_cumulSnow)
m1.d.3 <- gratia::fderiv(mod0_cumulSnow)

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")


#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Positive rate of change in cumulative snow both the early part of the record and late part of the record

#Add a column for periods of time when the trend is accelerating
cumulSnowPred <- cbind(cumulSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                       data.frame(decr=unlist(m1.dsig$decr)))


GAMS_AnnualSnow <- cumulSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=AnnualSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative annual snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualSnow

ggsave(plot=last_plot(), "Figures/MS/GAMS_AnnualSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Annual min temps --------------------------------------------------------------


mod0_MinTemp <- gam(AnnualMin ~ s(water_year),
                    data = yellowstone_full,
                    # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                    #specifies the correlation argument of gam
                    method = "REML")
summary(mod0_MinTemp)
draw(mod0_MinTemp)

# acf(residuals(mod0_MinTemp),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MinTempPred <- cbind(years,
                     data.frame(predict(
                       mod0_MinTemp, years,
                       type = "response",
                       se.fit = TRUE
                     )))

### Calculate upper and lower bounds
MinTempPred <- transform(MinTempPred,
                         upper = fit + (2 * se.fit),
                         lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MinTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MinTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend

#Add a column for periods of time when the trend is accelerating
MinTempPred <- cbind(MinTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                     data.frame(decr=unlist(m1.dsig$decr)))


GAMS_AnnualMin <- MinTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=AnnualMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative annual minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualMin

ggsave(plot=last_plot(), "Figures/MS/GAMS_AnnualMin.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Annual max temp. ----------------------------------------------------------------


mod0_maxTemp <- gam(AnnualMax ~ s(water_year),
                    family=Gamma(link="log"),
                    data = yellowstone_full,
                    correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                    #specifies the correlation argument of gam
                    method = "REML")
summary(mod0_maxTemp)
draw(mod0_maxTemp)

# acf(residuals(mod0_maxTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
maxTempPred <- cbind(years,
                     data.frame(predict(
                       mod0_maxTemp, years,
                       type = "response",
                       se.fit = TRUE
                     )))

### Calculate upper and lower bounds
maxTempPred <- transform(maxTempPred,
                         upper = fit + (2 * se.fit),
                         lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_maxTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend

#Add a column for periods of time when the trend is accelerating
maxTempPred <- cbind(maxTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                     data.frame(decr=unlist(m1.dsig$decr)))


GAMS_AnnualMax <- maxTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=AnnualMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative annual maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualMax

ggsave(plot=last_plot(), "Figures/MS/GAMS_AnnualMax.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Annual rain--------------------------------------------------------------

mod0_AnnualRain <- gam(AnnualRain ~ s(water_year),
                       data = yellowstone_full,
                       # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_AnnualRain)
draw(mod0_AnnualRain)

# acf(residuals(mod0_AnnualRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
AnnualRainPred <- cbind(years,
                        data.frame(predict(
                          mod0_AnnualRain, years,
                          type = "response",
                          se.fit = TRUE
                        )))

### Calculate upper and lower bounds
AnnualRainPred <- transform(AnnualRainPred,
                            upper = fit + (2 * se.fit),
                            lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_AnnualRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(AnnualRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
AnnualRainPred <- cbind(AnnualRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))

GAMS_AnnualRain <- AnnualRainPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=AnnualRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative annual rainfall (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualRain

ggsave(plot=last_plot(), "Figures/MS/GAMS_AnnualRain.png",
       dpi=600, width = 6, height = 5, units = 'in')





# ~ -> FigS1 - annual -----------------------------------------------------

GAMS_AnnualMin + GAMS_AnnualMax + GAMS_AnnualRain + GAMS_AnnualSnow +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Annual climatic trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/MS/FigureS1.AnnualTrends.png",
       dpi=600, width = 6, height = 6, units = 'in')




# ~ Winter snow ----------------------------------------------------------------

mod0_cumulWinterSnow <- gam(WinterSnow ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = yellowstone_full,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulWinterSnow)
draw(mod0_cumulWinterSnow)

# acf(residuals(mod0_cumulWinterSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulWinterSnowPred <- cbind(years,
                             data.frame(predict(
                               mod0_cumulWinterSnow, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
cumulWinterSnowPred <- transform(cumulWinterSnowPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulWinterSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulWinterSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
cumulWinterSnowPred <- cbind(cumulWinterSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterSnow <- cumulWinterSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=WinterSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative winter snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterSnow

ggsave(plot=last_plot(), "Figures/MS/GAMS_WinterSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Winter rain--------------------------------------------------------------

mod0_WinterRain <- gam(WinterRain ~ s(water_year),
                       data = yellowstone_full,
                       # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_WinterRain)
draw(mod0_WinterRain)

# acf(residuals(mod0_WinterRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
WinterRainPred <- cbind(years,
                        data.frame(predict(
                          mod0_WinterRain, years,
                          type = "response",
                          se.fit = TRUE
                        )))

### Calculate upper and lower bounds
WinterRainPred <- transform(WinterRainPred,
                            upper = fit + (2 * se.fit),
                            lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_WinterRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(WinterRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend


#Add a column for periods of time when the trend is accelerating
WinterRainPred <- cbind(WinterRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterRain <- WinterRainPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=WinterRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative winter rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterRain

ggsave(plot=last_plot(), "Figures/MS/GAMS_WinterRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Winter min temp. ----------------------------------------------------------

mod0_minWinterTemp <- gam(WinterMin ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minWinterTemp)
draw(mod0_minWinterTemp)

# acf(residuals(mod0_minWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
minWinterTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_minWinterTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
minWinterTempPred <- transform(minWinterTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_minWinterTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(minWinterTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Winter temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
minWinterTempPred <- cbind(minWinterTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterMin <- minWinterTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=WinterMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative winter minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterMin

ggsave(plot=last_plot(), "Figures/GAMS_WinterMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Winter max temp. ----------------------------------------------------------

mod0_MaxWinterTemp <- gam(WinterMax ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_MaxWinterTemp)
draw(mod0_MaxWinterTemp)

# acf(residuals(mod0_MaxWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MaxWinterTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_MaxWinterTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
MaxWinterTempPred <- transform(MaxWinterTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MaxWinterTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MaxWinterTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum Winter temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
MaxWinterTempPred <- cbind(MaxWinterTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterMax <- MaxWinterTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=WinterMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative winter maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterMax

ggsave(plot=last_plot(), "Figures/GAMS_WinterMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Winter cumul. temp. ----------------------------------------------------------

mod0_WinterTempSum <- gam(WinterTempSum ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_WinterTempSum)
draw(mod0_WinterTempSum)

# acf(residuals(mod0_MaxWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                           max(water_year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SumWinterTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_WinterTempSum, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
SumWinterTempPred <- transform(SumWinterTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_WinterTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumWinterTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative
plot.Deriv(m1.d)
## Maximum Winter temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumWinterTempPred <- cbind(SumWinterTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterTempSum <- SumWinterTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=WinterTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative winter mean temperatures (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  scale_y_continuous(breaks=seq(-1500,0,500),
                     limits = c(-1500,0))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterTempSum

ggsave(plot=last_plot(), "Figures/GAMS_WinterTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS2 - Winter -----------------------------------------------------

GAMS_WinterMin + GAMS_WinterMax + GAMS_WinterRain + GAMS_WinterSnow + GAMS_WinterTempSum + patchwork::plot_spacer() +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Winter (Dec-Feb) climatic trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/MS/FigureS2.WinterTrends.png",
       dpi=600, width = 6, height = 6, units = 'in')




# ~ Spring snow ----------------------------------------------------------------

mod0_cumulSpringSnow <- gam(SpringSnow ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = yellowstone_full,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulSpringSnow)
draw(mod0_cumulSpringSnow)

# acf(residuals(mod0_cumulSpringSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulSpringSnowPred <- cbind(years,
                             data.frame(predict(
                               mod0_cumulSpringSnow, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
cumulSpringSnowPred <- transform(cumulSpringSnowPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulSpringSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulSpringSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
cumulSpringSnowPred <- cbind(cumulSpringSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringSnow <- cumulSpringSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringSnow

ggsave(plot=last_plot(), "Figures/GAMS_SpringSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Spring rain--------------------------------------------------------------

mod0_SpringRain <- gam(SpringRain ~ s(water_year),
                       data = yellowstone_full,
                       # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_SpringRain)
draw(mod0_SpringRain)

# acf(residuals(mod0_SpringRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SpringRainPred <- cbind(years,
                        data.frame(predict(
                          mod0_SpringRain, years,
                          type = "response",
                          se.fit = TRUE
                        )))

### Calculate upper and lower bounds
SpringRainPred <- transform(SpringRainPred,
                            upper = fit + (2 * se.fit),
                            lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_SpringRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SpringRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend


#Add a column for periods of time when the trend is accelerating
SpringRainPred <- cbind(SpringRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringRain <- SpringRainPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringRain

ggsave(plot=last_plot(), "Figures/GAMS_SpringRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Spring min temp. ----------------------------------------------------------

mod0_minSpringTemp <- gam(SpringMin ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minSpringTemp)
draw(mod0_minSpringTemp)

# acf(residuals(mod0_minSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
minSpringTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_minSpringTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
minSpringTempPred <- transform(minSpringTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_minSpringTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(minSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Spring temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
minSpringTempPred <- cbind(minSpringTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringMin <- minSpringTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringMin

ggsave(plot=last_plot(), "Figures/GAMS_SpringMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Spring max temp. ----------------------------------------------------------

mod0_MaxSpringTemp <- gam(SpringMax ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_MaxSpringTemp)
draw(mod0_MaxSpringTemp)

# acf(residuals(mod0_MaxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MaxSpringTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_MaxSpringTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
MaxSpringTempPred <- transform(MaxSpringTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MaxSpringTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MaxSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum Spring temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
MaxSpringTempPred <- cbind(MaxSpringTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringMax <- MaxSpringTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringMax

ggsave(plot=last_plot(), "Figures/GAMS_SpringMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Spring cumul. temp. ----------------------------------------------------------

mod0_SpringTempSum <- gam(SpringTempSum ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_SpringTempSum)
draw(mod0_SpringTempSum)

# acf(residuals(mod0_MaxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SumSpringTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_SpringTempSum, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
SumSpringTempPred <- transform(SumSpringTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_SpringTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative
plot.Deriv(m1.d)
## Maximum Spring temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumSpringTempPred <- cbind(SumSpringTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringTempSum <- SumSpringTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SpringTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative spring mean temperatures (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-1500,0,500),
  #                    limits = c(-15))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringTempSum

ggsave(plot=last_plot(), "Figures/GAMS_SpringTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS3 - Spring -----------------------------------------------------
GAMS_SpringMin + GAMS_SpringMax + GAMS_SpringRain + GAMS_SpringSnow + GAMS_SpringTempSum +
  patchwork::plot_spacer() +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Spring (Mar-May) climatic trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/MS/FigureS3.SpringTrends.png",
       dpi=600, width = 6, height = 6, units = 'in')



# ~ Summer snow ----------------------------------------------------------------

mod0_cumulSummerSnow <- gam(SummerSnow ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = yellowstone_full,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulSummerSnow)
draw(mod0_cumulSummerSnow)

# acf(residuals(mod0_cumulSummerSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulSummerSnowPred <- cbind(years,
                             data.frame(predict(
                               mod0_cumulSummerSnow, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
cumulSummerSnowPred <- transform(cumulSummerSnowPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulSummerSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulSummerSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
cumulSummerSnowPred <- cbind(cumulSummerSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerSnow <- cumulSummerSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SummerSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative summer snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerSnow

ggsave(plot=last_plot(), "Figures/GAMS_SummerSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Summer rain--------------------------------------------------------------

mod0_SummerRain <- gam(SummerRain ~ s(water_year),
                       data = yellowstone_full,
                       # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_SummerRain)
draw(mod0_SummerRain)

# acf(residuals(mod0_SummerRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SummerRainPred <- cbind(years,
                        data.frame(predict(
                          mod0_SummerRain, years,
                          type = "response",
                          se.fit = TRUE
                        )))

### Calculate upper and lower bounds
SummerRainPred <- transform(SummerRainPred,
                            upper = fit + (2 * se.fit),
                            lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_SummerRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SummerRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend


#Add a column for periods of time when the trend is accelerating
SummerRainPred <- cbind(SummerRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerRain <- SummerRainPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SummerRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative summer rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerRain

ggsave(plot=last_plot(), "Figures/GAMS_SummerRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Summer min temp. ----------------------------------------------------------

mod0_minSummerTemp <- gam(SummerMin ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minSummerTemp)
draw(mod0_minSummerTemp)

# acf(residuals(mod0_minSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
minSummerTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_minSummerTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
minSummerTempPred <- transform(minSummerTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_minSummerTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(minSummerTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Summer temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
minSummerTempPred <- cbind(minSummerTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerMin <- minSummerTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SummerMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative summer minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerMin

ggsave(plot=last_plot(), "Figures/GAMS_SummerMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Summer max temp. ----------------------------------------------------------

mod0_MaxSummerTemp <- gam(SummerMax ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_MaxSummerTemp)
draw(mod0_MaxSummerTemp)

# acf(residuals(mod0_MaxSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MaxSummerTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_MaxSummerTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
MaxSummerTempPred <- transform(MaxSummerTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MaxSummerTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MaxSummerTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum Summer temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
MaxSummerTempPred <- cbind(MaxSummerTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerMax <- MaxSummerTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SummerMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative summer maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerMax

ggsave(plot=last_plot(), "Figures/GAMS_SummerMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Summer cumul. temp. ----------------------------------------------------------

mod0_SummerTempSum <- gam(SummerTempSum ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_SummerTempSum)
draw(mod0_SummerTempSum)

# acf(residuals(mod0_MaxSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SumSummerTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_SummerTempSum, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
SumSummerTempPred <- transform(SumSummerTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_SummerTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumSummerTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative
plot.Deriv(m1.d)
## Maximum Summer temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumSummerTempPred <- cbind(SumSummerTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerTempSum <- SumSummerTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=SummerTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative summer mean temperatures (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-1500,0,500),
  #                    limits = c(-15))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerTempSum

ggsave(plot=last_plot(), "Figures/GAMS_SummerTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS4 - Summer -----------------------------------------------------

  
  
GAMS_SummerMin + GAMS_SummerMax + GAMS_SummerRain + GAMS_SummerSnow + GAMS_SummerTempSum +
  patchwork::plot_spacer()+
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Summer (Jun-Sep) climatic trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/MS/FigureS4.SummerTrends.png",
       dpi=600, width = 6, height = 6, units = 'in')




# ~ Fall snow ----------------------------------------------------------------

mod0_cumulFallSnow <- gam(FallSnow ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = yellowstone_full,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulFallSnow)
draw(mod0_cumulFallSnow)

# acf(residuals(mod0_cumulFallSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
cumulFallSnowPred <- cbind(years,
                             data.frame(predict(
                               mod0_cumulFallSnow, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
cumulFallSnowPred <- transform(cumulFallSnowPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_cumulFallSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulFallSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
cumulFallSnowPred <- cbind(cumulFallSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallSnow <- cumulFallSnowPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallSnow

ggsave(plot=last_plot(), "Figures/GAMS_FallSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Fall rain--------------------------------------------------------------

mod0_FallRain <- gam(FallRain ~ s(water_year),
                       data = yellowstone_full,
                       # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_FallRain)
draw(mod0_FallRain)

# acf(residuals(mod0_FallRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
FallRainPred <- cbind(years,
                        data.frame(predict(
                          mod0_FallRain, years,
                          type = "response",
                          se.fit = TRUE
                        )))

### Calculate upper and lower bounds
FallRainPred <- transform(FallRainPred,
                            upper = fit + (2 * se.fit),
                            lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_FallRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(FallRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the trend


#Add a column for periods of time when the trend is accelerating
FallRainPred <- cbind(FallRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallRain <- FallRainPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallRain

ggsave(plot=last_plot(), "Figures/GAMS_FallRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Fall min temp. ----------------------------------------------------------

mod0_minFallTemp <- gam(FallMin ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minFallTemp)
draw(mod0_minFallTemp)

# acf(residuals(mod0_minFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
minFallTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_minFallTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
minFallTempPred <- transform(minFallTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_minFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(minFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum Fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
minFallTempPred <- cbind(minFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallMin <- minFallTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall minimum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallMin

ggsave(plot=last_plot(), "Figures/GAMS_FallMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Fall max temp. ----------------------------------------------------------

mod0_MaxFallTemp <- gam(FallMax ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_MaxFallTemp)
draw(mod0_MaxFallTemp)

# acf(residuals(mod0_MaxFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
MaxFallTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_MaxFallTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
MaxFallTempPred <- transform(MaxFallTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_MaxFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MaxFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum Fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
MaxFallTempPred <- cbind(MaxFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallMax <- MaxFallTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall maximum temperature (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallMax

ggsave(plot=last_plot(), "Figures/GAMS_FallMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Fall cumul. temp. ----------------------------------------------------------

mod0_FallTempSum <- gam(FallTempSum ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = yellowstone_full,
                          correlation = corARMA(form = ~ 1 | water_year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_FallTempSum)
draw(mod0_FallTempSum)

# acf(residuals(mod0_MaxFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(yellowstone_full, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                            max(water_year, na.rm=TRUE),
                                                            length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and water_year, on the response scale.
SumFallTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_FallTempSum, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
SumFallTempPred <- transform(SumFallTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_FallTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative
plot.Deriv(m1.d)
## Maximum Fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumFallTempPred <- cbind(SumFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallTempSum <- SumFallTempPred %>%
  ggplot(aes(x=water_year,y=fit))+
  geom_point(data=yellowstone_full, aes(x=water_year, y=FallTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  # geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  # geom_line(aes(x=water_year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_line(aes(x=water_year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  # geom_ribbon(aes(ymin = (lower), ymax = (upper), x = water_year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Year",y="Cumulative fall mean temperatures (°C)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-1500,0,500),
  #                    limits = c(-15))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallTempSum

ggsave(plot=last_plot(), "Figures/GAMS_FallTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS5 - Fall -----------------------------------------------------



GAMS_FallMin + GAMS_FallMax + GAMS_FallRain + GAMS_FallSnow + GAMS_FallTempSum +
  patchwork::plot_spacer()+
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Fall (Oct-Nov) climatic trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/MS/FigureS5.FallTrends.png",
       dpi=600, width = 6, height = 6, units = 'in')


