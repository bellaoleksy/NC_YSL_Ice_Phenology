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

#Look into QDO, PDO, ENSO
source("0_functions.R")
summarize <- dplyr::summarize

# Read in ice phenology data ----------------------------------------------
Timing<-read.csv(file="Data/R/YSL_Ice.csv",header=T,sep=",")
Weather<-read.csv(file="Data/R/Yellowstone_Snow_Rain.csv",header=T,sep=",")

#Summarize weather data
#Annual
AnnualWeather<-ddply(Weather,.(Year),summarise,AnnualMax=max(max.C,na.rm=T),
                     AnnualMin=min(min.C,na.rm=T),AnnualRain=sum(rain.mm,na.rm=T),
                     AnnualSnow=sum(snow.mm,na.rm=T),SnowDepth=max(SnowDepth.mm,na.rm=T))
#Remove last row
AnnualWeather<-AnnualWeather[-nrow(AnnualWeather),]
#inf because all NAs

#Merge dataframes (Ice on)
YSL1<-merge(x=Timing,y=AnnualWeather,by="Year",all.x=T)
#YSL2<-merge(x=YSL1,y=MonthWeatherOn,by="Year",all.x=T)

#Make seasonal dataframe for ice on and off

#####Ice ON!
#Spring
Spring<-Weather[Weather$Month=="Mar"|Weather$Month=="Apr"|Weather$Month=="May",]
Summer<-Weather[Weather$Month=="Jun"|Weather$Month=="Jul"|Weather$Month=="Aug",]
Fall<-Weather[Weather$Month=="Sep"|Weather$Month=="Oct"|Weather$Month=="Nov",]
Winter<-Weather[Weather$Month=="Dec"|Weather$Month=="Jan"|Weather$Month=="Feb",]

Spring1 <-
  ddply(
    Spring,
    .(Year),
    summarise,
    SpringMax = max(max.C, na.rm = T),
    SpringMin = min(min.C, na.rm = T),
    SpringRain = sum(rain.mm, na.rm = T),
    SpringSnow = sum(snow.mm, na.rm = T),
    SpringTempSum = sum(mean.C, na.rm = T),
    SpringSnowDepth = max(SnowDepth.mm, na.rm = T)
  )
Summer1 <-
  ddply(
    Summer,
    .(Year),
    summarise,
    SummerMax = max(max.C, na.rm = T),
    SummerMin = min(min.C, na.rm = T),
    SummerRain = sum(rain.mm, na.rm = T),
    SummerSnow = sum(snow.mm, na.rm = T),
    SummerTempSum = sum(mean.C, na.rm = T)
    # SummerWind = mean(Knots, na.rm = T)
  )
Fall1 <-
  ddply(
    Fall,
    .(Year),
    summarise,
    FallMax = max(max.C, na.rm = T),
    FallMin = min(min.C, na.rm = T),
    FallRain = sum(rain.mm, na.rm = T),
    FallSnow = sum(snow.mm, na.rm = T),
    FallTempSum = sum(mean.C, na.rm = T)
    # FallWind = mean(Knots, na.rm = T)
  )
Winter1 <-
  ddply(
    Winter,
    .(Year),
    summarise,
    WinterMax = max(max.C, na.rm = T),
    WinterMin = min(min.C, na.rm = T),
    WinterRain = sum(rain.mm, na.rm = T),
    WinterSnow = sum(snow.mm, na.rm = T),
    WinterSnowDepth = max(SnowDepth.mm, na.rm = T),
    WinterTempSum = sum(mean.C, na.rm=T)
  )
#Bind seasonal dataframes together
# SeasonWeather1<-cbind(Spring1,Summer1,Fall1,Winter1)
# #Remove the extra year columns
# SeasonWeather<-subset(SeasonWeather1,select=-c(7,14,21))

SeasonWeather <- full_join(Spring1, Summer1) %>%
  full_join(Fall1) %>%
  full_join(Winter1)

YSLon<-merge(x=YSL1,y=SeasonWeather,by="Year",all.x=T)
#Replace Inf with NA
YSLon=do.call(data.frame,lapply(YSLon,function(value) replace(value,is.infinite(value),NA)))


##ICE-OFF was already compiled by Lusha, here:
YSLoff<-read.csv(file="Data/R/YSLoff_Shifted.csv",header=T,sep=",") %>%
  select(-1)

hist(YSLon$IceOffJulian)
hist(YSLon$IceOnJulian)


## Non-Yellowstone lakes
non_ysl <- read_csv(here::here("Data/other_phenology.txt"))

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
ysl_ice <- YSLon %>%
  select(Year, IceOnDate, IceOnJulian, IceOffDate, IceOffJulian) %>%
  mutate(IceOnJulian_new = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                     TRUE ~ IceOnJulian)) %>%
  mutate(IceOn = ymd(parse_date_time(paste(Year, IceOnJulian_new), orders = "yj")),
         IceOff = ymd(parse_date_time(paste(Year, IceOffJulian), orders = "yj")),
         j_on_wy = hydro.day(IceOn),
         j_off_wy = hydro.day(IceOff),
         ice_days = j_off_wy - j_on_wy,
         start_year = Year,
         water_year = calcWaterYear(IceOn),
         lake = "yellowstone") %>% 
  select(-c(IceOnDate:IceOnJulian_new)) 
  

#IAO WORK HERE!! 
ysl_iceon <- YSLon %>%
  select(Year, IceOnDate, IceOnJulian) %>%
  drop_na(IceOnJulian) %>%
  mutate(IceOnJulian_new = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                     TRUE ~ IceOnJulian)) %>%
  mutate(iceOn = ymd(parse_date_time(paste(Year, IceOnJulian_new), orders = "yj")),
         j_on_wy = hydro.day(iceOn),
         water_year = calcWaterYear(iceOn),
         lake = "yellowstone") %>% 
  select(lake, water_year, iceOn) %>%
  pivot_longer(iceOn) %>%
  group_by(water_year) %>%
  slice(which.max(value)) #use the 2nd date if there are duplicates
  # drop_na(water_year)


ysl_iceoff <- YSLoff %>%
  select(Year, IceOffDate, IceOffJulian) %>%
  drop_na(IceOffJulian) %>%
  mutate(iceOff = ymd(parse_date_time(paste(Year, IceOffJulian), orders = "yj")),
         j_off_wy = hydro.day(iceOff),
         water_year = calcWaterYear(iceOff),
         lake = "yellowstone") %>% 
  select(lake, water_year, iceOff) %>%
  pivot_longer(iceOff)

ysl_bind <- bind_rows(ysl_iceon, ysl_iceoff) %>%
  pivot_wider(values_from = value, names_from = name) %>%
  mutate(j_on_wy = hydro.day(iceOn),
         j_off_wy = hydro.day(iceOff),
         ice_days = j_off_wy - j_on_wy,
         # duration_d = difftime(IceOn,IceOff,units="days"),
         start_year = year(iceOn))

str(ysl_bind)
str(non_ysl)

# combine data sets
full_data <- bind_rows(non_ysl, ysl_bind) %>%
  arrange(lake, start_year)
  # mutate(start_y_c = start_year - mean(start_year, na.rm = TRUE))


# GAMs all lakes ----------------------------------------------------------

# for-loop fitting models ####
lake_names <- full_data %>%
  pull(lake) %>%
  unique()

# ice-on fits ####
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

# save list of response lists ####
ice_response_list = list(duration = ice_days_list,
                         off_date = ice_off_list,
                         on_date = ice_on_list)

# gam summaries ####
# ice-on summary ####
on_stats <- gam_summary(full_data, j_on_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_on") 
on_stats

# ice-off summary ####
off_stats <- gam_summary(full_data, j_off_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_off") 
off_stats

# ice duration summary ####
duration_stats <- gam_summary(full_data, ice_days~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_duration") 
duration_stats

# write csv of gam summaries ####
bind_rows(on_stats, off_stats, duration_stats) %>%
  write_csv("Figures/Tables/gam_stat_table.csv")




# ~ ANNUAL snow ----------------------------------------------------------------


mod0_cumulSnow <- gam(AnnualSnow ~ s(Year) ,
                      # family=Gamma(link="log"),
                      data = YSLon,
                      # correlation = corCAR1(form = ~ water_year),
                      method = "REML")
summary(mod0_cumulSnow)
draw(mod0_cumulSnow)
# report(mod0_cumulSnow)
# acf(residuals(mod0_cumulSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_cumulSnow) #in theory gratia::derivatives should work here
m1.d.2 <- gratia::derivatives(mod0_cumulSnow)
m1.d.3 <- gratia::fderiv(mod0_cumulSnow)

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")


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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=AnnualSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative annual snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualSnow

ggsave(plot=last_plot(), "Figures/GAMS_AnnualSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ ANNUAL min temps --------------------------------------------------------------


mod0_MinTemp <- gam(AnnualMin ~ s(Year),
                    data = YSLon,
                    # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                    #specifies the correlation argument of gam
                    method = "REML")
summary(mod0_MinTemp)
draw(mod0_MinTemp)

# acf(residuals(mod0_MinTemp),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_MinTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=AnnualMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Annual minimum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualMin

ggsave(plot=last_plot(), "Figures/GAMS_AnnualMin.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ ANNUAL max temp. ----------------------------------------------------------------


mod0_maxTemp <- gam(AnnualMax ~ s(Year),
                    family=Gamma(link="log"),
                    data = YSLon,
                    correlation = corARMA(form = ~ 1 | Year, p = 1), 
                    #specifies the correlation argument of gam
                    method = "REML")
summary(mod0_maxTemp)
draw(mod0_maxTemp)

# acf(residuals(mod0_maxTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_maxTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=AnnualMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Annual maximum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualMax

ggsave(plot=last_plot(), "Figures/GAMS_AnnualMax.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ ANNUAL rain--------------------------------------------------------------

mod0_AnnualRain <- gam(AnnualRain ~ s(Year),
                       data = YSLon,
                       # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_AnnualRain)
draw(mod0_AnnualRain)

# acf(residuals(mod0_AnnualRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_AnnualRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=AnnualRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Annual rainfall (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_AnnualRain

ggsave(plot=last_plot(), "Figures/GAMS_AnnualRain.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ ANNUAL Max snow depth ----------------------------------------------------------------


mod0_maxSnowDepth <- gam(SnowDepth ~ s(Year),
                         family=Gamma(link="log"),
                         data = YSLon,
                         # correlation = corCAR1(form = ~ Year),
                         method = "REML")
summary(mod0_maxSnowDepth)
draw(mod0_maxSnowDepth)

# acf(residuals(mod0_maxSnowDepth), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
maxSnowDepthPred <- cbind(years,
                          data.frame(predict(
                            mod0_maxSnowDepth, years,
                            type = "response",
                            se.fit = TRUE
                          )))

### Calculate upper and lower bounds
maxSnowDepthPred <- transform(maxSnowDepthPred,
                              upper = fit + (2 * se.fit),
                              lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "Year"
m1.d <- Deriv(mod0_maxSnowDepth) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxSnowDepthPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Positive rate of change in cumulative snow both the early part of the record and late part of the record

#Add a column for periods of time when the trend is accelerating
maxSnowDepthPred <- cbind(maxSnowDepthPred, data.frame(incr=unlist(m1.dsig$incr)),
                          data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SnowDepth <- maxSnowDepthPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SnowDepth),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Maximum snow depth (mm)")+
  coord_cartesian(xlim=c(1925,2025),
                  ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SnowDepth

ggsave(plot=last_plot(), "Figures/GAMS_SnowDepth.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ -> FigS1 - annual -----------------------------------------------------

GAMS_SnowDepth + GAMS_AnnualMin + GAMS_AnnualMax + GAMS_AnnualRain + GAMS_AnnualSnow +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Annual trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/FigureS1.AnnualTrends.png",
       dpi=600, width = 6, height = 5, units = 'in')




# ~ Winter snow ----------------------------------------------------------------

mod0_cumulWinterSnow <- gam(WinterSnow ~ s(Year),
                            # family=Gamma(link="log"),
                            data = YSLon,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulWinterSnow)
draw(mod0_cumulWinterSnow)

# acf(residuals(mod0_cumulWinterSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_cumulWinterSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=WinterSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative Winter snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterSnow

ggsave(plot=last_plot(), "Figures/GAMS_WinterSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Winter snow depth--------------------------------------------------------------

mod0_WinterSnowDepth <- gam(WinterSnow ~ s(Year),
                            data = YSLoff,
                            # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                            #specifies the correlation argument of gam
                            method = "REML")
summary(mod0_WinterSnowDepth)
draw(mod0_WinterSnowDepth)

# acf(residuals(mod0_WinterSnowDepth),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
WinterSnowDepthPred <- cbind(years,
                             data.frame(predict(
                               mod0_WinterSnowDepth, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
WinterSnowDepthPred <- transform(WinterSnowDepthPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "Year"
m1.d <- Deriv(mod0_WinterSnowDepth) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(WinterSnowDepthPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
WinterSnowDepthPred <- cbind(WinterSnowDepthPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_WinterSnowDepth <- WinterSnowDepthPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLoff, aes(x=Year, y=WinterSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Maximum Winter snow depth")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterSnowDepth

ggsave(plot=last_plot(), "Figures/GAMS_WinterSnowDepth.png",
       dpi=600, width = 6, height = 5, units = 'in')





# ~ Winter rain--------------------------------------------------------------

mod0_WinterRain <- gam(WinterRain ~ s(Year),
                       data = YSLon,
                       # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_WinterRain)
draw(mod0_WinterRain)

# acf(residuals(mod0_WinterRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_WinterRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=WinterRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative Winter rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterRain

ggsave(plot=last_plot(), "Figures/GAMS_WinterRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Winter min temp. ----------------------------------------------------------
# Started here since minimum temps probably control ice formation more so than mean or max?

mod0_minWinterTemp <- gam(WinterMin ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minWinterTemp)
draw(mod0_minWinterTemp)

# acf(residuals(mod0_minWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_minWinterTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=WinterMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Winter minimum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterMin

ggsave(plot=last_plot(), "Figures/GAMS_WinterMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Winter max temp. ----------------------------------------------------------

mod0_MaxWinterTemp <- gam(WinterMax ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_MaxWinterTemp)
draw(mod0_MaxWinterTemp)

# acf(residuals(mod0_MaxWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_MaxWinterTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=WinterMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Winter Maximum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterMax

ggsave(plot=last_plot(), "Figures/GAMS_WinterMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Winter cumul. temp. ----------------------------------------------------------

mod0_WinterTempSum <- gam(WinterTempSum ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          correlation = corARMA(form = ~ 1 | Year, p = 1),
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_WinterTempSum)
draw(mod0_WinterTempSum)

# acf(residuals(mod0_MaxWinterTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_WinterTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=WinterTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Winter cumulative temperatures")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-200,400,200))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_WinterTempSum

ggsave(plot=last_plot(), "Figures/GAMS_WinterTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS2 - Winter -----------------------------------------------------

GAMS_WinterSnowDepth + GAMS_WinterMin + GAMS_WinterMax + GAMS_WinterRain + GAMS_WinterSnow +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Winter trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/FigureS2.WinterTrends.png",
       dpi=600, width = 6, height = 5, units = 'in')




# ~ Spring snow ----------------------------------------------------------------

mod0_cumulSpringSnow <- gam(SpringSnow ~ s(Year),
                            # family=Gamma(link="log"),
                            data = YSLoff,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulSpringSnow)
draw(mod0_cumulSpringSnow)

# acf(residuals(mod0_cumulSpringSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_cumulSpringSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(cumulSpringSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## No acceleration in the change

#Add a column for periods of time when the trend is accelerating
cumulSpringSnowPred <- cbind(cumulSpringSnowPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringSnow <- cumulSpringSnowPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLoff, aes(x=Year, y=SpringSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=1, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative spring snow (mm)")+
  coord_cartesian(xlim=c(1925,2025),
                  ylim=c(0,200))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringSnow

ggsave(plot=last_plot(), "Figures/GAMS_SpringSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')





# ~ Spring snow depth--------------------------------------------------------------

mod0_SpringSnowDepth <- gam(SpringSnowDepth ~ s(Year),
                            data = YSLon,
                            # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                            #specifies the correlation argument of gam
                            method = "REML")
summary(mod0_SpringSnowDepth)
draw(mod0_SpringSnowDepth)

# acf(residuals(mod0_SpringSnowDepth),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
SpringSnowDepthPred <- cbind(years,
                             data.frame(predict(
                               mod0_SpringSnowDepth, years,
                               type = "response",
                               se.fit = TRUE
                             )))

### Calculate upper and lower bounds
SpringSnowDepthPred <- transform(SpringSnowDepthPred,
                                 upper = fit + (2 * se.fit),
                                 lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "Year"
m1.d <- Deriv(mod0_SpringSnowDepth) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SpringSnowDepthPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
SpringSnowDepthPred <- cbind(SpringSnowDepthPred, data.frame(incr=unlist(m1.dsig$incr)),
                             data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringSnowDepth <- SpringSnowDepthPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SpringSnowDepth),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Maximum spring snow depth")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringSnowDepth

ggsave(plot=last_plot(), "Figures/GAMS_SpringSnowDepth.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Spring rain--------------------------------------------------------------

mod0_SpringRain <- gam(SpringRain ~ s(Year),
                       data = YSLoff,
                       # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_SpringRain)
draw(mod0_SpringRain)

# acf(residuals(mod0_SpringRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_SpringRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SpringRainPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)

#Add a column for periods of time when the trend is accelerating
SpringRainPred <- cbind(SpringRainPred, data.frame(incr=unlist(m1.dsig$incr)),
                        data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringRain <- SpringRainPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLoff, aes(x=Year, y=SpringRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative spring rain (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringRain

ggsave(plot=last_plot(), "Figures/GAMS_SpringRain.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Spring min temp. ----------------------------------------------------------
# 

mod0_minSpringTemp <- gam(SpringMin ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minSpringTemp)
draw(mod0_minSpringTemp)

# acf(residuals(mod0_minSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_minSpringTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SpringMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Spring minimum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringMin

ggsave(plot=last_plot(), "Figures/GAMS_SpringMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Spring max temp. ----------------------------------------------------------
# 

mod0_maxSpringTemp <- gam(SpringMax ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_maxSpringTemp)
draw(mod0_maxSpringTemp)

# acf(residuals(mod0_maxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
maxSpringTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_maxSpringTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))

### Calculate upper and lower bounds
maxSpringTempPred <- transform(maxSpringTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "Year"
m1.d <- Deriv(mod0_maxSpringTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## maximum Spring temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
maxSpringTempPred <- cbind(maxSpringTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SpringMax <- maxSpringTempPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SpringMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Spring maximum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringMax

ggsave(plot=last_plot(), "Figures/GAMS_SpringMax.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Spring cumul. temp. ----------------------------------------------------------

mod0_SpringTempSum <- gam(SpringTempSum ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_SpringTempSum)
draw(mod0_SpringTempSum)

# acf(residuals(mod0_MaxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_SpringTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SpringTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Spring cumulative temperatures")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-200,400,200))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SpringTempSum

ggsave(plot=last_plot(), "Figures/GAMS_SpringTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ -> FigS3 - Spring -----------------------------------------------------

GAMS_SpringSnowDepth + GAMS_SpringMin + GAMS_SpringMax + GAMS_SpringRain + GAMS_SpringSnow + GAMS_SpringTempSum +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Spring trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/FigureS3.SpringTrends.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Summer min temp. ----------------------------------------------------------
# 

mod0_minSummerTemp <- gam(SummerMin ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_minSummerTemp)
draw(mod0_minSummerTemp)

# acf(residuals(mod0_minSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_minSummerTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SummerMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Summer minimum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerMin

ggsave(plot=last_plot(), "Figures/GAMS_SummerMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Summer max temp. ----------------------------------------------------------
# 

mod0_maxSummerTemp <- gam(SummerMax ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_maxSummerTemp)
draw(mod0_maxSummerTemp)

# acf(residuals(mod0_maxSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
maxSummerTempPred <- cbind(years,
                           data.frame(predict(
                             mod0_maxSummerTemp, years,
                             type = "response",
                             se.fit = TRUE
                           )))


### Calculate upper and lower bounds
maxSummerTempPred <- transform(maxSummerTempPred,
                               upper = fit + (2 * se.fit),
                               lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "Year"
m1.d <- Deriv(mod0_maxSummerTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxSummerTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## maximum Summer temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
maxSummerTempPred <- cbind(maxSummerTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                           data.frame(decr=unlist(m1.dsig$decr)))


GAMS_SummerMax <- maxSummerTempPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SummerMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Summer maximum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerMax

ggsave(plot=last_plot(), "Figures/GAMS_SummerMax.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Summer cumul. temp. ----------------------------------------------------------

mod0_SummerTempSum <- gam(SummerTempSum ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          correlation = corARMA(form = ~ 1 | Year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_SummerTempSum)
draw(mod0_SummerTempSum)

# acf(residuals(mod0_MaxSummerTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_SummerTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SummerTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Summer cumulative temperatures")+
  coord_cartesian(xlim=c(1925,2025))+
  # ylim=c(-200,400))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(-200,400,200))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerTempSum

ggsave(plot=last_plot(), "Figures/GAMS_SummerTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ Summer snow ----------------------------------------------------------------

mod0_cumulSummerSnow <- gam(SummerSnow ~ s(Year),
                            # family=Gamma(link="log"),
                            data = YSLon,
                            # correlation = corCAR1(form = ~ Year),
                            method = "REML")
summary(mod0_cumulSummerSnow)
draw(mod0_cumulSummerSnow)

# acf(residuals(mod0_cumulSummerSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_cumulSummerSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SummerSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative Summer snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerSnow

ggsave(plot=last_plot(), "Figures/GAMS_SummerSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')





# ~ Summer rain--------------------------------------------------------------

mod0_SummerRain <- gam(SummerRain ~ s(Year),
                       data = YSLon,
                       # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                       #specifies the correlation argument of gam
                       method = "REML")
summary(mod0_SummerRain)
draw(mod0_SummerRain)

# acf(residuals(mod0_SummerRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_SummerRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=SummerRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative Summer rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_SummerRain

ggsave(plot=last_plot(), "Figures/GAMS_SummerRain.png",
       dpi=600, width = 6, height = 5, units = 'in')



# ~ -> FigS4 - Summer -----------------------------------------------------

GAMS_SummerMin + GAMS_SummerMax + GAMS_SummerRain + GAMS_SummerSnow + GAMS_SummerTempSum +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Summer trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/FigureS4.SummerTrends.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Fall snow ----------------------------------------------------------------


mod0_cumulFallSnow <- gam(FallSnow ~ s(Year),
                          # family=Gamma(link="log"),
                          data = YSLon,
                          # correlation = corCAR1(form = ~ Year),
                          method = "REML")
summary(mod0_cumulFallSnow)
draw(mod0_cumulFallSnow)

# acf(residuals(mod0_cumulFallSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_cumulFallSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=FallSnow),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative fall snow (mm)")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallSnow

ggsave(plot=last_plot(), "Figures/GAMS_FallSnow.png",
       dpi=600, width = 6, height = 5, units = 'in')





# ~ Fall rain--------------------------------------------------------------

mod0_FallRain <- gam(FallRain ~ s(Year),
                     data = YSLon,
                     # correlation = corARMA(form = ~ 1 | Year, p = 1), 
                     #specifies the correlation argument of gam
                     method = "REML")
summary(mod0_FallRain)
draw(mod0_FallRain)

# acf(residuals(mod0_FallRain),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_FallRain) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

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
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=FallRain),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Cumulative Fall rain (mm)")+
  # coord_cartesian(xlim=c(1925,2025),
  #                 ylim=c(900,2100))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallRain

ggsave(plot=last_plot(), "Figures/GAMS_FallRain.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Fall min temp. ----------------------------------------------------------
# Started here since minimum temps probably control ice formation more so than mean or max?

mod0_minFallTemp <- gam(FallMin ~ s(Year),
                        # family=Gamma(link="log"),
                        data = YSLon,
                        correlation = corARMA(form = ~ 1 | Year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_minFallTemp)
draw(mod0_minFallTemp)

# acf(residuals(mod0_minFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_minFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(minFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Minimum fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
minFallTempPred <- cbind(minFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                         data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallMin <- minFallTempPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=FallMin),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", linewidth=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", linewidth=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Fall minimum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  # scale_y_continuous(breaks=seq(22,30,2))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallMin

ggsave(plot=last_plot(), "Figures/GAMS_FallMin.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ Fall max temp. ----------------------------------------------------------

mod0_MaxFallTemp <- gam(FallMax ~ s(Year),
                        # family=Gamma(link="log"),
                        data = YSLon,
                        correlation = corARMA(form = ~ 1 | Year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_MaxFallTemp)
draw(mod0_MaxFallTemp)

# acf(residuals(mod0_MaxFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_MaxFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(MaxFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
MaxFallTempPred <- cbind(MaxFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                         data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallMax <- MaxFallTempPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=FallMax),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Fall Maximum temperature")+
  coord_cartesian(xlim=c(1925,2025))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  scale_y_continuous(breaks=seq(5,30,5))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallMax

ggsave(plot=last_plot(), "Figures/GAMS_FallMax.png",
       dpi=600, width = 6, height = 5, units = 'in')

# ~ Fall cumul. temp. ----------------------------------------------------------

mod0_FallTempSum <- gam(FallTempSum ~ s(Year),
                        # family=Gamma(link="log"),
                        data = YSLon,
                        correlation = corARMA(form = ~ 1 | Year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_FallTempSum)
draw(mod0_FallTempSum)

# acf(residuals(mod0_MaxFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(YSLon, data.frame(Year = seq(min(Year, na.rm=TRUE),
                                           max(Year, na.rm=TRUE),
                                           length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
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
Term = "Year"
m1.d <- Deriv(mod0_FallTempSum) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "Year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(SumFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Maximum fall temps have been getting warmer in the last few decades!

#Add a column for periods of time when the trend is accelerating
SumFallTempPred <- cbind(SumFallTempPred, data.frame(incr=unlist(m1.dsig$incr)),
                         data.frame(decr=unlist(m1.dsig$decr)))


GAMS_FallTempSum <- SumFallTempPred %>%
  ggplot(aes(x=Year,y=fit))+
  geom_point(data=YSLon, aes(x=Year, y=FallTempSum),
             shape=21,fill="grey50", alpha=0.5)+ #Plot raw data
  geom_line(size=0.5, alpha=0.8)+ #Plot fitted trend
  geom_line(aes(x=Year, y=incr), color="red", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_line(aes(x=Year, y=decr), color="blue", size=2, alpha=0.8)+ #Highlight period of increasing trend
  geom_ribbon(aes(ymin = (lower), ymax = (upper), x = Year), alpha = 0.5, inherit.aes = FALSE)+ #Plot CI around fitted trend
  labs(x="Water year",y="Fall cumulative temperatures")+
  coord_cartesian(xlim=c(1925,2025),
                  ylim=c(-200,400))+
  scale_x_continuous(breaks=seq(1930, 2020, 15))+
  scale_y_continuous(breaks=seq(-200,400,200))+
  theme_pubr(base_size=8, border=TRUE)
GAMS_FallTempSum

ggsave(plot=last_plot(), "Figures/GAMS_FallTempSum.png",
       dpi=600, width = 6, height = 5, units = 'in')


# ~ -> FigS5 - Fall -----------------------------------------------------

GAMS_FallMin + GAMS_FallMax + GAMS_FallRain + GAMS_FallSnow + GAMS_FallTempSum +
  patchwork::plot_annotation(tag_levels = 'a',
                             title="Fall trends")  & scale_x_continuous(breaks=c(1930,1960,1990,2020))

ggsave(plot=last_plot(), "Figures/FigureS5.FallTrends.png",
       dpi=600, width = 6, height = 5, units = 'in')



# Ice Phenology -----------------------------------------------------------


# ~ MODELS - ICE ON ----------------------------------------------------------------

### I added Family Gamma here since observations are always > 0
mod0_iceOn <- gam(j_on_wy ~ s(start_year),
                  family=Gamma(link="log"),
                  data = ysl_ice,
                  correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod0_iceOn)
draw(mod0_iceOn)
appraise(mod0_iceOn)

#PLOT Autocorrelation function of residuals from the additive model with AR(1) errors
ACF <- acf(resid(mod0_iceOn, type = "response"), plot = FALSE)
ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
ggplot(ACF, aes(x = Lag, y = ACF)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = Lag, yend = 0))
#Suggests that an AR(1) model isn't necessary


### I added Family Gamma here since observations are always > 0
mod1_iceOn <- gam(IceOnJulian ~ s(FallMin)+s(FallMax)+s(FallRain)+s(FallSnow)+
                    s(FallTempSum),
                  family=Gamma(link="log"),
                  data = YSLon,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod1_iceOn)
draw(mod1_iceOn)
appraise(mod1_iceOn)


### I added Family Gamma here since observations are always > 0
mod2_iceOn <- gam(IceOnJulian ~ s(FallMin)+s(FallMax)+s(FallSnow)+
                    s(FallTempSum),
                  family=Gamma(link="log"),
                  data = YSLon,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod2_iceOn)
draw(mod2_iceOn)
appraise(mod2_iceOn)

### I added Family Gamma here since observations are always > 0
mod3_iceOn <- gam(IceOnJulian ~ s(FallMin)+s(FallSnow)+
                    s(FallTempSum),
                  family=Gamma(link="log"),
                  data = YSLon,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod3_iceOn)
draw(mod3_iceOn)
appraise(mod3_iceOn)

### I added Family Gamma here since observations are always > 0
mod4_iceOn <- gam(IceOnJulian ~ s(FallMin)+s(FallSnow),
                  family=Gamma(link="log"),
                  data = YSLon,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod4_iceOn)
draw(mod4_iceOn)
appraise(mod4_iceOn)


#How to compare the fits of multiple GAMs models?
#Would be worth digging into more but found this as a solution:
#https://rdrr.io/cran/itsadug/man/compareML.html
#From the documentation: "This method is preferred over other functions such as AIC for models that include an AR1 model or random effects (especially nonlinear random smooths using bs='fs'). CompareML also reports the AIC difference, but that value should be treated with care."

compareML(mod3_iceOn, mod4_iceOn) #Very similar, mod4 might slightly better than 3
#but deviance explained much higher in model 3

hist(mod3_iceOn$residuals)


# Ice On Figure -----------------------------------------------------------
#annotate panel letters inside plot
panelLetter.normal <- data.frame(
  xpos = c(-Inf),
  ypos =  c(Inf),
  hjustvar = c(-0.5) ,
  vjustvar = c(1.5))

summary(mod3_iceOn)

# ... Panel A -- Ice On vs. Fall Min --------------------------------------

new_data <-
  with(YSLon,
       expand.grid(
         FallMin = seq(
           min(FallMin, na.rm = TRUE),
           max(FallMin, na.rm =
                 TRUE),
           length = 200
         ),
         FallSnow = median(FallSnow, na.rm =
                             TRUE),
         FallTempSum = median(FallTempSum, na.rm=TRUE)
       ))

ilink <- family(mod3_iceOn)$linkinv
pred_FallMin <- predict(mod3_iceOn, new_data, type = "link", se.fit = TRUE)
pred_FallMin <- cbind(pred_FallMin, new_data)
pred_FallMin <- transform(pred_FallMin, lwr_ci = ilink(fit - (2 * se.fit)),
                          upr_ci = ilink(fit + (2 * se.fit)),
                          fitted = ilink(fit))

pred_FallMin <- pred_FallMin %>%
  select(FallMin, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_FallMin = lwr_ci,
                upr_ci_FallMin = upr_ci,
                fitted_FallMin = fitted)

#Modify axis labels: fed DOY -> dates
# labels_IceOnDayofYear_fed<-c(50,70,90,110,130)
# labels_IsoMax_fed<-c(70,90,110)
#To figure out conversion of FED DOY to Date use something like: as.Date(274+110, origin="2014-01-02")
breaks_IceOn<-c(340,350,360,370,380,390)
# figure out what DOY corresponds to as a date
# as.Date(390-365, origin = "2016-01-01")


IceOn_FallMin<-
  ggplot(pred_FallMin, aes(x = FallMin, y = fitted_FallMin)) +
  geom_ribbon(aes(ymin = lwr_ci_FallMin, ymax = upr_ci_FallMin), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLon, aes(x=FallMin,
                             y=IceOnJulian,
                             fill=Year), 
             shape=21, alpha=0.9)+
  labs(x="Fall min. temp (°C)",
       # x=expression(Iso["max,"]["17day,"]["0°C"]),
       # x="SpringSnow Formula: TempMax in degC, 17 day window, 0 degC threshold",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,labels=c("06-Dec","16-Dec","26-Dec","06-Jan","16-Jan","26-Jan"))+
  scale_x_continuous(breaks=c(-35,-30,-25,-20,-15),limits=c(-35,-15))+
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


# ... Panel B -- Ice On vs. Fall Snow -------------------------------------

new_data <-
  with(YSLon,
       expand.grid(
         FallSnow = seq(
           min(FallSnow, na.rm = TRUE),
           max(FallSnow, na.rm =
                 TRUE),
           length = 200
         ),
         FallMin = median(FallMin, na.rm =
                            TRUE),
         FallTempSum = median(FallTempSum, na.rm=TRUE)
       ))

ilink <- family(mod3_iceOn)$linkinv
pred_FallSnow <- predict(mod3_iceOn, new_data, type = "link", se.fit = TRUE)
pred_FallSnow <- cbind(pred_FallSnow, new_data)
pred_FallSnow <- transform(pred_FallSnow, lwr_ci = ilink(fit - (2 * se.fit)),
                           upr_ci = ilink(fit + (2 * se.fit)),
                           fitted = ilink(fit))
pred_FallSnow <- pred_FallSnow %>%
  select(FallSnow, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_FallSnow = lwr_ci,
                upr_ci_FallSnow = upr_ci,
                fitted_FallSnow = fitted)

#Modify axis labels: fed DOY -> dates
# labels_IceOnDayofYear_fed<-c(50,70,90,110,130)
# labels_IsoMax_fed<-c(70,90,110)
#To figure out conversion of FED DOY to Date use something like: as.Date(274+110, origin="2014-01-02")


IceOn_FallSnow<-
  ggplot(pred_FallSnow, aes(x = FallSnow, y = fitted_FallSnow)) +
  geom_ribbon(aes(ymin = lwr_ci_FallSnow, ymax = upr_ci_FallSnow), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLon, aes(x=FallSnow,
                             y=IceOnJulian,
                             fill=Year), 
             shape=21,alpha=0.9)+
  labs(x="Fall snow (mm)",
       # x=expression(Iso["max,"]["17day,"]["0°C"]),
       # x="SpringSnow Formula: TempMax in degC, 17 day window, 0 degC threshold",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,labels=c("06-Dec","16-Dec","26-Dec","06-Jan","16-Jan","26-Jan"))+
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
                label="b",
                fontface="bold"))



# ... Panel C -- Ice On vs. Fall Cumul. Temp  -------------------------------------

new_data <-
  with(YSLon,
       expand.grid(
         FallTempSum = seq(
           min(FallTempSum, na.rm = TRUE),
           max(FallTempSum, na.rm =
                 TRUE),
           length = 200
         ),
         FallMin = median(FallMin, na.rm =
                            TRUE),
         FallSnow = median(FallSnow, na.rm=TRUE)
       ))

ilink <- family(mod3_iceOn)$linkinv
pred_FallTempSum <- predict(mod3_iceOn, new_data, type = "link", se.fit = TRUE)
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
# labels_IceOnDayofYear_fed<-c(50,70,90,110,130)
# labels_IsoMax_fed<-c(70,90,110)
#To figure out conversion of FED DOY to Date use something like: as.Date(274+110, origin="2014-01-02")


IceOn_FallTempSum<-
  ggplot(pred_FallTempSum, aes(x = FallTempSum, y = fitted_FallTempSum)) +
  geom_ribbon(aes(ymin = lwr_ci_FallTempSum, ymax = upr_ci_FallTempSum), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLon, aes(x=FallTempSum,
                             y=IceOnJulian,
                             fill=Year), 
             shape=21,alpha=0.9)+
  labs(x="Fall cumulative temperatures (°C)",
       # x=expression(Iso["max,"]["17day,"]["0°C"]),
       # x="SpringSnow Formula: TempMax in degC, 17 day window, 0 degC threshold",
       y="Ice-on date")+
  scale_y_continuous(breaks=breaks_IceOn,labels=c("06-Dec","16-Dec","26-Dec","06-Jan","16-Jan","26-Jan"))+
  # scale_x_continuous(breaks=labels_IsoMax_fed,labels=c("12-Dec","01-Jan","21-Jan"),limits=c(70,120))+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  theme_pubr(border=TRUE, base_size=8)+
  theme(
    # axis.text.y=element_blank(),
    # axis.ticks.y=element_blank(),
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



#...Ice On  predictor trends----------------------------------------


# 
GAMS_FallMin_new <-GAMS_FallMin +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="a",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(title="Ice-on")

GAMS_FallSnow_new <-GAMS_FallSnow +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="b",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


GAMS_FallTempSum_new <-GAMS_FallTempSum +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="c",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"))

GAMS_SummerTempSum_new <- GAMS_SummerTempSum +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="d",
                fontface="bold"))
# 
# 
# IceOn_preds <- (GAMS_FallMin_new+GAMS_FallSnow_new+GAMS_FallTempSum_new)+
#   patchwork::plot_layout(ncol = 3)
# IceOn_preds & scale_x_continuous(breaks=c(1930,1960,1990,2020))
# 
# 
# ggsave("Figures/GAMS_IceOn_Predictors.png", plot=IceOn_preds, width=10, height=3.5,units="in", dpi=300)
# 
# 
# IceOn / IceOn_preds
# 
# ggsave("Figures/GAMS_IceOn_WithPredictors.png", width=10, height=5,units="in", dpi=300)


# ~ MODELS - ICE OFF----------------------------------------------------------------

#Distribution of y
hist(YSLon$IceOffJulian)

### I added Family Gamma here since observations are always > 0
# mod0_iceOff <- gam(IceOffJulian ~ s(Year),
#                    family=Gamma(link="log"),
#                    data = YSLoff,
#                    # correlatiOff = corCAR1(form = ~ Year),
#                    method = "REML")
# summary(mod0_iceOff)
# draw(mod0_iceOff)
# appraise(mod0_iceOff)

mod0_iceOff <- gam(j_off_wy ~ s(start_year),
                   family=Gamma(link="log"),
                   data = ysl_ice,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceOff)
draw(mod0_iceOff)
appraise(mod0_iceOff)

#PLOT AutocorrelatiOff functiOff of residuals from the additive model with AR(1) errors
# ACF <- acf(resid(mod0_iceOff, type = "respOffse"), plot = FALSE)
# ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
# ggplot(ACF, aes(x = Lag, y = ACF)) +
#   geom_hline(aes(yintercept = 0)) +
#   geom_segment(mapping = aes(xend = Lag, yend = 0))
#Suggests that an AR(1) model isn't necessary

## I added Family Gamma here since observatiOffs are always > 0
mod1_iceOff <- gam(IceOffJulian ~ s(WinterSnow)+ s(SnowDepth) + s(SpringSnow)+ s(SpringRain) + s(WinterMin)+ s(SpringMin) + s(SpringMax),
                   family=Gamma(link="log"),
                   data = YSLoff,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod1_iceOff)
draw(mod1_iceOff)
appraise(mod1_iceOff)


### I added Family Gamma here since observatiOffs are always > 0
mod2_iceOff <- gam(IceOffJulian ~  s(WinterSnow)+ s(SnowDepth) + s(SpringRain) +s(SpringSnow)+ s(WinterMin)+ s(SpringMax),
                   family=Gamma(link="log"),
                   data = YSLoff,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod2_iceOff)
draw(mod2_iceOff)
appraise(mod2_iceOff)

### I added Family Gamma here since observatiOffs are always > 0
mod3_iceOff <- gam(IceOffJulian ~  s(SnowDepth) + s(SpringRain) +s(SpringSnow)+ s(WinterMin)+ s(SpringMax),
                   family=Gamma(link="log"),
                   data = YSLoff,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod3_iceOff)
draw(mod3_iceOff)
appraise(mod3_iceOff)

### I added Family Gamma here since observatiOffs are always > 0
mod4_iceOff <- gam(IceOffJulian ~  s(SnowDepth) + s(SpringRain) +s(SpringSnow)+ s(WinterMin),
                   family=Gamma(link="log"),
                   data = YSLoff,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod4_iceOff)
draw(mod4_iceOff)
appraise(mod4_iceOff)

### I added Family Gamma here since observatiOffs are always > 0
mod5_iceOff <- gam(IceOffJulian ~  s(SnowDepth) + s(SpringRain) +s(SpringSnow),
                   family=Gamma(link="log"),
                   data = YSLoff,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod5_iceOff)
draw(mod5_iceOff)
appraise(mod5_iceOff)
# 
# #How to compare the fits of multiple GAMs models? 
# #Would be worth digging into more but found this as a solutiOff:
# #https://rdrr.io/cran/itsadug/man/compareML.html
# #From the documentatiOff: "This method is preferred over other functiOffs such as AIC for models that include an AR1 model or random effects (especially nOfflinear random smooths using bs='fs'). CompareML also reports the AIC difference, but that value should be treated with care."
# 
# compareML(mod0_iceOff, mod1_iceOff) 
# compareML(mod1_iceOff, mod2_iceOff)
# compareML(mod2_iceOff, mod3_iceOff)
# compareML(mod3_iceOff, mod4_iceOff)
# compareML(mod4_iceOff, mod5_iceOff)
# 
# 
# visreg::visreg2d(mod2_iceOff, xvar='SnowDepth', yvar='SpringSnow', scale='response')
# visreg::visreg2d(mod4_iceOff, xvar='SnowDepth', yvar='SpringSnow', scale='response')
# visreg::visreg2d(mod2_iceOff, xvar='SpringSnow', yvar='SpringRain', scale='response')
# 
# 
# #Top 3 models 
# 
# 
# #Effective degrees of freedom (as a metric for model complexity)
# sum(influence(mod2_iceOff))
# sum(influence(mod3_iceOff))
# sum(influence(mod4_iceOff))
# 
# #terms
# mod2_iceOff$terms
# mod3_iceOff$terms
# mod4_iceOff$terms
# 
# #AIC
# mod1_iceOff$aic
# mod2_iceOff$aic
# mod3_iceOff$aic
# mod4_iceOff$aic
# 
# #Dev explained
# summary(mod2_iceOff)$dev.expl
# summary(mod3_iceOff)$dev.expl
# summary(mod4_iceOff)$dev.expl
# summary(mod5_iceOff)$dev.expl
# 


# ---- Ice Off Figure ---- ------------------------------------------------



#Final variables for paper--
summary(mod5_iceOff)




# ...Panel A -- Ice Off vs. SpringSnow ------------------------------------

mean(YSLoff$IceOffJulian, na.rm=TRUE)
as.Date(142, origin = "2016-01-01")
min(YSLoff$IceOffJulian, na.rm=TRUE)
as.Date(118, origin = "2016-01-01")
max(YSLoff$IceOffJulian, na.rm=TRUE)
as.Date(163, origin = "2016-01-01")


new_data <-
  with(YSLoff,
       expand.grid(
         SpringSnow = seq(
           min(SpringSnow, na.rm = TRUE),
           max(SpringSnow, na.rm =
                 TRUE),
           length = 200
         ),
         SpringRain = median(SpringRain, na.rm =
                               TRUE),
         SnowDepth = median(SnowDepth, na.rm=TRUE),
         WinterMin = median(WinterMin, na.rm=TRUE),
         SpringMax = median(SpringMax, na.rm=TRUE)
       ))

ilink <- family(mod3_iceOff)$linkinv
pred_SpringSnow <- predict(mod3_iceOff, new_data, type = "link", se.fit = TRUE)
pred_SpringSnow <- cbind(pred_SpringSnow, new_data)
pred_SpringSnow <- transform(pred_SpringSnow, lwr_ci = ilink(fit - (2 * se.fit)),
                             upr_ci = ilink(fit + (2 * se.fit)),
                             fitted = ilink(fit))
pred_SpringSnow <- pred_SpringSnow %>%
  select(SpringSnow, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_SpringSnow = lwr_ci,
                upr_ci_SpringSnow = upr_ci,
                fitted_SpringSnow = fitted)

#Modify axis labels: 
breaks_IceOff<-c(120,130,140,150,160)
# figure out what DOY corresponds to as a date
as.Date(160, origin = "2016-01-01")


IceOff_SpringSnow<-
  ggplot(pred_SpringSnow, aes(x = SpringSnow, y = fitted_SpringSnow)) +
  geom_ribbon(aes(ymin = lwr_ci_SpringSnow, ymax = upr_ci_SpringSnow), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLoff, aes(x=SpringSnow,
                              y=IceOffJulian,
                              fill=Year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative spring snow (mm)",
       y="Ice-off date")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,labels=c("30-Apr","10-May","20-May","30-May","09-Jun"),limits=c(120,160))+
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


# ...Panel B -- Ice Off vs. Spring Rain -----------------------------------


new_data <-
  with(YSLoff,
       expand.grid(
         SpringRain = seq(
           min(SpringRain, na.rm = TRUE),
           max(SpringRain, na.rm =
                 TRUE),
           length = 200
         ),
         SpringSnow = median(SpringSnow, na.rm =
                               TRUE),
         SnowDepth = median(SnowDepth, na.rm=TRUE),
         WinterMin = median(WinterMin, na.rm=TRUE),
         SpringMax = median(SpringMax, na.rm=TRUE)
       ))


ilink <- family(mod3_iceOff)$linkinv
pred_SpringRain <- predict(mod3_iceOff, new_data, type = "link", se.fit = TRUE)
pred_SpringRain <- cbind(pred_SpringRain, new_data)
pred_SpringRain <- transform(pred_SpringRain, lwr_ci = ilink(fit - (2 * se.fit)),
                             upr_ci = ilink(fit + (2 * se.fit)),
                             fitted = ilink(fit))
pred_SpringRain <- pred_SpringRain %>%
  select(SpringRain, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_SpringRain = lwr_ci,
                upr_ci_SpringRain = upr_ci,
                fitted_SpringRain = fitted)


IceOff_SpringRain<-
  ggplot(pred_SpringRain, aes(x = SpringRain, y = fitted_SpringRain)) +
  geom_ribbon(aes(ymin = lwr_ci_SpringRain, ymax = upr_ci_SpringRain), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLoff, aes(x=SpringRain,
                              y=IceOffJulian,
                              fill=Year), 
             shape=21, alpha=0.9)+
  labs(x="Cumulative spring rain (mm)",
       y="Ice off")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,labels=c("30-Apr","10-May","20-May","30-May","09-Jun"),limits=c(120,160))+
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



# ...Panel C - Ice off versus Snow Depth ----------------------------------
new_data <-
  with(YSLoff,
       expand.grid(
         SnowDepth = seq(
           min(SnowDepth, na.rm = TRUE),
           max(SnowDepth, na.rm =
                 TRUE),
           length = 200
         ),
         SpringSnow = median(SpringSnow, na.rm =
                               TRUE),
         SpringRain = median(SpringRain, na.rm=TRUE),
         WinterMin = median(WinterMin, na.rm=TRUE),
         SpringMax = median(SpringMax, na.rm=TRUE)
       ))

ilink <- family(mod3_iceOff)$linkinv
pred_SnowDepth <- predict(mod3_iceOff, new_data, type = "link", se.fit = TRUE)
pred_SnowDepth <- cbind(pred_SnowDepth, new_data)
pred_SnowDepth <- transform(pred_SnowDepth, lwr_ci = ilink(fit - (2 * se.fit)),
                            upr_ci = ilink(fit + (2 * se.fit)),
                            fitted = ilink(fit))
pred_SnowDepth <- pred_SnowDepth %>%
  select(SnowDepth, lwr_ci:fitted) %>%
  dplyr::rename(lwr_ci_SnowDepth = lwr_ci,
                upr_ci_SnowDepth = upr_ci,
                fitted_SnowDepth = fitted)


IceOff_SnowDepth<-
  ggplot(pred_SnowDepth, aes(x = SnowDepth, y = fitted_SnowDepth)) +
  geom_ribbon(aes(ymin = lwr_ci_SnowDepth, ymax = upr_ci_SnowDepth), alpha = 0.2) +
  geom_line() +
  geom_point(data=YSLoff, aes(x=SnowDepth,
                              y=IceOffJulian,
                              fill=Year), 
             shape=21, alpha=0.9)+
  labs(x="Maximum snow depth (mm)",
       y="Ice off")+
  grafify::scale_fill_grafify(palette = "grey_conti", name = "Year")+ #yellow_conti scheme
  scale_y_continuous(breaks=breaks_IceOff,labels=c("30-Apr","10-May","20-May","30-May","09-Jun"),limits=c(120,160))+
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


#...Composite Ice Off Figure----------------------------------------

# 
# 
# IceOff<-(IceOff_SpringSnow+IceOff_SpringRain+IceOff_SnowDepth) + 
#   patchwork::plot_layout(ncol = 3, guides="collect") &
#   theme(legend.position="top")
# IceOff
# 
# ggsave("Figures/GAMS_IceOff.png", plot=IceOff, width=8, height=4,units="in", dpi=300)
# 

# ...Compile predictors of ice-off TS -------------------------------------

# 
GAMS_SpringSnow_new <-GAMS_SpringSnow +
  # theme(plot.margin=unit(c(0.5,0,0.5,0.5), "lines"),
  #       axis.text.y = element_text(angle = 90, hjust=0.5)) +
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="d",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(title="Ice-off")


GAMS_SpringRain_new <- GAMS_SpringRain +
  # theme(axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       plot.margin=unit(c(0.5,0,0.5,0), "lines"),
  #       axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="e",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

GAMS_SnowDepth_new <- GAMS_SnowDepth +
  # theme(axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank(),
  #       axis.title.y=element_blank(),
  #       plot.margin=unit(c(0.5,0,0.5,0), "lines"),
  #       axis.ticks.length.y = unit(0, "pt"))+
  geom_text(data=panelLetter.normal,
            aes(x=xpos,
                y=ypos,
                hjust=hjustvar,
                vjust=vjustvar,
                label="f",
                fontface="bold"))+
  theme(plot.margin=unit(c(0,0,0,0), "lines"))
# 
# #...Ice Off  predictor trends----------------------------------------
# 
# IceOff_preds <- (GAMS_SpringSnow_new+GAMS_SpringRain_new+GAMS_SnowDepth_new+patchwork::plot_spacer()) +  patchwork::plot_layout(ncol = 4)
# IceOff_preds
# 
# 
# ggsave("Figures/GAMS_IceOff_Predictors.png", plot=IceOff_preds, width=8, height=3,units="in", dpi=300)
# 

# MODELS - ICE DURATION ---------------------------------------------------

# YSLoff2 <- YSLoff %>%
#   mutate(IceDuration = IceOnJulian-IceOffJulian)
# 
# 
# ### I added Family Gamma here since observatiOffs are always > 0
# mod0_iceDuration <- gam(IceDuration ~ s(Year),
#                         family=Gamma(link="log"),
#                         data = YSLoff2,
#                         # correlatiOff = corCAR1(form = ~ Year),
#                         method = "REML")
# summary(mod0_iceDuration)
# draw(mod0_iceDuration)
# appraise(mod0_iceDuration)


mod0_iceDuration <- gam(ice_days ~ s(start_year),
                        family=Gamma(link="log"),
                        data = ysl_ice,
                        # correlatiOff = corCAR1(form = ~ Year),
                        method = "REML")
summary(mod0_iceDuration)
draw(mod0_iceDuration)
appraise(mod0_iceDuration)

#PLOT AutocorrelatiOff function of residuals from the additive model with AR(1) errors
ACF <- acf(resid(mod0_iceDuration, type = "response"), plot = FALSE)
ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
ggplot(ACF, aes(x = Lag, y = ACF)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = Lag, yend = 0))
#Suggests that an AR(1) model isn't necessary

min(ysl_ice$ice_days, na.rm=TRUE)
max(ysl_ice$ice_days, na.rm=TRUE)
mean(ysl_ice$ice_days, na.rm=TRUE)

# FIGURE 2  - TS Ice On, Off, Duration ----------------------------------------------------------------


hi_res <- full_data %>% 
  select(lake, start_year, ice_days, j_on_wy, j_off_wy) %>%
  pivot_longer(ice_days:j_off_wy) %>%
  mutate(name = factor(name,
                       levels=c("j_on_wy",
                                "j_off_wy",
                                "ice_days"),
                       labels=c("Ice on (days since Oct 1)",
                                "Ice off (days since Oct 1)",
                                "Ice duration (days)"))) %>%
  ggplot(aes(x = start_year,
             y = value)) +
  geom_point(shape=21, fill="grey50", size=3, alpha=0.9) +
  theme_bw() +
  facet_grid(name ~ lake, scales = "free_y", switch="y") +
  geom_smooth(method = "gam", color="black", formula = y ~s(x, k=3)) + 
  scale_x_continuous(breaks=c(1930,1960,1990,2020),
                     limits=c(1930,2020))+
  labs(x="Year")+
  theme_pubr(border=TRUE, base_size=16)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold"))
hi_res

#wy_day 150 = Feb 28
#wy_day 100 = Jan 10
#wy_day 50 = Nov 20
#wy_day 240 = May 29
#wy_day 220 = May 9
#wy_day 200 = Apr 19


ggsave(filename = here::here("Figures/Figure2_full_fig_hi-res.pdf"),
       plot = hi_res,
       dpi = 600)
ggsave(filename = here::here("Figures/Figure2_full_fig_hi-res.png"),
       plot = hi_res,
       dpi = 600)

# #How many complete ice-on observations?
# length((YSLon %>%
#           drop_na(IceOnJulian) %>%
#           pull(IceOnJulian)))
# #How many complete ice-off obsersvations?
# length((YSLoff %>%
#           drop_na(IceOffJulian) %>%
#           pull(IceOffJulian)))
# 
# 
# IceOnTS <-
#   ggplot() + 
#   geom_point(data=YSLon, aes(x=Year,
#                              y=IceOnJulian), 
#              shape=21, fill="grey50", alpha=0.9)+
#   labs(x="Year",
#        y="Ice-on date")+
#   # scale_y_continuous(breaks=breaks_IceOn,labels=c("06-Dec","16-Dec","26-Dec","06-Jan","16-Jan","26-Jan"))+
#   coord_cartesian(xlim=c(1925,2025))+
#   scale_x_continuous(breaks=seq(1930, 2020, 15))+
#   theme_pubr(border=TRUE, base_size=8)+
#   theme(plot.margin=unit(c(0.5,0.5,0,0.5), "lines"),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.x=element_blank())+
#   # axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2)) +
#   geom_text(data=panelLetter.normal,
#             aes(x=xpos,
#                 y=ypos,
#                 hjust=hjustvar,
#                 vjust=vjustvar,
#                 label="a",
#                 fontface="bold"))
#   
# IceOffTS <-
#   ggplot() + 
#   geom_point(data=YSLoff, aes(x=Year,
#                               y=IceOffJulian), 
#              shape=21, fill="grey50", alpha=0.9)+
#   labs(x="Year",
#        y="Ice-off date")+
#   scale_y_continuous(breaks=breaks_IceOff,labels=c("30-Apr","10-May","20-May","30-May","09-Jun"),limits=c(120,160))+
#   coord_cartesian(xlim=c(1925,2025))+
#   scale_x_continuous(breaks=seq(1930, 2020, 15))+
#   theme_pubr(border=TRUE, base_size=8)+
#   theme(plot.margin=unit(c(0,0.5,0,0.5), "lines"),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.x=element_blank())+
#   # axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2)) +
#   geom_text(data=panelLetter.normal,
#             aes(x=xpos,
#                 y=ypos,
#                 hjust=hjustvar,
#                 vjust=vjustvar,
#                 label="b",
#                 fontface="bold"))
# 
# IceDuration <-
#   ggplot() + 
#   geom_point(data=YSLoff, aes(x=Year,
#                               y=IceOnJulian-IceOffJulian), 
#              shape=21, fill="grey50", alpha=0.9)+
#   labs(x="Year",
#        y="Ice duration (days)")+
#   # scale_y_continuous(breaks=breaks_IceOff,labels=c("30-Apr","10-May","20-May","30-May","09-Jun"),limits=c(120,160))+
#   # coord_cartesian(xlim=c(1925,2025))+
#   scale_x_continuous(breaks=seq(1930, 2020, 15))+
#   theme_pubr(border=TRUE, base_size=8)+
#   theme(plot.margin=unit(c(0,0.5,0.5,0.5), "lines"))+
#   # axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2)) +
#   geom_text(data=panelLetter.normal,
#             aes(x=xpos,
#                 y=ypos,
#                 hjust=hjustvar,
#                 vjust=vjustvar,
#                 label="c",
#                 fontface="bold"))
# 
# 
# (IceOnTS / IceOffTS / IceDuration) & scale_x_continuous(breaks=c(1930,1960,1990,2020))
# ggsave("Figures/Figure2_IcePhenology.png", width=3, height=6,units="in", dpi=300)



# FIGURE 3  - ice on and off ----------------------------------------------------------------
##Vertically aligned
IceOn_FallMin_vert <- IceOn_FallMin +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.title.y=element_blank(),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))+
  labs(title="Ice-on")
IceOn_FallSnow_vert <- IceOn_FallSnow +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))
IceOn_FallTempSum_vert <- IceOn_FallTempSum +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))

IceOff_SpringSnow_vert <- IceOff_SpringSnow +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.title.y=element_blank(),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))+
  labs(title="Ice-off")
IceOff_SpringRain_vert <- IceOff_SpringRain +
  theme(plot.margin=unit(c(0,0.1,0.3,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))
IceOff_SnowDepth_vert <- IceOff_SnowDepth +
  theme(plot.margin=unit(c(0,0.1,0,0.5), "lines"),
        axis.ticks.y=element_line(),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.y = element_text(angle = 45, hjust=0.5, vjust=1.2))

combined <- (IceOn_FallMin_vert+ #a
               IceOff_SpringSnow_vert+ #d
               IceOn_FallSnow_vert+ #b
               IceOff_SpringRain_vert+ #e
               IceOn_FallTempSum_vert+ #c
               IceOff_SnowDepth_vert) & #f
  scale_fill_gradient( low = "white", high = "black",
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              frame.colour = "black",
                                              ticks = TRUE, 
                                              nbin = 10,
                                              label.position = "bottom",
                                              barwidth = 13,
                                              barheight = 1.3, 
                                              direction = 'horizontal')) &
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=8))

combined <- combined  +
  patchwork::plot_layout(ncol = 2, guides="collect")

combined
ggsave("Figures/Figure3_GAMS_IceOn_IceOff.pdf", width=5, height=7,units="in", dpi=600)
ggsave("Figures/Figure3_GAMS_IceOn_IceOff.png", width=5, height=7,units="in", dpi=600)


# FIGURE 4  - predictors - ice on and off ----------------------------------------------------------------
IceOff_nolegend <- IceOff &
  theme(legend.position="none")

#Row 1
Row1 <- (
  (GAMS_FallMin_new + GAMS_FallSnow_new + GAMS_FallTempSum_new) +
    patchwork::plot_layout(ncol = 1) &
    scale_x_continuous(breaks = c(1930, 1960, 1990, 2020))
) +
  labs(x = "Year")

#Row 2
Row2 <-
  (
    (GAMS_SpringSnow_new + GAMS_SpringRain_new + GAMS_SnowDepth_new) +  patchwork::plot_layout(ncol = 1) &
      scale_x_continuous(breaks = c(1930, 1960, 1990, 2020))
  )+
  labs(x="Year")

cowplot::plot_grid(Row1, Row2)

ggsave("Figures/Figure4_GAMS_IceOn_IceOff_Predictors.pdf", width=5, height=7,units="in", dpi=600)
ggsave("Figures/Figure4_GAMS_IceOn_IceOff_Predictors.png", width=5, height=7,units="in", dpi=600)



# TABLE S1 - GAMs summaries -----------------------------------------------


mod1_iceOn_df <- tidy(mod1_iceOn) %>%
  mutate(mod="mod1",
         devExp=summary(mod1_iceOn)$dev.expl) %>%
  left_join(., glance(mod1_iceOn) %>% mutate(mod="mod1",
                                             y="iceOn"))

mod2_iceOn_df <- tidy(mod2_iceOn) %>%
  mutate(mod="mod2",
         devExp=summary(mod2_iceOn)$dev.expl) %>%
  left_join(., glance(mod2_iceOn) %>% mutate(mod="mod2",
                                             y="iceOn"))

mod3_iceOn_df <- tidy(mod3_iceOn) %>%
  mutate(mod="mod3",
         devExp=summary(mod3_iceOn)$dev.expl) %>%
  left_join(., glance(mod3_iceOn) %>% mutate(mod="mod3",
                                             y="iceOn"))

mod4_iceOn_df <- tidy(mod4_iceOn) %>%
  mutate(mod="mod4",
         devExp=summary(mod4_iceOn)$dev.expl) %>%
  left_join(., glance(mod4_iceOn) %>% mutate(mod="mod4",
                                             y="iceOn"))

mod1_iceOff_df <- tidy(mod1_iceOff) %>%
  mutate(mod="mod1",
         devExp=summary(mod1_iceOff)$dev.expl) %>%
  left_join(., glance(mod1_iceOff) %>% mutate(mod="mod1",
                                              y="iceOff"))

mod2_iceOff_df <- tidy(mod2_iceOff) %>%
  mutate(mod="mod2",
         devExp=summary(mod2_iceOff)$dev.expl) %>%
  left_join(., glance(mod2_iceOff) %>% mutate(mod="mod2",
                                              y="iceOff"))

mod3_iceOff_df <- tidy(mod3_iceOff) %>%
  mutate(mod="mod3",
         devExp=summary(mod3_iceOff)$dev.expl) %>%
  left_join(., glance(mod3_iceOff) %>% mutate(mod="mod3",
                                              y="iceOff"))

mod4_iceOff_df <- tidy(mod4_iceOff) %>%
  mutate(mod="mod4",
         devExp=summary(mod4_iceOff)$dev.expl) %>%
  left_join(., glance(mod4_iceOff) %>% mutate(mod="mod4",
                                              y="iceOff"))

mod5_iceOff_df <- tidy(mod5_iceOff) %>%
  mutate(mod="mod5",
         devExp=summary(mod5_iceOff)$dev.expl) %>%
  left_join(., glance(mod5_iceOff) %>% mutate(mod="mod5",
                                              y="iceOff"))


mod_iceOn_df <- bind_rows(mod1_iceOn_df,
                          mod2_iceOn_df,
                          mod3_iceOn_df,
                          mod4_iceOn_df) %>%
  arrange(AIC)%>%
  mutate_if(is.numeric,
            round,
            digits = 3)

mod_iceOff_df <- bind_rows(mod1_iceOff_df,
                           mod2_iceOff_df,
                           mod3_iceOff_df,
                           mod4_iceOff_df,
                           mod5_iceOff_df) %>%
  arrange(AIC) %>%
  mutate_if(is.numeric,
            round,
            digits = 3)

all_mods <- bind_rows(mod_iceOn_df, mod_iceOff_df) %>%
  relocate(mod, .before=term) %>%
  relocate(y, .before=term) %>%
  select(-c(BIC:nobs), -df) 

all_mods_hux <- 
  hux(all_mods) %>% 
  # add_colnames() %>% 
  set_bold(row = 1, col = everywhere, value = TRUE) %>% 
  set_all_borders(TRUE) 

theme_plain(all_mods_hux)
quick_docx(all_mods_hux)

