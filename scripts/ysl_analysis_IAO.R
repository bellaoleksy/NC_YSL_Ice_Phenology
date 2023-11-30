library(tidyverse)
library(ggpubr)
library(mgcv)
library(gratia)

source("0_functions.R")


  


# Read in ice phenology data ----------------------------------------------


ysl_ice <- read.csv("scripts/YSLoff.csv")%>% #UPDATE FILEPATH HERE
  select(Year:SnowDepth) %>%
  mutate(IceOnJulian_new = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                     TRUE ~ IceOnJulian)) %>%
  mutate(new_iceondate = ymd(parse_date_time(paste(Year, IceOnJulian_new), orders = "yj")),
         new_iceoffdate = ymd(parse_date_time(paste(Year, IceOffJulian), orders = "yj")),
         
         IceOn_fedDOY = hydro.day(new_iceondate),
         IceOff_fedDOY = hydro.day(new_iceoffdate)) %>%
  select(-contains("Julian"))
  

#Get water-year for each ice phenology indicator separately, then join back together. 
ysl_iceOn <- ysl_ice %>%
  select(new_iceondate, IceOn_fedDOY) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceondate)) %>%
  filter(!new_iceondate=="1952-12-24") %>%
  arrange(water_year) %>%
  group_by(water_year) %>%
  slice_tail(n=1)

ysl_iceOn_old <- ysl_ice %>%
  select(new_iceondate, IceOn_fedDOY) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceondate)) %>%
  filter(!new_iceondate=="1952-12-24")

#Which years are duplicated?
dups_iceOn <- ysl_iceOn_old %>%
  group_by(water_year) %>% 
  filter(n()>1) 

# How many years have a duplicate ice-on date?
length(unique(dups_iceOn$water_year)) # !!!! ??? 


ysl_iceOff <- ysl_ice %>%
  select(new_iceoffdate, IceOff_fedDOY, AnnualMax:SnowDepth) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceoffdate))

#Which years are duplicated?
dups_iceOff <- ysl_iceOff %>%
  drop_na(water_year) %>%
  group_by(water_year) %>% 
  filter(n()>1) 
#Not an issue here luckily... 


ysl_ice_clean <- full_join(ysl_iceOn, ysl_iceOff, by="water_year") %>%
  mutate(iceDuration = IceOff_fedDOY-IceOn_fedDOY) 

ysl_ice_clean %>%
  select(water_year, IceOn_fedDOY, IceOff_fedDOY, iceDuration) %>%
  pivot_longer(-water_year) %>%
  ggplot(aes(x=water_year, y=value, color=name))+
  geom_point()+
  geom_smooth(method="gam")+
  facet_wrap(~name, scales="free")+
  theme_pubr()

#Indeed, no trend. Cool.



# Read in climate data ----------------------------------------------------


median(ysl_ice_clean$IceOn_fedDOY, na.rm=TRUE)
#Approx December 26 or so
median(ysl_ice_clean$IceOff_fedDOY, na.rm=TRUE)
#Approx May 22 or so

#Could use median ice on and ice off days to make decisions about seasonal calculations

Weather<-read.csv("Data/R/Yellowstone_Snow_Rain.csv",header=T,sep=",") %>%
  select(-OriginalOrder)

Seasonal_wx <- Weather %>%
  mutate(Date=mdy(Date),
         fedDOY=hydro.day(Date),
         water_year=dataRetrieval::calcWaterYear(Date)) %>%
  group_by(water_year) %>%
  mutate(season = case_when(Month %in% c("Nov","Dec","Jan") ~ "Nov_Dec_Jan",
                            Month %in% c("Mar","Apr","May") ~ "Mar_Apr_May")) %>%
  drop_na(season) %>%
  group_by(water_year, season) %>%
  summarize(temp_maxC=mean(max.C, na.rm=TRUE),
            temp_minC=mean(min.C, na.rm=TRUE),
            temp_meanC=mean(mean.C, na.rm=TRUE),
            cumul_rain=sum(rain.mm, na.rm=TRUE),
            cumul_snow=sum(snow.mm, na.rm=TRUE)) %>%
  filter(!water_year %in% c("1919","1978")) %>% #No data for these years, some values coming up as zero. Drop.
  pivot_wider(names_from="season", values_from = c(temp_maxC:cumul_snow),
              names_sep = "_")


# plot(Seasonal_wx$water_year, Seasonal_wx$cumul_snow_Nov_Dec_Jan)
# plot(Seasonal_wx$water_year, Seasonal_wx$cumul_snow_Mar_Apr_May)
# plot(Seasonal_wx$water_year, Seasonal_wx$temp_minC_Mar_Apr_May)
# plot(Seasonal_wx$water_year, Seasonal_wx$temp_meanC_Mar_Apr_May)
# plot(Seasonal_wx$water_year, Seasonal_wx$temp_meanC_Nov_Dec_Jan)



Snow_wx <- Weather %>%
  mutate(Date=mdy(Date),
         fedDOY=hydro.day(Date),
         water_year=dataRetrieval::calcWaterYear(Date)) %>%
  drop_na(SnowDepth.mm) %>%
  group_by(water_year) %>%
  summarize(cumul_snow_mm=sum(snow.mm, na.rm=TRUE),
            max_snowDepth_mm=max(SnowDepth.mm, na.rm=TRUE),
            DateMaxSnowDepth = min(Date[which(SnowDepth.mm == max(SnowDepth.mm))]),
            fedDOYMaxSnowDepth = hydro.day(DateMaxSnowDepth),
            MonthDayMaxSnowDepth = format(as.Date(DateMaxSnowDepth), "%m-%d"))

plot(Snow_wx$water_year, Snow_wx$max_snowDepth_mm)
plot(Snow_wx$water_year, Snow_wx$fedDOYMaxSnowDepth)

Wx_combined <- full_join(Snow_wx,Seasonal_wx)

# Trends in climate variables? --------------------------------------------
master_df <- ysl_ice_clean %>%
  left_join(., Wx_combined)

plot(master_df$cumul_snow_mm, master_df$AnnualSnow)
plot(master_df$max_snowDepth_mm, master_df$SnowDepth)
#Makes me wonder what SnowDepth is actually representing?
#Or how we calculated these things differently?

# ~ CUMUL. annual snow (by water year) ----------------------------------------------------------------


mod0_cumulSnow <- gam(cumul_snow_mm ~ s(water_year),
                 # family=Gamma(link="log"),
                 data = master_df,
                 # correlation = corCAR1(form = ~ water_year),
                 method = "REML")
summary(mod0_cumulSnow)
draw(mod0_cumulSnow)

# acf(residuals(mod0_cumulSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                       max(water_year, na.rm=TRUE),
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
Term = "water_year"
m1.d <- Deriv(mod0_cumulSnow) #in theory gratia::derivatives should work here

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

# ~ ANNUAL snow (original) ----------------------------------------------------------------


mod0_AnnualSnow <- gam(AnnualSnow ~ s(water_year),
                      # family=Gamma(link="log"),
                      data = master_df,
                      # correlation = corCAR1(form = ~ water_year),
                      method = "REML")
summary(mod0_AnnualSnow)
draw(mod0_AnnualSnow)

# acf(residuals(mod0_AnnualSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
                                                     length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
annualSnowPred <- cbind(years,
                       data.frame(predict(
                         mod0_AnnualSnow, years,
                         type = "response",
                         se.fit = TRUE
                       )))

### Calculate upper and lower bounds
annualSnowPred <- transform(annualSnowPred,
                           upper = fit + (2 * se.fit),
                           lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_AnnualSnow) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(annualSnowPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Positive rate of change in cumulative snow both the early part of the record and late part of the record


# ~ Max snow depth ----------------------------------------------------------------


mod0_maxSnowDepth <- gam(max_snowDepth_mm ~ s(water_year),
                      family=Gamma(link="log"),
                      data = master_df,
                      # correlation = corCAR1(form = ~ water_year),
                      method = "REML")
summary(mod0_maxSnowDepth)
draw(mod0_maxSnowDepth)

# acf(residuals(mod0_maxSnowDepth), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
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
Term = "water_year"
m1.d <- Deriv(mod0_maxSnowDepth) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxSnowDepthPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## Positive rate of change in cumulative snow both the early part of the record and late part of the record

# ~ DOY max snow depth ----------------------------------------------------------------


mod0_DOYmaxSnowDepth <- gam(fedDOYMaxSnowDepth ~ s(water_year),
                         family=Gamma(link="log"),
                         data = master_df,
                         # correlation = corCAR1(form = ~ water_year),
                         method = "REML")
summary(mod0_DOYmaxSnowDepth)
draw(mod0_DOYmaxSnowDepth)
#No trend in timing of maximum snow

# acf(residuals(mod0_DOYmaxSnowDepth), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation

# lm_mod0_DOYmaxSnowDepth <- lm(fedDOYMaxSnowDepth~water_year,
#                               data=master_df)



# ~ Cumul spring snow ----------------------------------------------------------------
hist(master_df$cumul_snow_Mar_Apr_May)

mod0_cumulSpringSnow <- gam(cumul_snow_Mar_Apr_May ~ s(water_year),
                         # family=Gamma(link="log"),
                         data = master_df,
                         # correlation = corCAR1(form = ~ water_year),
                         method = "REML")
summary(mod0_cumulSpringSnow)
draw(mod0_cumulSpringSnow)

# acf(residuals(mod0_cumulSpringSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
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
## No acceleration in the change

# ~ Cumul fall snow ----------------------------------------------------------------
hist(master_df$cumul_snow_Nov_Dec_Jan)

mod0_cumulFallSnow <- gam(cumul_snow_Nov_Dec_Jan ~ s(water_year),
                            # family=Gamma(link="log"),
                            data = master_df,
                            # correlation = corCAR1(form = ~ water_year),
                            method = "REML")
summary(mod0_cumulFallSnow)
draw(mod0_cumulFallSnow)

# acf(residuals(mod0_cumulFallSnow), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
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


# ~ Annual max temp. ----------------------------------------------------------------


mod0_maxTemp <- gam(AnnualMax ~ s(water_year),
                 family=Gamma(link="log"),
                 data = master_df,
                 correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                 #specifies the correlation argument of gam
                 method = "REML")
summary(mod0_maxTemp)
draw(mod0_maxTemp)

# acf(residuals(mod0_maxTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                         max(water_year, na.rm=TRUE),
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


# ~ Fall min temp. ----------------------------------------------------------
# Started here since minimum temps probably control ice formation more so than mean or max?

mod0_minFallTemp <- gam(temp_minC_Nov_Dec_Jan ~ s(water_year),
                 # family=Gamma(link="log"),
                 data = master_df,
                 correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                 #specifies the correlation argument of gam
                 method = "REML")
summary(mod0_minFallTemp)
draw(mod0_minFallTemp)

# acf(residuals(mod0_minFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                         max(water_year, na.rm=TRUE),
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
## Minimum fall temps have been getting warmer in the last few decades!



# ~ Fall max temp. ----------------------------------------------------------
# Not much to see here

mod0_maxFallTemp <- gam(temp_maxC_Nov_Dec_Jan ~ s(water_year),
                        # family=Gamma(link="log"),
                        data = master_df,
                        correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_maxFallTemp)
draw(mod0_maxFallTemp)

# acf(residuals(mod0_maxFallTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(max(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
                                                     length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
maxFallTempPred <- cbind(years,
                         data.frame(predict(
                           mod0_maxFallTemp, years,
                           type = "response",
                           se.fit = TRUE
                         )))

### Calculate upper and lower bounds
maxFallTempPred <- transform(maxFallTempPred,
                             upper = fit + (2 * se.fit),
                             lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_maxFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## maximum fall temps have been getting warmer in the last few decades!


# ~ Fall mean temp. ----------------------------------------------------------
# Not much to see here

mod0_meanFallTemp <- gam(temp_meanC_Nov_Dec_Jan ~ s(water_year),
                        # family=Gamma(link="log"),
                        data = master_df,
                        # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_meanFallTemp)
draw(mod0_meanFallTemp)

# acf(residuals(mod0_meanFallTemp), lag.mean = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(mean(water_year, na.rm=TRUE),
                                                     mean(water_year, na.rm=TRUE),
                                                     length.out = 200)))

#Create a dataframe with predicted ("fitted") values from the GAM and year, on the response scale.
meanFallTempPred <- cbind(years,
                         data.frame(predict(
                           mod0_meanFallTemp, years,
                           type = "response",
                           se.fit = TRUE
                         )))

### Calculate upper and lower bounds
meanFallTempPred <- transform(meanFallTempPred,
                             upper = fit + (2 * se.fit),
                             lower = fit - (2 * se.fit))

#Extract first derivative of the trend
Term = "water_year"
m1.d <- Deriv(mod0_meanFallTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(meanFallTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## mean fall temps accelerated 80s to 2000s



# ~ Spring min temp. ----------------------------------------------------------
# 

mod0_minSpringTemp <- gam(temp_minC_Mar_Apr_May ~ s(water_year),
                        # family=Gamma(link="log"),
                        data = master_df,
                        # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                        #specifies the correlation argument of gam
                        method = "REML")
summary(mod0_minSpringTemp)
draw(mod0_minSpringTemp)

# acf(residuals(mod0_minSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
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


# ~ Spring max temp. ----------------------------------------------------------
# 

mod0_maxSpringTemp <- gam(temp_maxC_Mar_Apr_May ~ s(water_year),
                          # family=Gamma(link="log"),
                          data = master_df,
                          # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                          #specifies the correlation argument of gam
                          method = "REML")
summary(mod0_maxSpringTemp)
draw(mod0_maxSpringTemp)

# acf(residuals(mod0_maxSpringTemp), lag.max = 10, plot = TRUE, main = "ACF of Residuals")
#Lag 1 temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(max(water_year, na.rm=TRUE),
                                                     max(water_year, na.rm=TRUE),
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
Term = "water_year"
m1.d <- Deriv(mod0_maxSpringTemp) #in theory gratia::derivatives should work here

#Calculate confidence intervals around the first derivative
m1.dci <- confint(m1.d, term = "water_year")

#Extract periods of increasing or decreasing trends
m1.dsig <- signifD(maxSpringTempPred$fit,
                   d = m1.d[[Term]]$deriv,
                   m1.dci[[Term]]$upper,
                   m1.dci[[Term]]$lower)

#Plot the first derivative 
plot.Deriv(m1.d)
## maximum Spring temps have been getting warmer in the last few decades!



# ~ ANNUAL min temps --------------------------------------------------------------


mod0_MinTemp <- gam(AnnualMin ~ s(water_year),
                    data = master_df,
                    # correlation = corARMA(form = ~ 1 | water_year, p = 1), 
                    #specifies the correlation argument of gam
                    method = "REML")
summary(mod0_MinTemp)
draw(mod0_MinTemp)

# acf(residuals(mod0_MinTemp),  plot = TRUE, main = "ACF of Residuals")
#No temporal autocorrelation


#Test for a change in the rate of change using 1st derivatives
#Extract years for next step
years <- with(master_df, data.frame(water_year = seq(min(water_year, na.rm=TRUE),
                                                         max(water_year, na.rm=TRUE),
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

# Ice Phenology -----------------------------------------------------------


# ~ Ice on ----------------------------------------------------------------


#Distribution of y
hist(master_df$IceOn_fedDOY)

### I added Family Gamma here since observations are always > 0
mod0_iceOn <- gam(IceOn_fedDOY ~ s(water_year),
                        family=Gamma(link="log"),
                        data = master_df,
                        # correlation = corCAR1(form = ~ water_year),
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
mod1_iceOn <- gam(IceOn_fedDOY ~ s(water_year) + s(temp_minC_Nov_Dec_Jan),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod1_iceOn)
draw(mod1_iceOn)
appraise(mod1_iceOn)


### I added Family Gamma here since observations are always > 0
mod2_iceOn <- gam(IceOn_fedDOY ~ s(water_year) + s(temp_meanC_Nov_Dec_Jan),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod2_iceOn)
draw(mod2_iceOn)
appraise(mod2_iceOn)

### I added Family Gamma here since observations are always > 0
mod3_iceOn <- gam(IceOn_fedDOY ~ s(water_year) + s(temp_maxC_Nov_Dec_Jan),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod3_iceOn)
draw(mod3_iceOn)
appraise(mod3_iceOn)


anova(mod2_iceOn, mod3_iceOn, test="F")

### I added Family Gamma here since observations are always > 0
mod4_iceOn <- gam(IceOn_fedDOY ~ s(water_year) + s(temp_maxC_Nov_Dec_Jan) + s(cumul_snow_Nov_Dec_Jan),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod4_iceOn)
draw(mod4_iceOn)
appraise(mod4_iceOn)


### I added Family Gamma here since observations are always > 0
mod5_iceOn <- gam(IceOn_fedDOY ~ s(water_year) + s(temp_maxC_Nov_Dec_Jan) + s(cumul_snow_Nov_Dec_Jan) + s(cumul_rain_Nov_Dec_Jan),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod5_iceOn)
draw(mod5_iceOn)
appraise(mod5_iceOn)


# ~ Ice off ----------------------------------------------------------------


#Distribution of y
hist(master_df$IceOff_fedDOY)

### I added Family Gamma here since observations are always > 0
mod0_iceOff <- gam(IceOff_fedDOY ~ s(water_year),
                  family=Gamma(link="log"),
                  data = master_df,
                  # correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod0_iceOff)
draw(mod0_iceOff)
appraise(mod0_iceOff)

#PLOT Autocorrelation function of residuals from the additive model with AR(1) errors
ACF <- acf(resid(mod0_iceOff, type = "response"), plot = FALSE)
ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
ggplot(ACF, aes(x = Lag, y = ACF)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = Lag, yend = 0))
#Suggests that an AR(1) model isn't necessary


mod1_iceOff <- gam(IceOff_fedDOY ~ s(water_year) + s(cumul_snow_mm),
                   family=Gamma(link="log"),
                   data = master_df,
                   # correlation = corCAR1(form = ~ water_year),
                   method = "REML")
summary(mod1_iceOff)
draw(mod1_iceOff)
appraise(mod1_iceOff)

mod2_iceOff <- gam(IceOff_fedDOY ~ s(water_year) + s(max_snowDepth_mm),
                   family=Gamma(link="log"),
                   data = master_df,
                   # correlation = corCAR1(form = ~ water_year),
                   method = "REML")
summary(mod2_iceOff)
draw(mod2_iceOff)
appraise(mod2_iceOff)

mod3_iceOff <- gam(IceOff_fedDOY ~ s(water_year) + s(fedDOYMaxSnowDepth),
                   family=Gamma(link="log"),
                   data = master_df,
                   # correlation = corCAR1(form = ~ water_year),
                   method = "REML")
summary(mod3_iceOff)
draw(mod3_iceOff)
appraise(mod3_iceOff)

mod4_iceOff <- gam(IceOff_fedDOY ~ s(water_year) + s(cumul_snow_Mar_Apr_May),
                   family=Gamma(link="log"),
                   data = master_df,
                   # correlation = corCAR1(form = ~ water_year),
                   method = "REML")
summary(mod4_iceOff)
draw(mod4_iceOff)
appraise(mod4_iceOff)

mod5_iceOff <- gam(IceOff_fedDOY ~ s(water_year) + s(cumul_snow_Mar_Apr_May) + s(temp_meanC_Mar_Apr_May),
                   family=Gamma(link="log"),
                   data = master_df,
                   # correlation = corCAR1(form = ~ water_year),
                   method = "REML")
summary(mod5_iceOff)
draw(mod5_iceOff)
appraise(mod5_iceOff)

# ~ Ice duration ----------------------------------------------------------


#Distribution of y
hist(master_df$iceDuration)

### I added Family Gamma here since observations are always > 0
mod0_iceDuration <- gam(iceDuration ~ s(water_year),
                          family=Gamma(link="log"),
                          data = master_df,
                          correlation = corCAR1(form = ~ water_year),
                          method = "REML")
summary(mod0_iceDuration)

#PLOT Autocorrelation function of residuals from the additive model with AR(1) errors
ACF <- acf(resid(mod0_iceDuration, type = "response"), plot = FALSE)
ACF <- setNames(data.frame(unclass(ACF)[c("acf", "lag")]), c("ACF","Lag"))
ggplot(ACF, aes(x = Lag, y = ACF)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = Lag, yend = 0))
#Suggests that an AR(1) model isn't necessary

mod1_iceDuration <- gam(iceDuration ~ s(water_year) + s(AnnualSnow),
                  family=Gamma(link="log"),
                  data = master_df,
                  correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod1_iceDuration)
appraise(mod1_iceDuration)
draw(mod1_iceDuration)
gam.check(mod1_iceDuration)


mod2_iceDuration <- gam(iceDuration ~ s(water_year) + s(AnnualSnow) + s(AnnualMin),
                  family=Gamma(link="log"),
                  data = master_df,
                  correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod2_iceDuration)
appraise(mod2_iceDuration)
draw(mod2_iceDuration)
gam.check(mod2_iceDuration)

cor.test(master_df$AnnualMin, master_df$AnnualMax)
mod3_iceDuration <- gam(iceDuration ~ s(water_year) + s(AnnualSnow) + s(AnnualMin) + s(AnnualMax),
                  family=Gamma(link="log"),
                  data = master_df,
                  correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod3_iceDuration)
appraise(mod3_iceDuration)
draw(mod3_iceDuration)
gam.check(mod3_iceDuration)

mod4_iceDuration <- gam(iceDuration ~ s(water_year) + s(AnnualSnow) + s(AnnualMin) + s(AnnualRain),
                  family=Gamma(link="log"),
                  data = master_df,
                  correlation = corCAR1(form = ~ water_year),
                  method = "REML")
summary(mod4_iceDuration)
appraise(mod4_iceDuration)
draw(mod4_iceDuration)
gam.check(mod4_iceDuration)

mod5_iceDuration <- gam(iceDuration ~ s(water_year) + s(cumul_snow_mm) + s(AnnualMin) + s(AnnualRain),
                        family=Gamma(link="log"),
                        data = master_df,
                        correlation = corCAR1(form = ~ water_year),
                        method = "REML")
summary(mod5_iceDuration)
appraise(mod5_iceDuration)
draw(mod5_iceDuration)
gam.check(mod5_iceDuration)

mod6_iceDuration <- gam(iceDuration ~ s(water_year) + s(cumul_snow_mm) + s(AnnualMin) + s(AnnualRain),
                        family=Gamma(link="log"),
                        data = master_df,
                        correlation = corCAR1(form = ~ water_year),
                        method = "REML")
summary(mod5_iceDuration)
appraise(mod5_iceDuration)
draw(mod5_iceDuration)
gam.check(mod5_iceDuration)
