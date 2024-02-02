

#Ice On

Timing<-read.csv(file="Data/R/YSL_Ice.csv",header=T,sep=",")
Weather<-read.csv(file="Data/R/Yellowstone_Snow_Rain.csv",header=T,sep=",")

yellowstone_on <- Timing %>%
  select(Year, IceOnDate, IceOnJulian) %>%
  mutate(IceOnJulian = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                 TRUE ~ IceOnJulian),
         IceOnDate = ymd(parse_date_time(paste(Year, IceOnJulian), orders = "yj")),
         water_year = calcWaterYear(IceOnDate)) %>%
  dplyr::rename(Year_on = Year) %>%
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
mean(yellowstone_phenology$ice_days, na.rm=TRUE)


#Annual weather by water-year
yellowstone_wx_wy <- Weather %>%
  mutate(
    Date = mdy(Date),
    year = year(Date),
    month = month(Date),
    month_name = month(Date, label = TRUE),
    water_year = dataRetrieval::calcWaterYear(Date)) %>% 
  group_by(water_year) %>%
  summarize(AnnualMax=max(max.C,na.rm=T),
            AnnualMin=min(min.C,na.rm=T),
            AnnualRain=sum(rain.mm,na.rm=T),
            AnnualSnow=sum(snow.mm,na.rm=T),
            SnowDepth=max(SnowDepth.mm,na.rm=T))

#Get rid of Inf values
na_strings <- c(-Inf,Inf)
yellowstone_wx_wy <- yellowstone_wx_wy %>% naniar::replace_with_na_all(condition = ~.x %in% na_strings)

Weather_2007 <- Weather %>%
  mutate(Date=mdy(Date), 
         water_year=calcWaterYear(Date)) %>%
  filter(water_year==2007)
#There's barely any data for this water_year... some November dates but nothing into the spring. 
#How much more data could be missing?

naniar::vis_miss(Weather)
#There's a LOT of missing data for SnowDepth.mm particularly toward the end of the dataset
#It makes me not want to rely on this too heavily, especially since we see an *apparent*
#decline in max snow depth later in the record. It's probably just because the data are spotty.

#Visualize the missing data over time
ggplot(Weather %>%
         mutate(Date=mdy(Date), 
                water_year=calcWaterYear(Date)), 
       aes(x = Date, 
           y = SnowDepth.mm)) + 
  naniar::geom_miss_point(alpha=0.05)

#Do we have that issue with say snow.mm?
ggplot(Weather %>%
         mutate(Date=mdy(Date), 
                water_year=calcWaterYear(Date)), 
       aes(x = Date, 
           y = snow.mm)) + 
  naniar::geom_miss_point(alpha=0.05)
#Not quite as bad... 


#Season weather by water-year
yellowstone_wx_seasons <- Weather %>%
  mutate(
    Date = mdy(Date),
    year = year(Date),
    month = month(Date),
    month_name = month(Date, label = TRUE),
    water_year = dataRetrieval::calcWaterYear(Date),
    water_year_corrected = case_when(month == '9' ~ water_year + 1,
                                     TRUE ~ water_year),
    #This is a little janky but we want to pretend that Sept is part of the same
    #water year as the next month (Oct) for the purposes of summarizing the data below
    season = case_when(Month %in% c("Mar","Apr","May") ~ "Spring",
                       Month %in% c("Jun","Jul","Aug") ~ "Summer",
                       Month %in% c("Sep","Oct","Nov") ~ "Fall",
                       Month %in% c("Dec","Jan","Feb") ~ "Winter")) %>% 
  group_by(water_year_corrected, season) %>%
  dplyr::summarize(  
    Max = max(max.C, na.rm = T),
    Min = min(min.C, na.rm = T),
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
  ) %>%
  dplyr::rename(water_year=water_year_corrected) 


#Join them all together
yellowstone_full <- yellowstone_phenology %>%
  left_join(., yellowstone_wx_wy) %>%
  left_join(., yellowstone_wx_seasons)
#Get rid of Inf values again
yellowstone_full <- yellowstone_full %>% naniar::replace_with_na_all(condition = ~.x %in% na_strings)


# ~ MODELS - ICE ON ----------------------------------------------------------------

### I added Family Gamma here since observations are always > 0
mod0_iceOn_new <- gam(j_on_wy ~ s(water_year),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod0_iceOn_new)
summary(mod0_iceOn)
draw(mod0_iceOn_new)
draw(mod0_iceOn)
appraise(mod0_iceOn_new)
appraise(mod0_iceOn)
#These results are different as with Lusha's code

### Start with all fall variables
mod1_iceOn_new <- gam(IceOnJulian ~ s(FallMin)+s(FallMax)+s(FallRain)+s(FallSnow)+
                    s(FallTempSum),
                  family=Gamma(link="log"),
                  data = yellowstone_full,
                  # correlation = corCAR1(form = ~ Year),
                  method = "REML")
summary(mod1_iceOn_new)
summary(mod1_iceOn)
draw(mod1_iceOn_new)
draw(mod1_iceOn)
appraise(mod1_iceOn_new)
appraise(mod1_iceOn_new)



# ~ MODELS - ICE OFF----------------------------------------------------------------

mod0_iceOff_new <- gam(j_off_wy ~ s(water_year),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod0_iceOff_new)
summary(mod0_iceOff)
draw(mod0_iceOff_new)
draw(mod0_iceOff_new)
appraise(mod0_iceOff_new)
appraise(mod0_iceOff_new)
#These are the same as with Lusha's code

mod1_iceOff_new <- gam(IceOffJulian ~ s(WinterSnow)+ s(SnowDepth) + s(SpringSnow)+ s(SpringRain) + s(WinterMin)+ s(SpringMin) + s(SpringMax),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod1_iceOff_new)
summary(mod1_iceOff)
draw(mod1_iceOff_new)
draw(mod1_iceOff)
appraise(mod1_iceOff_new)
appraise(mod1_iceOff)
#Identical. Great.

#Mod 5 is the one we went with in the paper
mod5_iceOff <- gam(IceOffJulian ~  s(SnowDepth) + s(SpringRain) +s(SpringSnow),
                   family=Gamma(link="log"),
                   data = yellowstone_full,
                   # correlatiOff = corCAR1(form = ~ Year),
                   method = "REML")
summary(mod5_iceOff)
draw(mod5_iceOff)
appraise(mod5_iceOff)

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
