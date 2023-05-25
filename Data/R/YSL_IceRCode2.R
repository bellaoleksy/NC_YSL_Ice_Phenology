##Yellowstone Lake Ice On and Off timing
#Lusha's desktop
setwd("G:/Dropbox (UW WYNDD)/Proj_NPS_YellowstoneLakeFoodWeb/Manuscript/IceOnOff/Data/R")
#Lusha's laptop
setwd("C:/Users/tronstad/Dropbox (UW WYNDD)/NC_YSL_Ice_Phenology/Data/R")

#Load data
Timing<-read.csv(file="YSL_Ice.csv",header=T,sep=",")
Weather<-read.csv(file="Yellowstone_Snow_Rain.csv",header=T,sep=",")
#MonthWeatherOn<-read.csv(file="MonthWeatherIceOn.csv",header=T,sep=",")
#MonthWeatherOff<-read.csv(file="MonthWeatherIceOff.csv",header=T,sep=",")
#SeasonWeather<-read.csv(file="SeasonWeather.csv",header=T,sep=",")

#Load packages
library(plyr)
library(vegan)
library(Matrix)
library(corrplot)


#########################################
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

Spring1<-ddply(Spring,.(Year),summarise,SpringMax=max(max.C,na.rm=T),SpringMin=min(min.C,na.rm=T),
                        SpringRain=sum(rain.mm,na.rm=T),SpringSnow=sum(snow.mm,na.rm=T))
Summer1<-ddply(Summer,.(Year),summarise,SummerMax=max(max.C,na.rm=T),SummerMin=min(min.C,na.rm=T),
               SummerRain=sum(rain.mm,na.rm=T),SummerSnow=sum(snow.mm,na.rm=T))
Fall1<-ddply(Fall,.(Year),summarise,FallMax=max(max.C,na.rm=T),FallMin=min(min.C,na.rm=T),
             FallRain=sum(rain.mm,na.rm=T),FallSnow=sum(snow.mm,na.rm=T))
Winter1<-ddply(Winter,.(Year),summarise,WinterMax=max(max.C,na.rm=T),WinterMin=min(min.C,na.rm=T),
               WinterRain=sum(rain.mm,na.rm=T),WinterSnow=sum(snow.mm,na.rm=T))
#Bind seasonal dataframes together
SeasonWeather1<-cbind(Spring1,Summer1,Fall1,Winter1)
#Remove the extra year columns
SeasonWeather<-subset(SeasonWeather1,select=-c(6,11,16))

YSLon<-merge(x=YSL1,y=SeasonWeather,by="Year",all.x=T)
#Replace Inf with NA
YSLon=do.call(data.frame,lapply(YSLon,function(value) replace(value,is.infinite(value),NA)))





#####Ice OFF!
#Spring
SpringF<-Weather[Weather$Month=="June"|Weather$Month=="Apr"|Weather$Month=="May",]
SummerF<-Weather[Weather$Month=="Sep"|Weather$Month=="Jul"|Weather$Month=="Aug",]
FallF<-Weather[Weather$Month=="Dec"|Weather$Month=="Oct"|Weather$Month=="Nov",]
WinterF<-Weather[Weather$Month=="Mar"|Weather$Month=="Jan"|Weather$Month=="Feb",]

Spring1F<-ddply(SpringF,.(Year),summarise,SpringMax=max(max.C,na.rm=T),SpringMin=min(min.C,na.rm=T),
               SpringRain=sum(rain.mm,na.rm=T),SpringSnow=sum(snow.mm,na.rm=T))
Summer1F<-ddply(SummerF,.(Year),summarise,SummerMax=max(max.C,na.rm=T),SummerMin=min(min.C,na.rm=T),
               SummerRain=sum(rain.mm,na.rm=T),SummerSnow=sum(snow.mm,na.rm=T))
Fall1F<-ddply(FallF,.(Year),summarise,FallMax=max(max.C,na.rm=T),FallMin=min(min.C,na.rm=T),
             FallRain=sum(rain.mm,na.rm=T),FallSnow=sum(snow.mm,na.rm=T))
Winter1F<-ddply(WinterF,.(Year),summarise,WinterMax=max(max.C,na.rm=T),WinterMin=min(min.C,na.rm=T),
               WinterRain=sum(rain.mm,na.rm=T),WinterSnow=sum(snow.mm,na.rm=T))
#Bind seasonal dataframes together
SeasonWeatherF<-cbind(Spring1F,Summer1F,Fall1F,Winter1F)
#Remove the extra year columns
SeasonWeatherF<-subset(SeasonWeatherF,select=-c(6,11,16))

YSLoff<-merge(x=YSL1,y=SeasonWeatherF,by="Year",all.x=T)
#Replace Inf with NA
YSLoff=do.call(data.frame,lapply(YSLon,function(value) replace(value,is.infinite(value),NA)))
#Export to shift summer and fall values
#Winter and spring occur before ice off of the current year
#Summer and fall occur after ice off happen during that calendar year
#So I shifted summer and winte fall values to the next year
#That way the 12 month before ice off occurs are represented
write.csv(YSLoff,"YSLoff.csv")
#The new file with summer and winter values shifted 1 year is called YSLoff_shifted
#USe this one for modeling!
YSLoff<-read.csv(file="YSLoff_Shifted.csv",header=T,sep=",")



#ICE ON
#Mins, maxs and means
min(YSLon$IceOnJulian,na.rm=T)
max(YSLon$IceOnJulian,na.rm=T)
mean(YSLon$IceOnJulian,na.rm=T)

min(YSLon$IceOffJulian,na.rm=T)
max(YSLon$IceOffJulian,na.rm=T)
mean(YSLon$IceOffJulian,na.rm=T)

min(YSLon$AnnualMin,na.rm=T)
max(YSLon$AnnualMin,na.rm=T)
mean(YSLon$AnnualMin,na.rm=T)

min(YSLon$AnnualMax,na.rm=T)
max(YSLon$AnnualMax,na.rm=T)
mean(YSLon$AnnualMax,na.rm=T)

min(YSLon$AnnualRain,na.rm=T)
max(YSLon$AnnualRain,na.rm=T)
mean(YSLon$AnnualRain,na.rm=T)

min(YSLon$AnnualSnow,na.rm=T)
max(YSLon$AnnualSnow,na.rm=T)
mean(YSLon$AnnualSnow,na.rm=T)

min(YSLon$SnowDepth,na.rm=T)
max(YSLon$SnowDepth,na.rm=T)
mean(YSLon$SnowDepth,na.rm=T)


#########################################
#Correlation
#Remove non-numeric columns
YSLM<-subset(YSLon,select=-c(IceOnDate,IceOffDate,IceOffJulian,IceOnJulian))
YSLM<-na.omit(YSLM)
cor(YSLM)
corrplot(YSLM,method="circle")




###########################################

#Time Series

library(prais)
library(lmtest)
citation(package ="prais")

#MIn temps through time
fit.lm=lm(AnnualMin~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(AnnualMin~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

dwtest(AnnualMin~Year,data=YSLon)

FitM<-lm(AnnualMin~Year,data=YSLon)
summary(FitM)


#Max temps through time
fit.lm=lm(AnnualMax~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(AnnualMax~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

dwtest(AnnualMax~Year,data=YSLon)

FitM<-lm(AnnualMax~Year,data=YSLon)
summary(FitM)


#Rain through time
fit.lm=lm(AnnualRain~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(AnnualRain~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

dwtest(AnnualRain~Year,data=YSLon)

FitM<-lm(AnnualRain~Year,data=YSLon)
summary(FitM)

#Snow through time
fit.lm=lm(AnnualSnow~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(AnnualSnow~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

dwtest(AnnualSnow~Year,data=YSLon)

FitM<-lm(AnnualSnow~Year,data=YSLon)
summary(FitM)

#Max snow depth through time
fit.lm=lm(SnowDepth~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(SnowDepth~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

dwtest(SnowDepth~Year,data=YSLon)

FitM<-lm(SnowDepth~Year,data=YSLon)
summary(FitM)


#Plot (export as 1000 wide)
Y5<-expression(paste("Maximum temperature (" , degree,   "C)"))
Y6<-expression(paste("Minimum temperature (" , degree,   "C)"))
par(mar=c(4,5,2,5)+.1)
plot(YSLon$Year,YSLon$AnnualMin,ylab=Y6,xlab="",cex.axis=1.5,cex.lab=2,
     ylim=c(-75,-20),xlim=c(1920,2020),pch=1)
abline(-158.7,0.06186,lty=1)
par(new=T)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
plot(YSLon$Year,YSLon$AnnualMax,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
     pch=15,col='black',axes=F,ylim=c(20,40),xlim=c(1920,2020))
#abline(17.0668,0.005498,lty=2)
axis(side=4,cex.axis=1.5)
mtext(side=4,line=3,text=Y5,cex=2)
legend("bottomleft",c("Minimum","Maximum"),pch=c(1,15),lty=c(1,0),cex=1.5)
title(main="a",adj=0,line=1)

#Next plot
Y7<-expression(paste("Precipitation (mm)"))
plot(YSLon$Year,YSLon$SnowDepth,ylab="Snow depth (mm)",xlab="",cex.axis=1.5,
     cex.lab=2,xlim=c(1920,2020),ylim=c(500,2100),pch=2)
title(main="b",adj=0,line=1)
#abline(-2731.3,2.195)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
#par(new=T)
#plot(YSLon$Year,YSLon$AnnualRain,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
 #    pch=15,col='black',axes=F,xlim=c(1920,2020),ylim=c(0,500))
#axis(side=4,cex.axis=1.5)
#mtext(side=4,line=3,text=Y7,cex=2)
#Third plot
plot(YSLon$Year,YSLon$AnnualRain,ylab=Y7,xlab="Year",cex.axis=1.5,cex.lab=2,
    pch=15,col='black',axes=T,xlim=c(1920,2020),ylim=c(0,500))
#abline(-78.8272,0.1439,lty=3)
points(YSLon$Year,YSLon$AnnualSnow,pch=8)
abline(-2483.6,1.3863,lty=2)
legend("topleft",c("Rain","Snow"),pch=c(15,8),lty=c(0,2),cex=1.5,horiz=F)
title(main="c",adj=0,line=1)



##############################################################
#Ice-on models
##AIC values in fit.ar
#Annual temp and precip
fit.lm=lm(IceOnJulian~Year+AnnualMax+AnnualMin+AnnualRain+AnnualSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+AnnualMax+AnnualMin+AnnualRain+AnnualSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Annual Min and Snow
fit.lm=lm(IceOnJulian~Year+AnnualMin+AnnualSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+AnnualMin+AnnualSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Year only
fit.lm=lm(IceOnJulian~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar


#Mins Model
fit.lm=lm(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Model  Max
fit.lm=lm(IceOnJulian~Year+SpringMax+SummerMax+FallMax+WinterMax,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringMax+SummerMax+FallMax+WinterMax,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Rain
fit.lm=lm(IceOnJulian~Year+SpringRain+SummerRain+FallRain+WinterRain,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringRain+SummerRain+FallRain+WinterRain,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Snow
fit.lm=lm(IceOnJulian~Year+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Summer and Autumn mins and snow
fit.lm=lm(IceOnJulian~Year+SummerMin+FallMin+FallSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SummerMin+FallMin+FallSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

# Autumn mins and snow
fit.lm=lm(IceOnJulian~Year+FallMin+FallSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+FallMin+FallSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar



######
#Full MODEL with everything
fit.lm=lm(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
          SpringRain+SummerRain+FallRain+WinterRain+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
            SpringRain+SummerRain+FallRain+WinterRain+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

######


#Full model without rain
fit.lm=lm(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
              SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
                          SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#All season Snow and Mins
fit.lm=lm(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#All variables for summer and fall
fit.lm=lm(IceOnJulian~Year+SummerMin+FallMin+SummerMax+FallMax+SummerRain+FallRain+SummerSnow+FallSnow,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+SummerMin+FallMin+SummerMax+FallMax+SummerRain+FallRain+SummerSnow+FallSnow,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar






####################################################
#Ice-off models
##AIC values in fit.ar

#Annual temp and precip
fit.lm=lm(IceOffJulian~Year+AnnualMax+AnnualMin+AnnualRain+AnnualSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+AnnualMax+AnnualMin+AnnualRain+AnnualSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Annual Min and Snow
fit.lm=lm(IceOffJulian~Year+AnnualMin+AnnualSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+AnnualMin+AnnualSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Year only
fit.lm=lm(IceOffJulian~Year,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Mins Model
fit.lm=lm(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Model  Max
fit.lm=lm(IceOffJulian~Year+SpringMax+SummerMax+FallMax+WinterMax,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMax+SummerMax+FallMax+WinterMax,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Rain
fit.lm=lm(IceOffJulian~Year+SpringRain+SummerRain+FallRain+WinterRain,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringRain+SummerRain+FallRain+WinterRain,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Snow
fit.lm=lm(IceOffJulian~Year+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#Winter and Spring mins and snow
fit.lm=lm(IceOffJulian~Year+WinterMin+SpringMin+SpringSnow+WinterSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+WinterMin+SpringMin+SpringSnow+WinterSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

# Spring mins and snow
fit.lm=lm(IceOffJulian~Year+SoringMin+SpringSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMin+SpringSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar



######
#Full MODEL with everything
fit.lm=lm(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
              SpringRain+SummerRain+FallRain+WinterRain+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
                          SpringRain+SummerRain+FallRain+WinterRain+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

######


#Full model without rain
fit.lm=lm(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
              SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringMax+SummerMax+FallMax+WinterMax+
                          SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#All season Snow and Mins
fit.lm=lm(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringMin+SummerMin+FallMin+WinterMin+SpringSnow+SummerSnow+FallSnow+WinterSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

#All variables for spring and winter
fit.lm=lm(IceOffJulian~Year+WinterMin+SpringMin+WinterMax+SpringMax+WinterRain+SpringRain+WinterSnow+SpringSnow,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+WinterMin+SpringMin+WinterMax+SpringMax+WinterRain+SpringRain+WinterSnow+SpringSnow,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar

# Spring rain and snow, winter & Spring mins
fit.lm=lm(IceOffJulian~Year+SpringRain+SpringSnow+WinterMin+SpringMin,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+SpringRain+SpringSnow+WinterMin+SpringMin,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)
fit2.ar





#################################
#Plot data
#Ice on and off
par(mar=c(5,5,2,5)+.1)
Y6<-expression(paste("Ice-off Julian day"))
plot(YSLon$Year,YSLon$IceOnJulian,ylab="Ice-on Julian day",xlab="Year",cex.axis=1.5,cex.lab=2)
abline(440,-.041)
par(new=T)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
plot(YSLoff$Year,YSLoff$IceOffJulian,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
     pch=15,col='black',axes=F)
abline(169,-.01,lty=2)
axis(side=4,cex.axis=1.5)
mtext(side=4,line=3,text=Y6,cex=2)
legend("topleft",c("Ice-on","Ice-off"),pch=c(1,15),lty=c(1,2),cex=1.5,horiz = T)
title(main="a",adj=0,line=1)

#Plot with out ablines
par(mar=c(5,5,2,5)+.1)
Y6<-expression(paste("Ice-off Julian day"))
plot(YSLon$Year,YSLon$IceOnJulian,ylab="Ice-on Julian day",xlab="Year",cex.axis=1.5,cex.lab=2)
#abline(440,-.041)
par(new=T)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
plot(YSLoff$Year,YSLoff$IceOffJulian,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
     pch=15,col='black',axes=F)
#abline(169,-.01,lty=2)
axis(side=4,cex.axis=1.5)
mtext(side=4,line=3,text=Y6,cex=2)
legend("topleft",c("Ice-on","Ice-off"),pch=c(1,15),cex=1.5,horiz = T)
title(main="a",adj=0,line=1)

plot(YSLoff$Year,YSLoff$IceOffJulian,ylab="Ice-off Julian day",xlab="Year",cex.axis=1.5,cex.lab=2)
#abline(302.8246,0.0258)


FitON<-lm(Year~IceOnJulian,data=YSLon)
library(lmtest)
dwtest(Year~IceOnJulian,data=YSLon)

FitON<-lm(Year~IceOffJulian,data=YSLoff)
library(lmtest)
dwtest(Year~IceOffJulian,data=YSLoff)

