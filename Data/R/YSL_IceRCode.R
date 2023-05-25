##Yellowstone Lake Ice On and Off timing
#Lusha's laptop
setwd("C:/Users/tronstad/Dropbox (UW WYNDD)/NC_YSL_Ice_Phenology/Data/R")

#Load data
Timing<-read.csv(file="YSL_Ice.csv",header=T,sep=",")
Weather<-read.csv(file="Yellowstone_Snow_Rain.csv",header=T,sep=",")
MonthWeatherOn<-read.csv(file="Data/R/Older/MonthlyWeather2IceOn.csv",header=T,sep=",")
MonthWeatherOff<-read.csv(file="Data/R/Older/MonthWeatherIceOff.csv",header=T,sep=",")
SeasonWeather<-read.csv(file="Data/R/Older/SeasonWeather.csv",header=T,sep=",")

#Load packages
library(plyr)
library(vegan)
library(Matrix)
library(corrplot)


#########################################
#Summarize weather data
#Annual
AnnualWeather<-ddply(Weather,.(Year),summarise,AnnualMax=max(MaxTempC,na.rm=T),
                     AnnualMin=min(MinTempC,na.rm=T),Rain=sum(Rain_mm,na.rm=T),
                     Snow=sum(Snow_mm,na.rm=T))
#Merge dataframes (Ice on)
YSL1<-merge(x=Timing,y=AnnualWeather,by="Year",all.x=T)
YSL2<-merge(x=YSL1,y=MonthWeatherOn,by="Year",all.x=T)
YSLon<-merge(x=YSL2,y=SeasonWeather,by="Year",all.x=T)

#Merge dataframes (Ice off)
YSL1<-merge(x=Timing,y=AnnualWeather,by="Year",all.x=T)
YSL2<-merge(x=YSL1,y=MonthWeatherOff,by="Year",all.x=T)
YSLoff<-merge(x=YSL2,y=SeasonWeather,by="Year",all.x=T)

min(YSLon$IceOnJulian,na.rm=T)
max(YSLon$IceOnJulian,na.rm=T)
mean(YSLon$IceOnJulian,na.rm=T)

min(YSLon$IceOffJulian)
max(YSLon$IceOffJulian)
mean(YSLon$IceOffJulian)

min(YSLon$AnnualMin)
max(YSLon$AnnualMin)
mean(YSLon$AnnualMin)

min(YSLon$AnnualMax)
max(YSLon$AnnualMax)
mean(YSLon$AnnualMax)

min(YSLon$Precip)
max(YSLon$Precip)
mean(YSLon$Precip)

min(YSLon$MaxSnowDepth_cm)
max(YSLon$MaxSnowDepth_cm)
mean(YSLon$MaxSnowDepth_cm)


#########################################
#Correlation
#Remove non-numeric columns
YSLM<-subset(YSLon,select=-c(IceOnDate,IceOffDate,IceOffJulian,IceOnJulian))
cor(YSLM)
corrplot(YSLM,method="circle",na.rm=T)




###########################################

#Time Series

library(prais)
#citation(package ="prais")

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

FitM<-lm(AnnualMax~Year,data=YSLon)
summary(FitM)


#Precipitation through time
fit.lm=lm(Precip~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(Precip~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

FitM<-lm(Precip~Year,data=YSLon)
summary(FitM)

#Max snow depth through time
fit.lm=lm(MaxSnowDepth_cm~Year,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(MaxSnowDepth_cm~Year,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

FitM<-lm(MaxSnowDepth_cm~Year,data=YSLon)
summary(FitM)


#Plot
Y5<-expression(paste("Maximum temperature (" , degree,   "C)"))
Y6<-expression(paste("Minimum temperature (" , degree,   "C)"))
par(mar=c(5,5,4,5)+.1)
plot(YSLon$Year,YSLon$AnnualMin,ylab=Y6,xlab="Year",cex.axis=1.5,cex.lab=2,
     ylim=c(-75,-25))
abline(-263.99,0.10761)
par(new=T)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
plot(YSLon$Year,YSLon$AnnualMax,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
     pch=15,col='black',axes=F,ylim=c(50,80))
#abline(44.9,0.009902,lty=2)
axis(side=4,cex.axis=1.5)
mtext(side=4,line=3,text=Y5,cex=2)
legend("topright",c("Maximum annual","Minimum annual"),pch=c(1,17),cex=1)

plot(YSLon$Year,YSLon$MaxSnowDepth_cm,ylab="Precipitation (cm)",xlab="Year",cex.axis=1.5,
     cex.lab=2,ylim=c(0,210))
abline(-371.86,0.244)
points(YSLon$Year,YSLon$Precip,pch=17,col="black")
abline(-315.73695,0.18464,lty=3)

##############################################################
#Ice-on models
##AIC values in fit.ar
#Annual temp and precip
fit.lm=lm(IceOnJulian~Year+AnnualMax+AnnualMin+Precip,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+AnnualMax+AnnualMin+Precip,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Annual temp and precip
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

#Mins Model
fit.lm=lm(IceOnJulian~Year+MinJul+MinAug+MinSep+MinOct+MinNov+MinDec,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinJul+MinAug+MinSep+MinOct+MinNov+MinDec,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model Summer Max
fit.lm=lm(IceOnJulian~Year+MaxJul+MaxAug+MaxSep,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MaxJul+MaxAug+MaxSep,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model Summer & Fall Max
fit.lm=lm(IceOnJulian~Year+MaxJul+MaxAug+MaxSep+MaxOct+MaxNov+MaxDec,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MaxJul+MaxAug+MaxSep+MaxOct+MaxNov+MaxDec,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

######
#HUGE MODEL of Precip (#1 model)
fit.lm=lm(IceOnJulian~Year+PrecipJan+PrecipFeb+PrecipMar+PrecipApr+PrecipMay+PrecipJun+PrecipJul+PrecipAug+ 
            PrecipSep+PrecipOct+PrecipNov+PrecipDec,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+PrecipJan+PrecipFeb+PrecipMar+PrecipApr+PrecipMay+PrecipJun+PrecipJul+PrecipAug+ 
                        PrecipSep+PrecipOct+PrecipNov+PrecipDec,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

FitIO<-lm(IceOnJulian~Year+PrecipJan+PrecipFeb+PrecipMar+PrecipApr+PrecipMay+PrecipJun+PrecipJul+PrecipAug+ 
     PrecipSep+PrecipOct+PrecipNov+PrecipDec,data=YSLon)
summary(FitIO)
######


#MOnthly temp  & precip
#December
fit.lm=lm(IceOnJulian~Year+MinDec+PrecipDec+MaxDec,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinDec+PrecipDec+MaxDec,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#November
fit.lm=lm(IceOnJulian~Year+MinNov+MaxNov+PrecipNov,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinNov+MaxNov+PrecipNov,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#August
fit.lm=lm(IceOnJulian~Year+MinAug+MaxAug+PrecipAug,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinAug+MaxAug+PrecipAug,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 1 (Best Model!)
fit.lm=lm(IceOnJulian~Year+MinMay+PrecipDec+PrecipNov+PrecipJul+PrecipJun+PrecipApr,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinMay+PrecipDec+PrecipNov+PrecipJul+PrecipJun+PrecipApr,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 2
fit.lm=lm(IceOnJulian~Year+MinMay+PrecipDec+PrecipNov+PrecipJul+PrecipJun,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinMay+PrecipDec+PrecipNov+PrecipJul+PrecipJun,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 3
fit.lm=lm(IceOnJulian~Year+MinMay+PrecipJun,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinMay+PrecipJun,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 4
fit.lm=lm(IceOnJulian~Year+PrecipJun,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+PrecipJun,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)


#Model 5
fit.lm=lm(IceOnJulian~Year+MinMay,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinMay,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 6
fit.lm=lm(IceOnJulian~Year+MinJul+PrecipJul,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinJul+PrecipJul,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 7
fit.lm=lm(IceOnJulian~Year+MinFall+MaxFall+PrecipFall,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinFall+MaxFall+PrecipFall,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 8
fit.lm=lm(IceOnJulian~Year+MinMay+MinJul+PrecipNov+PrecipDec,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinMay+MinJul+PrecipNov+PrecipDec,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 9
fit.lm=lm(IceOnJulian~Year+MinFall+MaxFall+PrecipFall,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinFall+MaxFall+PrecipFall,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 10
fit.lm=lm(IceOnJulian~Year+MinFall+MaxSummer+PrecipSummer+PrecipFall,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinFall+MaxSummer+PrecipSummer+PrecipFall,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 10
fit.lm=lm(IceOnJulian~Year+MinFall+MaxFall+MinSummer+MaxSummer+PrecipSummer+PrecipFall+MinSpring+MaxSpring+PrecipSpring,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MinFall+MaxFall+MinSummer+MaxSummer+PrecipSummer+PrecipFall+MinSpring+MaxSpring+PrecipSpring,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 11
fit.lm=lm(IceOnJulian~Year+PrecipSep+PrecipOct+PrecipNov,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+PrecipSep+PrecipOct+PrecipNov,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 12
fit.lm=lm(IceOnJulian~Year+MaxSnowDepth_cm,data=YSLon)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOnJulian~Year+MaxSnowDepth_cm,data=YSLon,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)




####################################################
#Ice-off models

#Full model Max
fit.lm=lm(IceOffJulian~Year+MaxJan+MaxFeb+MaxMar+MaxApr+MaxMay+MaxJun+
            MaxJul+MaxAug+MaxSep+MaxOct+MaxNov+MaxDec,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MaxJan+MaxFeb+MaxMar+MaxApr+MaxMay+MaxJun+
                        MaxJul+MaxAug+MaxSep+MaxOct+MaxNov+MaxDec,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Full model Min
fit.lm=lm(IceOffJulian~Year+MinJan+MinFeb+MinMar+MinApr+MinMay+MinJun+
            MinJul+MinAug+MinSep+MinOct+MinNov+MinDec,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinJan+MinFeb+MinMar+MinApr+MinMay+MinJun+
                        MinJul+MinAug+MinSep+MinOct+MinNov+MinDec,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Full model Precip
fit.lm=lm(IceOffJulian~Year+PrecipJan+PrecipFeb+PrecipMar+PrecipApr+PrecipMay+PrecipJun+
            PrecipJul+PrecipAug+PrecipSep+PrecipOct+PrecipNov+PrecipDec,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+PrecipJan+PrecipFeb+PrecipMar+PrecipApr+PrecipMay+PrecipJun+
                        PrecipJul+PrecipAug+PrecipSep+PrecipOct+PrecipNov+PrecipDec,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Annual weather
fit.lm=lm(IceOffJulian~Year+AnnualMax+AnnualMin+Precip,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+AnnualMax+AnnualMin+Precip,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 1
fit.lm=lm(IceOffJulian~Year+MinApr+MinMay+PrecipMay+MaxMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinApr+MinMay+PrecipMay+MaxMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 2
fit.lm=lm(IceOffJulian~Year+MinApr+MinMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinApr+MinMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 3
fit.lm=lm(IceOffJulian~Year+MaxApr+MaxMar+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MaxApr+MaxMar+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 4
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

#Model 5
fit.lm=lm(IceOffJulian~Year+MinMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 6
fit.lm=lm(IceOffJulian~Year+MinMay+PrecipJan,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinMay+PrecipJan,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 7
fit.lm=lm(IceOffJulian~Year+MinSpring+MaxSpring+PrecipSpring,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinSpring+MaxSpring+PrecipSpring,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 8
fit.lm=lm(IceOffJulian~Year+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

######
#Model 9 (best Model!)
fit.lm=lm(IceOffJulian~Year+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

FitOF<-lm(IceOffJulian~Year+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
summary(FitOF)

####

#Model 9
fit.lm=lm(IceOffJulian~Year+MinApr+MaxApr+PrecipApr,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinApr+MaxApr+PrecipApr,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

####

#Model 10
fit.lm=lm(IceOffJulian~Year+MaxSnowDepth_cm,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MaxSnowDepth_cm,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

####

#Model 11
fit.lm=lm(IceOffJulian~Year+MinNov+MinDec+MinJan+MinFeb+MinMar+MinApr+MinMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinNov+MinDec+MinJan+MinFeb+MinMar+MinApr+MinMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

####

#Model 12
fit.lm=lm(IceOffJulian~Year+MinApr+MinMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinApr+MinMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 13
fit.lm=lm(IceOffJulian~Year+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 14
fit.lm=lm(IceOffJulian~Year+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 15
fit.lm=lm(IceOffJulian~Year+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 15
fit.lm=lm(IceOffJulian~Year+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 16
fit.lm=lm(IceOffJulian~Year+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 17
fit.lm=lm(IceOffJulian~Year+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 17
fit.lm=lm(IceOffJulian~Year+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 18
fit.lm=lm(IceOffJulian~Year+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 19
fit.lm=lm(IceOffJulian~Year+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 20
fit.lm=lm(IceOffJulian~Year+MinJul+MaxJul+PrecipJul+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinJul+MaxJul+PrecipJul+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)

#Model 21
fit.lm=lm(IceOffJulian~Year+MinJun+MaxJun+PrecipJun+MinJul+MaxJul+PrecipJul+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff)
rlm=residuals(fit.lm)
fit.ar=arima(rlm,order=c(1,0,0),include.mean=F)
phi=fit.ar$coef[1]
require(prais)
fit.pw1=prais_winsten(IceOffJulian~Year+MinJun+MaxJun+PrecipJun+MinJul+MaxJul+PrecipJul+MinAug+MaxAug+PrecipAug+MinSep+MaxSep+PrecipSep+MinOct+MaxOct+PrecipOct+MinNov+MaxNov+PrecipNov+MinDec+MaxDec+PrecipDec+MinJan+MaxJan+PrecipJan+MinFeb+MaxFeb+PrecipFeb+MinMar+MaxMar+PrecipMar+MinApr+MaxApr+PrecipApr+MinMay+MaxMay+PrecipMay,data=YSLoff,iter=50,rho=phi)
fit.pw1
rpw1=residuals(fit.pw1)
fit2.ar=arima(rpw1,order=c(1,0,0),include.mean=F)
summary(fit.pw1)




#################################
#Plot data
par(mar=c(5,5,4,5)+.1)
Y6<-expression(paste("Ice-off Julian day"))
plot(YSLon$Year,YSLon$IceOnJulian,ylab="Ice-on Julian day",xlab="Year",cex.axis=1.5,cex.lab=2)
abline(295,.031)
par(new=T)
#points(YSLon$Year,YSLon$AnnualMax,pch=17,col="black")
plot(YSLoff$Year,YSLoff$IceOffJulian,ylab=NA,xlab=NA,cex.axis=1.5,cex.lab=2,
     pch=15,col='black',axes=F)
abline(148,-.0027,lty=2)
axis(side=4,cex.axis=1.5)
mtext(side=4,line=3,text=Y6,cex=2)

plot(YSLoff$Year,YSLoff$IceOffJulian,ylab="Ice-off Julian day",xlab="Year",cex.axis=1.5,cex.lab=2)
abline(302.8246,0.0258)


FitON<-lm(Year~IceOnJulian,data=YSLon)
library(lmtest)
dwtest(Year~IceOnJulian,data=YSLon)

FitON<-lm(Year~IceOffJulian,data=YSLoff)
library(lmtest)
dwtest(Year~IceOffJulian,data=YSLoff)

