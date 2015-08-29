# R Implementation of an integration of the microclimate model of Warren Porter's Niche Mapper system 
# Michael Kearney November 2013

# This version uses the Australia Water Availability Project (AWAP) daily 5km climate
# layers for Australia for air temperature, relative humidity, rainfall and cloud cover
# and uses monthly soil moisture estimates (splined) and the Australia Soils database to
# obtain soil properties, including their change through time due to soil moisture.
# Humidity is only from 1971 onward. Cloud cover is only from 1990 onward (and is based
# on daily solar layers relative to clear sky calculations from NicheMapR).
# It also uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
# Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
# Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
# by choosing the option 'rungads<-1' 

# required R packages
# raster
# sp
# ncdf
# XML
# dismo
# chron
# rgdal
# zoo
# RODBC

spatial<-"c:/Australian Environment/" # place where climate input files are kept

############## location and climatic data  ###################################

sitemethod <- 0 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
longlat <- c(147.2686,-36.901)#ITEXu c(147.25,-36.92) #Hotham c(147.14,-36.972), #Falls Creek c(147.274,-36.868), #Buller c(146.423,-37.146)       c(147.141,-36.984)
loc <- "Broken Hill, NSW, Australia" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
rungads<-1 # use the Global Aerosol Database?
dailywind<-1 # use daily windspeed database?
terrain<-1 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
snowmodel<-1 # run snow version? (slower!)
ystart <- 2004# start year for weather generator calibration dataset or AWAP database
yfinish <- 2012# end year for weather generator calibration dataset
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

############# microclimate model parameters ################################
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next for parameters are segmented velocity profiles due to bushes, rocks etc. on the surface, IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO!
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 1.25 # Integrator error for soil temperature calculations
DEP <- c(0., 1.5,  2.5, 5.,  10., 15.,  30., 50.,  100., 200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 1300 # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 20 # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
fieldcap<-40 # field capacity, mm
wilting<-5 # wilting point, mm
rainmult<-0.5 # rain multiplier for surface soil moisture (use to induce runoff), proportion
REFL<-0.10 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-0 # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise              										
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-50. # maximum available shade (%)
runshade<-1. # run the model twice, once for each shade level (1) or just for the first shade level (0)?
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
grasshade<-0 # this drives min shade value by the relative soil moisture multiplied by the maxshade parameter, above
Usrhyt <- 0.5# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.325 # snow density (mg/m3)
snowmelt<-1 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1.2 # undercatch multipier for converting rainfall to snow
rainmelt<-0.013 # paramter in equation that melts snow with rainfall as a function of air temp
write_input<-1 # write csv files of final input to working directory? 1=yes, 0=no.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value

# run the model
setwd('/git/micro_australia/')
  niche<-list(loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,fieldcap=fieldcap,wilting=wilting,rainmult=rainmult,runshade=runshade,grasshade=grasshade)
  source('NicheMapR_Setup_micro.R')
  nicheout<-NicheMapR(niche)
  
  # get output
  metout<-as.data.frame(nicheout$metout[1:(365*24*nyears),]) # above ground microclimatic conditions, min shade
  shadmet<-as.data.frame(nicheout$shadmet[1:(365*24*nyears),]) # above ground microclimatic conditions, max shade
  soil<-as.data.frame(nicheout$soil[1:(365*24*nyears),]) # soil temperatures, minimum shade
  shadsoil<-as.data.frame(nicheout$shadsoil[1:(365*24*nyears),]) # soil temperatures, maximum shade
  rainfall<-as.data.frame(nicheout$RAINFALL)
  MAXSHADES<-as.data.frame(nicheout$MAXSHADES)
  elev<-as.numeric(nicheout$ALTT)
  REFL<-as.numeric(nicheout$REFL)
  longlat<-as.matrix(nicheout$longlat)
  fieldcap<-as.numeric(nicheout$fieldcap)
  wilting<-as.numeric(nicheout$wilting)
  ectoin<-rbind(elev,REFL,longlat,fieldcap,wilting,ystart,yfinish)
  
  write.csv(metout,'metout.csv')
  write.csv(shadmet,'shadmet.csv')
  write.csv(soil,'soil.csv')
  write.csv(shadsoil,'shadsoil.csv')
  write.csv(rainfall,'rainfall.csv')
  write.csv(ectoin,'ectoin.csv')
  write.csv(DEP,'DEP.csv')
  write.csv(MAXSHADES,'MAXSHADES.csv')
  
  tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
  dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
  dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
  soil<-cbind(dates,soil)
  metout<-cbind(dates,metout)
  shadsoil<-cbind(dates,shadsoil)
  shadmet<-cbind(dates,shadmet)
  
  dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
  dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
  rainfall<-as.data.frame(cbind(dates2,rainfall))
  colnames(rainfall)<-c("dates","rainfall")
  
  years_snow<-c('90','91','92','93','94','95','96','97','98','99','00','01','02','03','04','05','06','07','08','09')
  metout<-as.data.frame(cbind(dates,metout[1:(8760*nyears),]))
  soil<-as.data.frame(cbind(dates,soil[1:(8760*nyears),]))
  shadmet<-as.data.frame(cbind(dates,shadmet[1:(8760*nyears),]))
  shadsoil<-as.data.frame(cbind(dates,shadsoil[1:(8760*nyears),]))
  

# test against ITEX data
load("C:/NicheMapR_Working/projects/snow/ITEXTempData_new.rda")
load("C:/NicheMapR_Working/projects/snow/TempData.rda")
load("C:/NicheMapR_Working/projects/snow/JCmicroclimateDate.rda")
load("C:/NicheMapR_Working/projects/snow/cover.rda")
TempData2<-ITEXlogger

yrs<-as.numeric(substr(TempData2$Date,1,4))
mths<-as.numeric(substr(TempData2$Date,6,7))
days<-as.numeric(substr(TempData2$Date,9,10))
hrs<-as.numeric(TempData2$Hour)
dateobs<-as.data.frame(cbind(yrs,mths,days,hrs))

dtime<-with(dateobs,ISOdatetime(yrs,mths,days,hrs,0,0))
TempData2<-cbind(dtime,TempData2)

yrs<-as.numeric(substr(microstation$Date,1,4))
mths<-as.numeric(substr(microstation$Date,6,7))
days<-as.numeric(substr(microstation$Date,9,10))
hrs<-as.numeric(microstation$Hour)
dateobs<-as.data.frame(cbind(yrs,mths,days,hrs))

dtime<-with(dateobs,ISOdatetime(yrs,mths,days,hrs,0,0))
microstation<-cbind(dtime,microstation)


ITEX3B<-subset(itex,SITE=="OZ3B" & Warming==0)
ITEX3B<-aggregate(ITEX3B$BG,by=list(ITEX3B$YEAR),mean)
ITEX4B<-subset(itex,SITE=="OZ4B" & Warming==0)
ITEX4B<-aggregate(ITEX4B$BG,by=list(ITEX4B$YEAR),mean)
ITEX1U<-subset(itex,SITE=="OZ1U" & Warming==0)
ITEX1U<-aggregate(ITEX1U$BG,by=list(ITEX1U$YEAR),mean)
ITEX2U<-subset(itex,SITE=="OZ2U" & Warming==0)
ITEX2U<-aggregate(ITEX2U$BG,by=list(ITEX2U$YEAR),mean)


addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


offset<-21
years<-seq(2004,2012,1)
#years<-seq(2008,2008,1)
month1s<-c("01","03","05","07","09","11")
month2s<-c("02","04","06","08","10","12")
TempData<-subset(TempData2,Treatment!='OTC')

TempData_unburnt<-subset(TempData,Site=="ITEX1U" | Site=="ITEX2U")
TempData_burnt<-subset(TempData,Site!="ITEX1U" & Site!="ITEX2U")

site_amb<-aggregate(TempData$CAmbient,by=list(TempData$dtime),mean,na.rm=TRUE)  
site_surf<-aggregate(TempData$CSurface,by=list(TempData$dtime),mean,na.rm=TRUE)  
site_5cm<-aggregate(TempData$CTemp_5cm,by=list(TempData$dtime),mean,na.rm=TRUE)  
site_10cm<-aggregate(TempData$CTemp_10cm,by=list(TempData$dtime),mean,na.rm=TRUE)

# site_amb_orig<-site_amb
# site_surf_orig<-site_surf 
# site_5cm_orig<-site_5cm  
# site_10cm_orig<-site_10cm

# site_amb_burnt<-aggregate(TempData_burnt$CAmbient,by=list(TempData_burnt$dtime),mean,na.rm=TRUE)  
# site_surf_burnt<-aggregate(TempData_burnt$CSurface,by=list(TempData_burnt$dtime),mean,na.rm=TRUE)  
# site_5cm_burnt<-aggregate(TempData_burnt$CTemp_5cm,by=list(TempData_burnt$dtime),mean,na.rm=TRUE)  
# site_10cm_burnt<-aggregate(TempData_burnt$CTemp_10cm,by=list(TempData_burnt$dtime),mean,na.rm=TRUE)
# 
# site_amb_unburnt<-aggregate(TempData_unburnt$CAmbient,by=list(TempData_unburnt$dtime),mean,na.rm=TRUE)  
# site_surf_unburnt<-aggregate(TempData_unburnt$CSurface,by=list(TempData_unburnt$dtime),mean,na.rm=TRUE)  
# site_5cm_unburnt<-aggregate(TempData_unburnt$CTemp_5cm,by=list(TempData_unburnt$dtime),mean,na.rm=TRUE)  
# site_10cm_unburnt<-aggregate(TempData_unburnt$CTemp_10cm,by=list(TempData_unburnt$dtime),mean,na.rm=TRUE)

# with(site_5cm,plot(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col="black",type='l'))
# with(site_5cm_burnt,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
# 
# with(site_5cm,plot(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col="black",type='l'))
# with(site_5cm_unburnt,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
# 
# with(site_5cm_unburnt,plot(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col="black",type='l'))
# with(site_5cm_burnt,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
# 
# with(site_5cm,plot(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col="black",type='l'))
# with(site_5cm_unburnt,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))


# 
# site_amb<-site_amb_burnt 
# site_surf<-site_surf_burnt  
# site_5cm<-site_5cm_burnt 
# site_10cm<-site_10cm_burnt
# 
# 
# site_amb<-site_amb_unburnt 
# site_surf<-site_surf_unburnt  
# site_5cm<-site_5cm_unburnt 
# site_10cm<-site_10cm_unburnt

site_amb<-site_amb_orig 
site_surf<-site_surf_orig  
site_5cm<-site_5cm_orig 
site_10cm<-site_10cm_orig

# # test differences
# years1<-c("05","06","07","08","09","10","11")
# for(i in 1:length(years1)){
# burnt<-subset(site_surf_burnt,format(site_surf_burnt$Group.1, "%y")==years1[i])
# unburnt<-subset(site_surf_unburnt,format(site_surf_unburnt$Group.1, "%y")==years1[i])
# merg<-merge(burnt,unburnt,by="Group.1")
# t<-t.test(merg$x.x,merg$x.y,paired=TRUE,na.rm=TRUE)
# ub<-mean(unburnt$x,na.rm=TRUE)
# b<-mean(burnt$x,na.rm=TRUE)
# if(i==1){
#   surf<-cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub)
# }else{
#   surf<-rbind(surf,cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub))
# }
# }
# colnames(surf)<-c("year","p","n","t","burnt","unburnt")
# 
# years1<-c("05","06","07","08","09","10","11")
# for(i in 1:length(years1)){
#   burnt<-subset(site_amb_burnt,format(site_amb_burnt$Group.1, "%y")==years1[i])
#   unburnt<-subset(site_amb_unburnt,format(site_amb_unburnt$Group.1, "%y")==years1[i])
#   merg<-merge(burnt,unburnt,by="Group.1")
#   t<-t.test(merg$x.x,merg$x.y,paired=TRUE,na.rm=TRUE)
#   ub<-mean(unburnt$x,na.rm=TRUE)
#   b<-mean(burnt$x,na.rm=TRUE)
#   if(i==1){
#     amb<-cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub)
#   }else{
#     amb<-rbind(amb,cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub))
#   }
# }
# colnames(amb)<-c("year","p","n","t","burnt","unburnt")
# 
# 
# years1<-c("05","06","07","08","09","10","11")
# for(i in 1:length(years1)){
#   burnt<-subset(site_5cm_burnt,format(site_5cm_burnt$Group.1, "%y")==years1[i])
#   unburnt<-subset(site_5cm_unburnt,format(site_5cm_unburnt$Group.1, "%y")==years1[i])
#   merg<-merge(burnt,unburnt,by="Group.1")
#   t<-t.test(merg$x.x,merg$x.y,paired=TRUE,na.rm=TRUE)
#   ub<-mean(unburnt$x,na.rm=TRUE)
#   b<-mean(burnt$x,na.rm=TRUE)
#   if(i==1){
#     d5cm<-cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub)
#   }else{
#     d5cm<-rbind(d5cm,cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub))
#   }
# }
# colnames(d5cm)<-c("year","p","n","t","burnt","unburnt")
# 
# 
# years1<-c("05","06","07","08","09","10","11")
# for(i in 1:length(years1)){
#   burnt<-subset(site_10cm_burnt,format(site_10cm_burnt$Group.1, "%y")==years1[i])
#   unburnt<-subset(site_10cm_unburnt,format(site_10cm_unburnt$Group.1, "%y")==years1[i])
#   merg<-merge(burnt,unburnt,by="Group.1")
#   t<-t.test(merg$x.x,merg$x.y,paired=TRUE,na.rm=TRUE)
#   ub<-mean(unburnt$x,na.rm=TRUE)
#   b<-mean(burnt$x,na.rm=TRUE)
#   if(i==1){
#     d10cm<-cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub)
#   }else{
#     d10cm<-rbind(d10cm,cbind(years1[i],t$p.value,t$parameter,t$statistic,b,ub))
#   }
# }
# colnames(d10cm)<-c("year","p","n","t","burnt","unburnt")

setwd("c:/NicheMapR_Working/projects/snow/ITEX")


for(i in 1:length(years)){
  doyear<-years[i]
  doyear2<-substr(doyear,3,4)
  #tzone<-paste("Etc/GMT+",10,sep="") # doing it this way ignores daylight savings!
  #dtime<-ISOdatetime(as.numeric(substr(TempData$Date, 7, 10)),as.numeric(substr(TempData$Date, 4, 5)),as.numeric(substr(TempData$Date, 1, 2)),TempData$Hour,0,0,tz=tzone)
  #TempData2<-cbind(dtime,TempData)
  #save(TempData2,file="C:/NicheMapR_Working/projects/snow/TempData2.rda")
  for(j in 1:length(month1s)){
    #plotsoil2<-subset(shadsoil,(format(shadsoil$dates, "%m")=='01' | format(shadsoil$dates, "%m")=='02') & format(shadsoil$dates, "%y")==doyear2)
    #plotmetout2<-subset(shadmet,(format(shadsoil$dates, "%m")=='01' | format(shadmet$dates, "%m")=='02') & format(shadmet$dates, "%y")==doyear2)
    #plotsoil2<-subset(shadsoil,format(shadsoil$dates, "%y")==doyear2)
    #plotmetout2<-subset(shadmet,format(shadmet$dates, "%y")==doyear2)
    plotsoil2<-subset(shadsoil,(format(shadsoil$dates, "%m")==month1s[j] | format(shadsoil$dates, "%m")==month2s[j]) & format(shadsoil$dates, "%y")==doyear2)
    plotmetout2<-subset(metout,(format(metout$dates, "%m")==month1s[j] | format(metout$dates, "%m")==month2s[j]) & format(metout$dates, "%y")==doyear2)
    plotshadmet2<-subset(metout,(format(shadmet$dates, "%m")==month1s[j] | format(shadmet$dates, "%m")==month2s[j]) & format(shadmet$dates, "%y")==doyear2)
    
    #plotsoil2<-subset(soil,format(soil$dates, "%y")==doyear2)
    #plotmetout2<-subset(metout,format(metout$dates, "%y")==doyear2)
    # remove heating treatments
    
    
    # site1<-subset(TempData,Site=="ITEX3B" & Plot !=700 & as.numeric(substr(TempData$Date, 1, 4)) == doyear)
    # site2<-subset(TempData,Site=="ITEX4B" & Plot !=700 & as.numeric(substr(TempData$Date, 1, 4)) == doyear)
    # site3<-subset(TempData,Site=="ITEX1U" & Plot !=700 & as.numeric(substr(TempData$Date, 1, 4)) == doyear)
    # site4<-subset(TempData,Site=="ITEX2U" & Plot !=700 & as.numeric(substr(TempData$Date, 1, 4)) == doyear)
    # dummy<-as.data.frame(cbind(as.POSIXct(site1[1:10,1]),rep(NA,10)))
    # dummy<-site1[1:10,1:2]
    # colnames(dummy)<-c('Group.1','x')
    # dummy$x<-NA
    # if(nrow(site1)>0){
    #   site1_amb<-aggregate(site1$CAmbient,by=list(site1$dtime),mean,na.rm=TRUE)  
    #   site1_surf<-aggregate(site1$CSurface,by=list(site1$dtime),mean,na.rm=TRUE)  
    #   site1_5cm<-aggregate(site1$CTemp_5cm,by=list(site1$dtime),mean,na.rm=TRUE)  
    #   site1_10cm<-aggregate(site1$CTemp_10cm,by=list(site1$dtime),mean,na.rm=TRUE)
    # }else{
    #   site1_amb<-dummy
    #   site1_surf<-dummy
    #   site1_5cm<-dummy
    #   site1_10cm<-dummy
    # }
    # if(nrow(site2)>0){
    #   site2_amb<-aggregate(site2$CAmbient,by=list(site2$dtime),mean,na.rm=TRUE)  
    #   site2_surf<-aggregate(site2$CSurface,by=list(site2$dtime),mean,na.rm=TRUE)  
    #   site2_5cm<-aggregate(site2$CTemp_5cm,by=list(site2$dtime),mean,na.rm=TRUE)  
    #   site2_10cm<-aggregate(site2$CTemp_10cm,by=list(site2$dtime),mean,na.rm=TRUE)
    # }else{
    #   site2_amb<-dummy
    #   site2_surf<-dummy
    #   site2_5cm<-dummy
    #   site2_10cm<-dummy
    # }
    # if(nrow(site3)>0){
    #   site3_amb<-aggregate(site3$CAmbient,by=list(site3$dtime),mean,na.rm=TRUE)  
    #   site3_surf<-aggregate(site3$CSurface,by=list(site3$dtime),mean,na.rm=TRUE)  
    #   site3_5cm<-aggregate(site3$CTemp_5cm,by=list(site3$dtime),mean,na.rm=TRUE)  
    #   site3_10cm<-aggregate(site3$CTemp_10cm,by=list(site3$dtime),mean,na.rm=TRUE)
    # }else{
    #   site3_amb<-dummy
    #   site3_surf<-dummy
    #   site3_5cm<-dummy
    #   site3_10cm<-dummy
    # }
    # if(nrow(site4)>0){
    #   site4_amb<-aggregate(site4$CAmbient,by=list(site4$dtime),mean,na.rm=TRUE)  
    #   site4_surf<-aggregate(site4$CSurface,by=list(site4$dtime),mean,na.rm=TRUE)  
    #   site4_5cm<-aggregate(site4$CTemp_5cm,by=list(site4$dtime),mean,na.rm=TRUE)  
    #   site4_10cm<-aggregate(site4$CTemp_10cm,by=list(site4$dtime),mean,na.rm=TRUE)
    # }else{
    #   site4_amb<-dummy
    #   site4_surf<-dummy
    #   site4_5cm<-dummy
    #   site4_10cm<-dummy
    # }
    # filename<-paste("pred_",doyear,"mon",month1s[j],"_",month2s[j],".pdf",sep="")
    # 
    # pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
    # par(mfrow = c(2,2)) # set up for 5 plots in 2 columns
    # par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
    # par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
    # 
    # with(plotshadmet2,plot(TALOC~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="0.5 cm air"))
    # with(plotsoil2,plot(D1.5cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="surface temp"))
    # with(plotsoil2,plot(D10cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="5cm soil 50% shade"))
    # with(plotsoil2,plot(D15cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="10cm soil 50% shade"))
    # 
    # dev.off()
    
    filename<-paste("pred_obs_",doyear,"mon",month1s[j],"_",month2s[j],".pdf",sep="")
    
    pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
    par(mfrow = c(2,2)) # set up for 5 plots in 2 columns
    par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
    par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
    
    
    #with(plotsoil,plot(DEP1~dates,type='l',col='1',ylab = "Snow Depth (cm)",xlab = "Year"))
    with(plotmetout2,plot(TALOC~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="0.5 cm air"))
    #with(plotsoil2,plot(D0cm~dates,type='l',col='1',ylab = "Snow Depth (cm)",xlab = "Year"))
    if(doyear>2009){
      with(site_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site1_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site2_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
      #   with(site3_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
      #   with(site4_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    }else{
      with(site_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site1_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site2_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
      #   with(site3_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
      #   with(site4_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    }
    #with(plotsoil,plot(DEP1~dates,type='l',col='1',ylab = "Snow Depth (cm)",xlab = "Year"))
    with(plotmetout2,plot(TS~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="surface temp"))
    #with(plotsoil2,plot(D0cm~dates,type='l',col='1',ylab = "Snow Depth (cm)",xlab = "Year"))
    if(doyear>2009){
      with(site_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site1_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site2_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
      #   with(site3_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
      #   with(site4_amb,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    }else{
      with(site_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site1_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
      #   with(site2_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
      #   with(site3_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
      #   with(site4_surf,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    }
    
    with(plotsoil2,plot(D5cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="5cm soil 50% shade"))
    with(site_5cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    # with(site1_5cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    # with(site2_5cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
    # with(site3_5cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
    # with(site4_5cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    
    with(plotsoil2,plot(D10cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="10cm soil 50% shade"))
    with(site_10cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    # with(site1_10cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    # with(site2_10cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("orange",150),type='l'))
    # with(site3_10cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark green",150),type='l'))
    # with(site4_10cm,points(x~as.POSIXct(Group.1+3600*offset),cex=0.1,col=addTrans("dark blue",150),type='l')) 
    dev.off()
  } #end months
} #end years






offset<-21
years<-seq(2010,2013,1)
#years<-seq(2012,2012,1)
month1s<-c("01","03","05","07","09","11")
month2s<-c("02","04","06","08","10","12")
setwd("c:/NicheMapR_Working/projects/snow/camac")
microstation2<-subset(microstation, Site=="ITEX1U")
Ambient_CTemp<-aggregate(microstation2$Ambient_CTemp,by=list(microstation2$dtime),mean,na.rm=TRUE)
CTemp_3cmBG<-aggregate(microstation2$CTemp_3cmBG,by=list(microstation2$dtime),mean,na.rm=TRUE)

for(i in 1:length(years)){
  doyear<-years[i]
  doyear2<-substr(doyear,3,4)
  
  filename<-paste("pred_",doyear,"_ambient_ITEX1U.pdf",sep="")
  
  pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,2)) # set up for 5 plots in 2 columns
  par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
  par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
  
  for(j in 1:length(month1s)){
    plotsoil2<-subset(soil,(format(soil$dates, "%m")==month1s[j] | format(soil$dates, "%m")==month2s[j]) & format(soil$dates, "%y")==doyear2)
    plotmetout2<-subset(metout,(format(metout$dates, "%m")==month1s[j] | format(metout$dates, "%m")==month2s[j]) & format(metout$dates, "%y")==doyear2)
    with(plotmetout2,plot(TALOC~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="10 cm air"))
    
    with(Ambient_CTemp,points(x~as.POSIXct(Ambient_CTemp$Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    
  } #end months
  dev.off()
  filename<-paste("pred_",doyear,"_3cm_ITEX1U.pdf",sep="")
  
  pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,2)) # set up for 5 plots in 2 columns
  par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
  par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
  
  for(j in 1:length(month1s)){
    plotsoil2<-subset(soil,(format(soil$dates, "%m")==month1s[j] | format(soil$dates, "%m")==month2s[j]) & format(soil$dates, "%y")==doyear2)
    plotmetout2<-subset(metout,(format(metout$dates, "%m")==month1s[j] | format(metout$dates, "%m")==month2s[j]) & format(metout$dates, "%y")==doyear2)
    with(plotsoil2,plot(D1.5cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="1.5cm soil"))
    with(CTemp_3cmBG,points(x~as.POSIXct(Ambient_CTemp$Group.1+3600*offset),cex=0.1,col=addTrans("red",150)))
    
    
  } #end months
  dev.off()
} #end years


years<-seq(2010,2013,1)
month1s<-c("01","03","05","07","09","11")
month2s<-c("02","04","06","08","10","12")
setwd("c:/NicheMapR_Working/projects/snow/camac")
microstation2<-subset(microstation, Site=="ITEX2.0")
Ambient_CTemp<-aggregate(microstation2$Ambient_CTemp,by=list(microstation2$dtime),mean,na.rm=TRUE)
CTemp_3cmBG<-aggregate(microstation2$CTemp_3cmBG,by=list(microstation2$dtime),mean,na.rm=TRUE)
for(i in 1:length(years)){
  doyear<-years[i]
  doyear2<-substr(doyear,3,4)
  
  filename<-paste("pred_",doyear,"_ambient_ITEX2.pdf",sep="")
  
  pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,2)) # set up for 5 plots in 2 columns
  par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
  par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
  
  for(j in 1:length(month1s)){
    plotsoil2<-subset(soil,(format(soil$dates, "%m")==month1s[j] | format(soil$dates, "%m")==month2s[j]) & format(soil$dates, "%y")==doyear2)
    plotmetout2<-subset(metout,(format(metout$dates, "%m")==month1s[j] | format(metout$dates, "%m")==month2s[j]) & format(metout$dates, "%y")==doyear2)
    with(plotmetout2,plot(TALOC~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="10 cm air"))
    #with(microstation2,points(Ambient_CTemp~as.POSIXct(dtime+3600*offset),cex=0.1,col=addTrans("red",150)))
    with(Ambient_CTemp,points(x~as.POSIXct(Ambient_CTemp$Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    
    
  } #end months
  dev.off()
  filename<-paste("pred_",doyear,"_3cm_ITEX2.pdf",sep="")
  
  pdf(filename,paper="A4r",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
  par(mfrow = c(3,2)) # set up for 5 plots in 2 columns
  par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
  par(mar = c(1,2,2,1) + 0.1) # margin spacing stuff
  
  for(j in 1:length(month1s)){
    plotsoil2<-subset(soil,(format(soil$dates, "%m")==month1s[j] | format(soil$dates, "%m")==month2s[j]) & format(soil$dates, "%y")==doyear2)
    plotmetout2<-subset(metout,(format(metout$dates, "%m")==month1s[j] | format(metout$dates, "%m")==month2s[j]) & format(metout$dates, "%y")==doyear2)
    with(plotsoil2,plot(D1.5cm~dates,type='l',col='1',ylab = "Soil Temp (deg C)",xlab = "Year",ylim=c(-10,60),main="1.5cm soil"))
    # with(microstation2,points(CTemp_3cmBG~as.POSIXct(dtime+3600*offset),cex=0.1,col=addTrans("red",150)))
    with(CTemp_3cmBG,points(x~as.POSIXct(Ambient_CTemp$Group.1+3600*offset),cex=0.1,col=addTrans("red",150),type='l'))
    
    
  } #end months
  dev.off()
} #end years


