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
longlat<-c(145.9, -27.56) # type a long/lat here in decimal degrees
loc <- "Broken Hill, NSW, Australia" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
rungads<-1 # use the Global Aerosol Database?
dailywind<-1 # use daily windspeed database?
terrain<-1 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
snowmodel<-1 # run snow version? (slower!)
ystart<-1993
yfinish<-2008
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

DEP <- c(0., 1.5,  2.5, 5.,  10., 15.,  30., 50.,  100., 200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature


soilpro<-read.csv("c:/git/Tiliqua_rugosa/microclimate/soilprops.txt",header=FALSE)
colnames(soilpro)<-c('i','site','long','lat','desc','blkdens','clay','silt','sand')
soilpro<-subset(soilpro,site==1)
soilpro<-soilpro[,5:9]
soilpro[,2]<-1.4
soilpro[,3]<-90
soilpro[,4]<-0.0
soilpro[,5]<-10
#    
soil_depths<-c(2.5,7.5,22.5,45,80,150)
# plot(soilpro$clay~soil_depths,ylim=c(0,100),col='red',type='l')
# points(soilpro$sand~soil_depths,ylim=c(0,100),col='orange',type='l')
# points(soilpro$silt~soil_depths,ylim=c(0,100),col='grey',type='l')
# title(main=loc)
# legend("topleft", inset=.05,
#        legend=round(soilpro[1,3:5],1),bty="n", 
#        horiz=TRUE, bg=NULL, cex=0.8)

DEP2<-rep(0,18)
j<-1
for(i in 1:length(DEP2)){
  if(i%%2==0){
    DEP2[i]<-DEP2[i-1]+(DEP[j]-DEP2[i-1])/2
  }else{
    DEP2[i]<-DEP[j]
    j<-j+1
  }
}
DEP2<-as.data.frame(floor(DEP2))
colnames(DEP2)<-"DEPTH"
 
par(mfrow=c(2,2))


CampNormTbl9_1<-read.csv('../micro_australia/CampNormTbl9_1.csv')
Richter<-read.csv('../micro_australia/Richter_Table1_SI.csv')
dclay<-0.001 #mm
dsilt<-0.026 #mm
dsand<-1.05 #mm
a<-(soilpro$clay/100)*log(dclay) + (soilpro$sand/100)*log(dsand) + (soilpro$silt/100)*log(dsilt)
b.1<-(((soilpro$clay/100)*log(dclay)^2+(soilpro$sand/100)*log(dsand)^2+(soilpro$silt/100)*log(dsilt)^2)-a^2)^(1/2)
dg<-exp(a)
sigma_g<-exp(b.1)
PES<-(0.5*dg^(-1/2))*-1
b<--2*PES+0.2*sigma_g
PE<-PES*(soilpro$blkdens/1.3)^(0.67*b)
KS<-0.004*(1.3/soilpro$blkdens)^(1.3*b)*exp(-6.9*soilpro$clay/100-3.7*soilpro$silt/100)
BD<-soilpro$blkdens
   

# plot(KS~soil_depths,xlim=c(-1,200),ylim=c(0.000017,0.0058))
KS_spline <-spline(soil_depths,KS,n=201,xmin=0,xmax=200,method='natural')
#points(KS_spline$y,col='red',type='l')
KS_spline<-as.data.frame(cbind(KS_spline$x,KS_spline$y))
colnames(KS_spline)<-c('DEPTH','VALUE')
KS<-merge(DEP2,KS_spline)
KS<-c(KS[1,2],KS[,2])
KS[KS<0.000017]<-0.000017
   
#plot(PE~soil_depths,xlim=c(-1,200),ylim=c(-15,0))
PE_spline <-spline(soil_depths,PE,n=201,xmin=0,xmax=200,method='natural')
#points(PE_spline$y,col='red',type='l')
PE_spline<-as.data.frame(cbind(PE_spline$x,PE_spline$y))
colnames(PE_spline)<-c('DEPTH','VALUE')
PE<-merge(DEP2,PE_spline)
PE<-c(-1*PE[1,2],-1*PE[,2])
   
#plot(b~soil_depths,xlim=c(-1,200),ylim=c(2,24))
b_spline <-spline(soil_depths,b,n=201,xmin=0,xmax=200,method='natural')
points(b_spline$y,col='red',type='l')
b_spline<-as.data.frame(cbind(b_spline$x,b_spline$y))
colnames(b_spline)<-c('DEPTH','VALUE')
b<-merge(DEP2,b_spline)
BB<-c(b[1,2],b[,2])
   
#plot(BD~soil_depths,xlim=c(-1,200),ylim=c(1,1.6))
BD_spline <-spline(soil_depths,BD,n=201,xmin=0,xmax=200,method='natural')
#points(BD_spline$y,col='red',type='l')
BD_spline<-as.data.frame(cbind(BD_spline$x,BD_spline$y))
colnames(BD_spline)<-c('DEPTH','VALUE')
BD<-merge(DEP2,BD_spline)
BD<-c(BD[1,2],BD[,2])

par(mfrow=c(1,1))

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
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 1300 # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 20 # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
rainmult<-1 # rain multiplier for surface soil moisture (use to induce runoff), proportion
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
SoilMoist_Init<-rep(0.3,10) # initial soil water content, m3/m3
evenrain<-1 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
maxpool<-500 # max depth for water pooling on the surface, mm (to account for runoff)
soiltype<-4
CampNormTbl9_1<-read.csv('../micro_australia/CampNormTbl9_1.csv')
fieldcap<-CampNormTbl9_1[soiltype,7] # field capacity, mm
wilting<-CampNormTbl9_1[soiltype,8]  # use value from digital atlas of Australian soils # wilting point, mm
# PE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg 
# KS<-rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
# BB<-rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameterPE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg 
# PE[1:2]<-CampNormTbl9_1[7,4] #air entry potential J/kg 
# KS[1:2]<-CampNormTbl9_1[7,6] #saturated conductivity, kg s/m3
# BB[1:2]<-CampNormTbl9_1[7,5] #soil 'b' parameter
# PE[10:19]<-CampNormTbl9_1[soiltype,4] #air entry potential J/kg 
# KS[10:19]<-CampNormTbl9_1[soiltype,6] #saturated conductivity, kg s/m3
# BB[10:19]<-CampNormTbl9_1[soiltype,5] #soil 'b' parameter
BD<-rep(1.4,19)# Mg/m3, soil bulk density for soil moisture calcs
L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
LAI<-0.0 # leaf area index, used to partition traspiration/evaporation from PET
REFL<-0.10 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-0. # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise            											
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-50. # maximum available shade (%)
runshade<-0. # run the model twice, once for each shade level (1) or just for the first shade level (0)?
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
Usrhyt <- 1# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.325 # snow density (mg/m3)
snowmelt<-1 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1. # undercatch multipier for converting rainfall to snow
rainmelt<-0.013 # paramter in equation that melts snow with rainfall as a function of air temp
write_input<-1 # write csv files of final input to working directory? 1=yes, 0=no.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value



siteslist<-read.csv('model tests/newsites.txt')

for(j in 1:nrow(siteslist)){
   longlat<-c(siteslist[j,2],siteslist[j,3])
   site<-siteslist[j,1]

 # run the model
  maindir<-getwd()
  setwd('/git/micro_australia/')
  niche<-list(L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
  source('NicheMapR_Setup_micro.R')
  nicheout<-NicheMapR(niche)
  setwd(maindir)

# get output
    dim<-nicheout$dim
    metout<-as.data.frame(nicheout$metout[1:(dim*24),]) # above ground microclimatic conditions, min shade
    shadmet<-as.data.frame(nicheout$shadmet[1:(dim*24),]) # above ground microclimatic conditions, max shade
    soil<-as.data.frame(nicheout$soil[1:(dim*24),]) # soil temperatures, minimum shade
    shadsoil<-as.data.frame(nicheout$shadsoil[1:(dim*24),]) # soil temperatures, maximum shade
    soilmoist<-as.data.frame(nicheout$soilmoist[1:(dim*24),]) # soil water content, minimum shade
    shadmoist<-as.data.frame(nicheout$shadmoist[1:(dim*24),]) # soil water content, maximum shade
    humid<-as.data.frame(nicheout$humid[1:(dim*24),]) # soil humidity, minimum shade
    shadhumid<-as.data.frame(nicheout$shadhumid[1:(dim*24),]) # soil humidity, maximum shade
    soilpot<-as.data.frame(nicheout$soilpot[1:(dim*24),]) # soil water potential, minimum shade
    shadpot<-as.data.frame(nicheout$shadpot[1:(dim*24),]) # soil water potential, maximum shade
    rainfall<-as.data.frame(nicheout$RAINFALL)
    MAXSHADES<-as.data.frame(nicheout$MAXSHADES)
    elev<-as.numeric(nicheout$ALTT)
    REFL<-as.numeric(nicheout$REFL)
    longlat<-as.matrix(nicheout$longlat)
    ectoin<-rbind(elev,REFL,longlat,0,0,1990,1990+nyears-1)

#     write.csv(metout,paste(microdir,'metout',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(soil,paste(microdir,'soil',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(soilpot,paste(microdir,'soilpot',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(humid,paste(microdir,'humid',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(soilmoist,paste(microdir,'soilmoist',ystart,'_',yfinish,'.csv',sep=""))
#     if(runshade==0){
#       write.csv(metout,paste(microdir,'shadmet',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(soil,paste(microdir,'shadsoil',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(humid,paste(microdir,'shadhumid',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(soilpot,paste(microdir,'shadpot',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(soilmoist,paste(microdir,'shadmoist',ystart,'_',yfinish,'.csv',sep=""))
#     }else{
#       write.csv(shadmet,paste(microdir,'shadmet',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(shadsoil,paste(microdir,'shadsoil',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(shadhumid,paste(microdir,'shadhumid',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(shadpot,paste(microdir,'shadpot',ystart,'_',yfinish,'.csv',sep=""))
#       write.csv(shadmoist,paste(microdir,'shadmoist',ystart,'_',yfinish,'.csv',sep=""))
#     }
#     write.csv(rainfall,paste(microdir,'rainfall',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(ectoin,paste(microdir,'ectoin',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(DEP,paste(microdir,'DEP',ystart,'_',yfinish,'.csv',sep=""))
#     write.csv(MAXSHADES,paste(microdir,'MAXSHADES',ystart,'_',yfinish,'.csv',sep=""))

tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
soil<-cbind(dates,soil)
metout<-cbind(dates,metout)
shadsoil<-cbind(dates,shadsoil)
shadmet<-cbind(dates,shadmet)
soilmoist<-cbind(dates,soilmoist)
shadmoist<-cbind(dates,shadmoist)
humid<-cbind(dates,humid)
shadhumid<-cbind(dates,shadhumid)
soilpot<-cbind(dates,soilpot)
shadpot<-cbind(dates,shadpot)

dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
rainfall<-as.data.frame(cbind(dates2,rainfall))
colnames(rainfall)<-c("dates","rainfall")

years_snow<-c('93','94','95','96','97','98','99','00','01','02','03','04','05','06','07','08')
metout<-as.data.frame(cbind(dates,metout))
soil<-as.data.frame(cbind(dates,soil))
shadmet<-as.data.frame(cbind(dates,shadmet))
shadsoil<-as.data.frame(cbind(dates,shadsoil))


 if(nyears==1){
  plotmetout<-metout
  plotsoil<-soil
  #with(plotmetout,plot(TS*10~dates,type='l',col='grey',ylim=c(-100,550),ylab = "Snow Depth (cm)",xlab = "Year"))
  #with(plotmetout,plot(T2*10~dates,type='l',col='cyan',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
  #with(plotsoil,points(D100cm*10~dates,type='l',col='dark green',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
  #with(plotsoil,points(D20cm*10~dates,type='l',col='dark blue',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
  with(plotmetout,plot(SNOWDEP~dates,type='l',ylim=c(-10,350),ylab = ""))

    snowobs<-read.csv(file=paste("c:/NicheMapR_Working/projects/snow/new sites/",site,".csv",sep=''))

  Date<-ISOdatetime(snowobs$Year,snowobs$Month,snowobs$Day,hour = 12, min = 0, sec = 0, tz = "GMT")
  snowobs<-as.data.frame(cbind(Date,snowobs$Depth))
  colnames(snowobs)<-c('date','depth')
  points(snowobs$depth~snowobs$date,type='p',col='red',cex=.5)
  snowfall<-subset(plotmetout,SNOWFALL>0)
  #points(snowfall$SNOWFALL*10~snowfall$dates,type='l',col='blue',cex=.5)
}

#site<-"Lake Mountain"
for(i in 1:16){
  if(i<=8){
  filename<-paste("model tests/",site,"1993_2000.pdf",sep="")
  }else{
  filename<-paste("model tests/",site,"2000_2008.pdf",sep="")
  }
  if(i==1){
    pdf(filename,paper="A4",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
    par(mfrow = c(4,2)) # set up for 4 plots in 2 columns
    par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
    par(mar = c(1,2,1,1) + 0.1) # margin spacing stuff  
  }
  if(i==9){
    pdf(filename,paper="A4",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
    par(mfrow = c(4,2)) # set up for 4 plots in 2 columns
    par(oma = c(2,2,2,0) + 0.1) # margin spacing stuff
    par(mar = c(1,2,1,1) + 0.1) # margin spacing stuff
  }
plotmetout<-subset(metout,format(metout$dates, "%y")==years_snow[i])
plotsoil<-subset(soil,format(soil$dates, "%y")==years_snow[i])
#with(plotmetout,plot(TS*10~dates,type='l',col='grey',ylim=c(-100,550),ylab = "Snow Depth (cm)",xlab = "Year"))
#with(plotmetout,plot(T2*10~dates,type='l',col='cyan',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
#with(plotsoil,points(D100cm*10~dates,type='l',col='dark green',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
#with(plotsoil,points(D20cm*10~dates,type='l',col='dark blue',ylim=c(-100,350),ylab = "Snow Depth (cm)",xlab = "Year"))
with(plotmetout,plot(SNOWDEP~dates,type='l',ylim=c(-10,350),ylab = ""))
snowobs<-read.csv(file=paste("c:/NicheMapR_Working/projects/snow/new sites/",site,".csv",sep=''))
# if(site=="Mount Buller"){
# snowobs<-read.csv(file="c:/NicheMapR_Working/projects/snow/buller.csv")
# }
# if(site=="Falls Creek"){
# snowobs<-read.csv(file="c:/NicheMapR_Working/projects/snow/fallscreek.csv")
# }
# if(site=="Mount Hotham"){
# snowobs<-read.csv(file="c:/NicheMapR_Working/projects/snow/hotham.csv")
# }
Date<-ISOdatetime(snowobs$Year,snowobs$Month,snowobs$Day,hour = 12, min = 0, sec = 0, tz = "GMT")
snowobs<-as.data.frame(cbind(Date,snowobs$Depth))
colnames(snowobs)<-c('date','depth')
points(snowobs$depth~snowobs$date,type='p',col='red',cex=.5)
snowfall<-subset(plotmetout,SNOWFALL>0)
#points(snowfall$SNOWFALL*10~snowfall$dates,type='l',col='blue',cex=.5)
yearplot<-1993+i-1
legend(min(plotmetout$dates), 340, yearplot,box.lwd = 0,box.col = "transparent",bg = "transparent")
if(i==8 | i==16){
  title(main=site,outer=T)
  dev.off()
}
}
} # end loop after plots
# summarize snow predictions

metsoil<-cbind(metout,soil[,3:12])
metsoil<-as.data.frame(cbind(dates,metsoil))

maxsnow_daily<-aggregate(metsoil$SNOWDEP,by=list((substr(metsoil$dates,1,10))),max)
maxsnow_annual<-aggregate(metsoil$SNOWDEP,by=list((substr(metsoil$dates,1,4))),max)
maxtsurf_daily<-aggregate(metsoil$D0cm,by=list((substr(metsoil$dates,1,10))),max)
maxtsurf_annual<-aggregate(metsoil$D0cm,by=list((substr(metsoil$dates,1,4))),max)
mintsurf_daily<-aggregate(metsoil$D0cm,by=list((substr(metsoil$dates,1,10))),min)
mintsurf_annual<-aggregate(metsoil$D0cm,by=list((substr(metsoil$dates,1,4))),min)
maxt5cm_daily<-aggregate(metsoil$D5cm,by=list((substr(metsoil$dates,1,10))),max)
maxt5cm_annual<-aggregate(metsoil$D5cm,by=list((substr(metsoil$dates,1,4))),max)
mint5cm_daily<-aggregate(metsoil$D5cm,by=list((substr(metsoil$dates,1,10))),min)
mint5cm_annual<-aggregate(metsoil$D5cm,by=list((substr(metsoil$dates,1,4))),min)

count <- function(x) {
  length(na.omit(x))
}
mod_thresh<-10.
metsoil_snow<-subset(metsoil,SNOWDEP>mod_thresh)

  meansnow_annual<-aggregate(metsoil_snow$SNOWDEP,by=list((substr(metsoil_snow$dates,1,4))),mean)
dummysnowdates<-as.data.frame(cbind(seq(1993,2008,1),rep(0,16)))
colnames(dummysnowdates)<-c('Group.1','x')
d<-merge(dummysnowdates,meansnow_annual,by="Group.1",all=TRUE)
meansnow_annual<-cbind(d$Group.1,d$x.y)
colnames(meansnow_annual)<-c('Group.1','x')
meansnow_annual[is.na(meansnow_annual)]<-0
meansnow_annual<-as.data.frame(meansnow_annual)
  snowdays<-aggregate(metsoil_snow$SNOWDEP,by=list((substr(metsoil_snow$dates,1,10))),mean)
  snowdays$x[snowdays$x<=mod_thresh]<-0 # selection of the number of hours of frost per day at which frost is observed
  snowdays<-subset(snowdays,x>0)
  snowdays_daily<-snowdays
  
  snowdays<-aggregate(snowdays$x,by=list((substr(snowdays$Group,1,4))),count)
  
  #frostdays_twodeg<-aggregate(frostdays_twodeg$x,by=list((substr(frostdays_twodeg$Group,1,4))),sum)
  dummysnowdates<-as.data.frame(cbind(seq(1993,2008,1),rep(0,16)))
  colnames(dummysnowdates)<-c('Group.1','x')
  label1<-paste("hours",snowdays[,1],sep="")
  label2<-paste("days",snowdays[,1],sep="")
  label3<-paste("first",snowdays[,1],sep="")
  label4<-paste("last",snowdays[,1],sep="")
  label5<-paste("period",snowdays[,1],sep="")

    firstsnow<-as.data.frame(aggregate(metsoil_snow$dates,by=list((substr(metsoil_snow$dates,1,4))),min))
    firstsnow[,2]<-format(as.Date(firstsnow[,2]), "%d/%m/%Y")
    firstsnow[,2]<-strptime(firstsnow$x, "%d/%m/%Y")$yday+1
    lastsnow<-as.data.frame(aggregate(metsoil_snow$dates,by=list((substr(metsoil_snow$dates,1,4))),max))
    lastsnow[,2]<-format(as.Date(lastsnow[,2]), "%d/%m/%Y") 
    lastsnow[,2]<-strptime(lastsnow$x, "%d/%m/%Y")$yday+1
    firstsnow<-merge(firstsnow,dummysnowdates,by='Group.1',all.y = TRUE)
    lastsnow<-merge(lastsnow,dummysnowdates,by='Group.1',all.y = TRUE)
    firstsnow<-firstsnow[,-3]
    lastsnow<-lastsnow[,-3]
    snowperiod<-lastsnow[,2]-firstsnow[,2]
    snowdays<-merge(snowdays,dummysnowdates,by='Group.1',all.y = TRUE)
    firstsnow[is.na(firstsnow)]<-0
    lastsnow[is.na(lastsnow)]<-0
    snowperiod[is.na(snowperiod)]<-0
    snowdays[is.na(snowdays)]<-0
    for(p in 1:length(snowperiod)){ #make sure that frost period is 1 for when there is only one day of frost
      if(snowdays[p,2]==1){
        snowperiod[p]<-1
      }
    }

# summarize snow observations

snowobs$date<-as.Date(Date, "%d/%m/%Y")
snowobs[is.na(snowobs)] <- 0
maxsnow_daily_obs<-aggregate(snowobs$depth,by=list((substr(snowobs$date,1,10))),max)
maxsnow_annual_obs<-aggregate(snowobs$depth,by=list((substr(snowobs$date,1,4))),max)


snowobs_snow<-subset(snowobs,depth>mod_thresh)

meansnow_annual_obs<-aggregate(snowobs_snow$depth,by=list((substr(snowobs_snow$date,1,4))),mean)

snowdays_obs<-aggregate(snowobs_snow$depth,by=list((substr(snowobs_snow$date,1,10))),mean)
snowdays_obs$x[snowdays_obs$x<=mod_thresh]<-0 # selection of the number of hours of frost per day at which frost is observed
snowdays_obs<-subset(snowdays_obs,x>0)
snowdays_daily_obs<-snowdays_obs

snowdays_obs<-aggregate(snowdays_obs$x,by=list((substr(snowdays_obs$Group,1,4))),count)

#frostdays_twodeg<-aggregate(frostdays_twodeg$x,by=list((substr(frostdays_twodeg$Group,1,4))),sum)
dummysnowdates<-as.data.frame(cbind(seq(1993,2008,1),rep(0,16)))
colnames(dummysnowdates)<-c('Group.1','x')
label1<-paste("hours",snowdays[,1],sep="")
label2<-paste("days",snowdays[,1],sep="")
label3<-paste("first",snowdays[,1],sep="")
label4<-paste("last",snowdays[,1],sep="")
label5<-paste("period",snowdays[,1],sep="")

firstsnow_obs<-as.data.frame(aggregate(snowobs_snow$date,by=list((substr(snowobs_snow$date,1,4))),min))
firstsnow_obs[,2]<-format(as.Date(firstsnow_obs[,2]), "%d/%m/%Y")
firstsnow_obs[,2]<-strptime(firstsnow_obs$x, "%d/%m/%Y")$yday+1
lastsnow_obs<-as.data.frame(aggregate(snowobs_snow$date,by=list((substr(snowobs_snow$date,1,4))),max))
lastsnow_obs[,2]<-format(as.Date(lastsnow_obs[,2]), "%d/%m/%Y") 
lastsnow_obs[,2]<-strptime(lastsnow_obs$x, "%d/%m/%Y")$yday+1
firstsnow_obs<-merge(firstsnow_obs,dummysnowdates,by='Group.1',all.y = TRUE)
lastsnow_obs<-merge(lastsnow_obs,dummysnowdates,by='Group.1',all.y = TRUE)
firstsnow_obs<-firstsnow_obs[,-3]
lastsnow_obs<-lastsnow_obs[,-3]
snowperiod_obs<-lastsnow_obs[,2]-firstsnow_obs[,2]
snowdays_obs<-merge(snowdays_obs,dummysnowdates,by='Group.1',all.y = TRUE)
firstsnow_obs[is.na(firstsnow_obs)]<-0
lastsnow_obs[is.na(lastsnow_obs)]<-0
snowperiod_obs[is.na(snowperiod_obs)]<-0
snowdays_obs[is.na(snowdays_obs)]<-0
for(p in 1:length(snowperiod_obs)){ #make sure that frost period is 1 for when there is only one day of frost
  if(snowdays_obs[p,2]==1){
    snowperiod_obs[p]<-1
  }
}



filename<-paste("model tests/",site,"stats.pdf",sep="")
pdf(filename,paper="A4",width=8,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
par(mfrow = c(3,2))
par(oma = c(2,2,2,0) + 0.1)
par(mar = c(2,5,2,2) + 0.1)
plot(maxsnow_annual_obs$x~maxsnow_annual$x)
abline(1,1)
plot(meansnow_annual_obs$x~meansnow_annual$x)
abline(1,1)
plot(snowdays_obs$x.x~snowdays$x.x)
abline(1,1)
plot(firstsnow_obs$x~firstsnow$x)
abline(1,1)
plot(lastsnow_obs$x~lastsnow$x)
abline(1,1)
plot(snowperiod_obs~snowperiod)
abline(1,1)
dev.off()

maxsnow<-lm(maxsnow_annual_obs$x~maxsnow_annual$x)
R2_maxsnow<-summary(maxsnow)$r.squared
rmsd_maxsnow<-sqrt(mean(((maxsnow_annual_obs$x-maxsnow_annual$x)^2),na.rm=TRUE))
rmsdp_maxsnow<-rmsd_maxsnow/(max(rbind(maxsnow_annual_obs$x,maxsnow_annual$x),na.rm=TRUE)-min(rbind(maxsnow_annual_obs$x,maxsnow_annual$x),na.rm=TRUE))

meansnow<-lm(meansnow_annual_obs$x~meansnow_annual$x)
R2_meansnow<-summary(meansnow)$r.squared
rmsd_meansnow<-sqrt(mean(((meansnow_annual_obs$x-meansnow_annual$x)^2),na.rm=TRUE))
rmsdp_meansnow<-rmsd_meansnow/(max(rbind(meansnow_annual_obs$x,meansnow_annual$x),na.rm=TRUE)-min(rbind(meansnow_annual_obs$x,meansnow_annual$x),na.rm=TRUE))

snowdays_mod<-lm(snowdays_obs$x.x~snowdays$x.x)
R2_snowdays<-summary(snowdays_mod)$r.squared
rmsd_snowdays<-sqrt(mean(((snowdays_obs$x.x-snowdays$x.x)^2),na.rm=TRUE))
rmsdp_snowdays<-rmsd_snowdays/(max(rbind(snowdays_obs$x.x,snowdays$x.x),na.rm=TRUE)-min(rbind(snowdays_obs$x.x,snowdays$x.x),na.rm=TRUE))

snowperiod_mod<-lm(snowperiod_obs~snowperiod)
R2_snowperiod<-summary(snowperiod_mod)$r.squared
rmsd_snowperiod<-sqrt(mean(((snowperiod_obs-snowperiod)^2),na.rm=TRUE))
rmsdp_snowperiod<-rmsd_snowperiod/(max(rbind(snowperiod_obs,snowperiod),na.rm=TRUE)-min(rbind(snowperiod_obs,snowperiod),na.rm=TRUE))


firstsnow_mod<-lm(firstsnow_obs$x~firstsnow$x)
R2_firstsnow<-summary(firstsnow_mod)$r.squared
rmsd_firstsnow<-sqrt(mean(((firstsnow_obs$x-firstsnow$x)^2),na.rm=TRUE))
rmsdp_firstsnow<-rmsd_firstsnow/(max(rbind(firstsnow_obs$x,firstsnow$x),na.rm=TRUE)-min(rbind(firstsnow_obs$x,firstsnow$x),na.rm=TRUE))

lastsnow_mod<-lm(lastsnow_obs$x~lastsnow$x)
R2_lastsnow<-summary(lastsnow_mod)$r.squared
rmsd_lastsnow<-sqrt(mean(((lastsnow_obs$x-lastsnow$x)^2),na.rm=TRUE))
rmsdp_lastsnow<-rmsd_lastsnow/(max(rbind(lastsnow_obs$x,lastsnow$x),na.rm=TRUE)-min(rbind(lastsnow_obs$x,lastsnow$x),na.rm=TRUE))

filename<-paste("model tests/",site,"stats.csv",sep="")
stats_sum1<-rbind(R2_maxsnow,R2_meansnow,R2_snowdays,R2_firstsnow,R2_lastsnow,R2_snowperiod)
stats_sum2<-rbind(rmsd_maxsnow,rmsd_meansnow,rmsd_snowdays,rmsd_firstsnow,rmsd_lastsnow,rmsd_snowperiod)
stats_sum3<-rbind(rmsdp_maxsnow,rmsdp_meansnow,rmsdp_snowdays,rmsdp_firstsnow,rmsdp_lastsnow,rmsdp_snowperiod)
stats_sum_titles<-c('maxsnow','meansnow','snowdays','firstsnow','lastsnow','snowperiod')
stats_sum_cols<-c('obs','r2','rmsd','rmsdp')
stats<-cbind(stats_sum_titles,stats_sum1,stats_sum2,stats_sum3)
colnames(stats)<-stats_sum_cols
write.csv(stats,filename)
R2_maxsnow
R2_meansnow
R2_snowdays
R2_firstsnow
R2_lastsnow
R2_snowperiod


} # end loop through alps sites