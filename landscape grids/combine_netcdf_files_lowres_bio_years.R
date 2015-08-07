
# code to create netcdf and kmz files from daily output of microclimate model snow summaries, 
# here using maximum daily snow, summarizing per two year block from 1990-2009


library(raster)
library(ncdf)

bioreg<-'11'

setwd('landscape grids/')

alldata<-read.csv('maxsnow_daily.csv',head=FALSE)
alldata<-na.omit(alldata)

res<-0.05
lat1<-min(alldata$V3)-res/2 # min latitude
lat2<-max(alldata$V3)+res/2 # max latitude
lon1<-min(alldata$V2)-res/2 # min longitude
lon2<-max(alldata$V2)+res/2 # max longitude
quadwid<-(lon2-lon1)/res
quadlen<-(lat2-lat1)/res
library('raster')
gridout <- raster(ncol=quadwid, nrow=quadlen, xmn=lon1, xmx=lon2, ymn=lat1, ymx=lat2)


tzone<-paste("Etc/GMT+",10,sep="") # doing it this way ignores daylight savings!
leap<-c(2008,2004,2000,1996,1992)
year<-seq(1990,2009,1)
x<-cbind(alldata[,2],alldata[,3]) # list of co-ordinates
for(j in 1:20){
  ystart<-year[j]
  for(k in 1:3){
    if(k==1){
      for(i in (1+(j-1)*365):(122+(j-1)*365)){
        vals <- cbind(alldata[,2],alldata[,3],alldata[,(i+3)])
        grid <- rasterize(x, gridout, vals[,3])
        #grid<-focal(grid,w=matrix(1,3,3),na.rm=TRUE,fun=mean,pad=TRUE,NAonly=TRUE)
        if(i==(1+(j-1)*365)){s<-grid}else{s<-stack(s,grid)}
        if((i-(j-1)*365)==28 & ystart%in%leap==TRUE){ # add in an extra layer for leap years
          s<-stack(s,grid)
          cat('leap','\n')
        } 
        cat(i,'\n')
      }
      dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart),5,3,tz=tzone)-3600*13, by="days") 
      
    }

  if(k==2){
    for(i in (123+(j-1)*365):(245+(j-1)*365)){
      vals <- cbind(alldata[,2],alldata[,3],alldata[,(i+3)])
      grid <- rasterize(x, gridout, vals[,3])
      #grid<-focal(grid,w=matrix(1,3,3),na.rm=TRUE,fun=mean,pad=TRUE,NAonly=TRUE)
      if(i==(123+(j-1)*365)){s<-grid}else{s<-stack(s,grid)}
      cat(i,'\n')
      dates<-seq(ISOdate(ystart,5,3,tz=tzone)-3600*12, ISOdate((ystart),9,3,tz=tzone)-3600*13, by="days") 
    }
  }
  if(k==3){
    for(i in (246+(j-1)*365):(365+(j-1)*365)){
      vals <- cbind(alldata[,2],alldata[,3],alldata[,(i+3)])
      grid <- rasterize(x, gridout, vals[,3])
      #grid<-focal(grid,w=matrix(1,3,3),na.rm=TRUE,fun=mean,pad=TRUE,NAonly=TRUE)
      if(i==(246+(j-1)*365)){s<-grid}else{s<-stack(s,grid)}
      cat(i,'\n')
      dates<-seq(ISOdate(ystart,9,3,tz=tzone)-3600*12, ISOdate((ystart+1),1,1,tz=tzone)-3600*13, by="days") 
    }
  }
  filename<-paste("output/maxsnow_daily_lowres_bio",bioreg,"_",year[j],"_",k,".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  cat(k,'\n')
  snow<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( snow, varid='z', vals=dates)
  close.ncdf(snow)
  }  
  cat(j,'\n')
}

r1<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1990_1.nc',sep=""))
r2<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1990_2.nc',sep=""))
r3<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1990_3.nc',sep=""))
r4<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1991_1.nc',sep=""))
r5<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1991_2.nc',sep=""))
r6<-brick(paste('output/maxsnow_daily_lowres_bio',bioreg,'_1991_3.nc',sep=""))

rall<-stack(r1,r2,r3,r4,r5,r6)
ystart <- 1990# start year for weather generator calibration dataset or AWAP database
yfinish <- 1991# end year for weather generator calibration dataset
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)
tzone<-paste("Etc/GMT+",10,sep="") # doing it this way ignores daylight savings!

dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((yfinish+1),1,1,tz=tzone)-3600*13, by="days") 

filenam<-paste("maxsnow_daily_lowres_bio",bioreg,"_",ystart,"_",yfinish,".nc",sep="") 
writeRaster(rall, filename=filenam, overwrite=TRUE)
snow<-open.ncdf( filenam, write=TRUE, readunlim=TRUE)
put.var.ncdf( snow, varid='z', vals=dates)
close.ncdf(snow)

a<-brick(paste('maxsnow_daily_lowres_bio',bioreg,"_",ystart,'_',yfinish,'.nc',sep=""))
tzone<-paste("Etc/GMT-",10,sep="")
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12,ISOdate(yfinish,12,31,tz=tzone)-3600*12,by="days")
colours<-terrain.colors(12)
colours[1]<-'#00A60000'
names(a)<-dates
brks1<-c(0,20,40,60,80,100,120,140,180,200,250,350,500)
KML(a,file = paste('snow',ystart,'_',yfinish,sep=''),breaks=brks1,col=colours,overwrite=TRUE)


