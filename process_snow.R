# Code to create netcdf and kmz files from daily output of microclimate model
# snow summaries, here using maximum daily snow, summarizing per two year block
# from 1990-2009.

process_snow <- function(infile, outfile, res=0.05, years=1990:2009, kmz=TRUE, 
                         overwrite=TRUE, ...) {
  # infile:     Path to the csv file containing the data to process.
  # outfile:    Path to the output netcdf file, without extension (also used
  #              for the kmz file if kmz=TRUE).
  # res:        Grid resolution.
  # years:      The years corresponding to the data in columns >=4 of infile.
  # kmz:        Should a kmz file be written?
  # overwrite:  Should existing .nc and .kmz files be overwritten?
  # ...:        Additional arguments accepted by KML (see ?raster::KML). Ignored
  #              if kmz=FALSE.
  require(raster)
  require(ncdf)
  require(lubridate)
  dat <- na.omit(read.csv(infile, header=FALSE))
  y1 <- min(dat[, 3]) - res/2
  y2 <- max(dat[, 3]) + res/2
  x1 <- min(dat[, 2]) - res/2
  x2 <- max(dat[, 2]) + res/2
  nc <- (x2 - x1)/res
  nr <- (y2 - y1)/res
  xy <- dat[, 2:3]
  r <- raster(ncol=nc, nrow=nr, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  cells <- cellFromXY(r, xy)
  tzone <- 'Etc/GMT+10' # doing it this way ignores daylight savings!  
  dates <- seq(ymd(paste0(min(years), '-01-01'), tz=tzone), 
               ymd(paste0(max(years), '-12-31'), tz=tzone), 
               by='days')
  
  # Make a list of empty rasters
  l <- replicate(length(dates), r)
  
  # Fill raster template for each column (day) of data, and assign to those
  # elements of l that aren't 29th Feb.
  l[!(day(dates) == 29 & month(dates) == 2)] <- 
    lapply(dat[, -(1:3)], function(x) {
      r[cells] <- x
      r
    })
  
  # Duplicate Feb 28th snow for Feb 29th
  l[day(dates) == 29 & month(dates) == 2] <- 
    l[which(day(dates) == 29 & month(dates) == 2) - 1]
  
  s <- stack(l)
  writeRaster(s, filename=extension(outfile, 'nc'), overwrite=overwrite)
  snow <- open.ncdf(extension(outfile, 'nc'), write=TRUE, readunlim=TRUE)
  put.var.ncdf(snow, varid='z', vals=dates)
  close.ncdf(snow)
  rm(s, l)
  gc()
  if(kmz) {
    b <- brick(extension(outfile, 'nc'))
    names(b) <- dates
    brks <- c(0, 20, 40, 60, 80, 100, 120, 140, 180, 200, 250, 350, 500)
    KML(x=b, filename=extension(outfile, 'kml'), overwrite=overwrite, ...)  
  }
  invisible(return(NULL))
}

# Example
process_snow('landscape grids/maxsnow_daily.csv', 
             'landscape grids/output/maxsnow_daily_lowres_bio12_1990_2009', 
             kmz=FALSE)
