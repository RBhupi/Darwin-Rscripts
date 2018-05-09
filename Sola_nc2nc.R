#!---------------------------------------------------------------------------
# @author Bhupendra Raut
# @brief This R script reads CPOL netCDF files (Solapur radar) for each scan
# and writes them back as daily netCDF-4 files with unlimited time axis.
# The data is compressed.
#
#@todo :
#==========================================================================================

library(ncdf4)
library(plyr)
library(stringr)

start_t <- proc.time()

#---------------------------------------------------------- FUNCTION DEFINITIONS
#reads time from a single netcdf file. Using with laply.
ncread_time <- function(filename){
    ncfile <- nc_open(filename)
    time <- ncvar_get(ncfile, varid = "time")
    return(time)
}



#reads longname for the vector of varnames and returns a dataframe
get_var_longnames <- function (infile, var_names) {
    var_longNames <- c(NULL)
    for(var in var_names){
        att<- unlist(ncatt_get(infile, varid=var, attname = "long_name"))
        var_longNames <- c(var_longNames, as.vector(att[2]))
    }
    
    # then combin and rename them to use in mlply
    varnames_df<-cbind(var_names, var_longNames)
    colnames(varnames_df) <- c("name", "longname")
    return(varnames_df)
}



#'Function computes WT of the 2d array up to given max_scale.
#'Negative WT are removed. Not tested for non-square data.
#'@param data_matrix
#'@param max_scale
wavelet <- function (data, max_scale){
    #make an output wt array
    dim_data <- dim(data)
    nlat <- dim_data[1]
    nlon <- dim_data[2]
    min_window_size <- 2^(max_scale-1)
    if(nlat<min_window_size || nlon<min_window_size){
        stop(paste("Error: data array is not large enough to comput ",
                  max_scale, "scales!"))
    }
    
    wt <- array(data=0.0, dim = c(max_scale, nlon, nlat))
    
    sf=c(0.0625, 0.25, 0.375) #weighing function
    
    temp1<-array(dim = dim(data))
    temp2<-array(dim = dim(data))
    
    #start Wavelet loop
    for(scale in 1:max_scale){
        x1 <- 2^(scale-1)
        x2 <- 2 * x1
        
        #Row-wise (longitude) smoothing
        for (i in 1:nlon){
            
            #find the indices for prev and next points on the line
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image, "mirror" them
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            if(next1 > nlon | next2 > nlon){
                next1 <- prev1
                next2 <- prev2
            }
            
            
            for (j in 1:(nlat)) {
                left2  <-  data[j, prev2]
                left1  <-  data[j, prev1]
                right1  <-  data[j, next1]
                right2  <-  data[j, next2]
                temp1[j, i]  <-  sf[1] * (left2+right2) +
                    sf[2] * (left1 + right1) + sf[3] * data[j, i]
            }
        }
        
        
        #column-wise (latitude) smoothing
        for(i in 1:nlat){
            
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image use next values
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            if(next1 > nlat | next2 > nlat){
                next1  <-  prev1
                next2 <- prev2
            }
            
            
            
            for(j in 1:nlon){
                top2  <-  temp1[prev2, j]
                top1  <-  temp1[prev1, j]
                bottom1  <-  temp1[next1, j]
                bottom2  <-  temp1[next2, j]
                temp2[i, j]  <-  sf[1] * (top2+bottom2) +
                    sf[2] * (top1 + bottom1) + sf[3] * temp1[i, j]
            }
        }
        
        wt[scale, , ] <- data - temp2
        data <- temp2
    }
    invisible(wt)
}

#'checks if rainFrac is sufficient and returns cumulative sum of all WT components
#' @export
wt_cumsum<-function(data, max_scale){
    min_rainFrac<-0.001
    #check that rainfrac is more than threshold
    rainFrac <- length(data[data>0])/length(data)
    if(rainFrac < min_rainFrac){
        wt <- array(data = 0.0, dim = dim(data))
        invisible(wt)
    }
    
    wt <- wavelet(data, max_scale)
    wt_sum <- apply(wt, MARGIN = c(2, 3), FUN = sum)
    wt_sum <- replace(wt_sum, wt_sum<0, 0.0)
    
    wt_sum <- replace(wt_sum, wt_sum < mean(wt_sum)+3*sd(wt_sum), 0.0)
    invisible(wt_sum)
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
indir <- "/home/bhupendra/data/netcdf_solapur"
setwd(indir)

lev=3
var_id <- "DBZc"
#make output directory
outdir <- paste(indir, "/2dNC/", sep="")
outfPath <- paste(outdir, "SolaRadar2d_lev",lev, "_dbz.nc", sep="")
if(!dir.exists(outdir)) dir.create(outdir)

#get all input file names
fpat <- "*.nc"
flist <- Sys.glob(fpat)
print(paste(length(flist), "file(s) in the folder."))

#open a netcdf file for reading dims
infile <- nc_open(filename = flist[1])

#read x y z dims
x1<- ncvar_get(nc = infile, varid = "x0")
y1<- ncvar_get(nc = infile, varid = "y0")
z1<- ncvar_get(nc = infile, varid = "z0")

#read time from all the files
time_seconds<-laply(flist, ncread_time)

#also read radar latlon
lat1 <- ncvar_get(nc = infile, varid = "lat0")
lon1 <- ncvar_get(nc = infile, varid = "lon0")

#------------------------- MAKING DIMS & DATA VARIABLES -----------------------#
#make output dim variables for x, y, z
x_dim <- ncdim_def(name = "x", units = "km", vals = x1, longname = "distance from radar in zonal direction")
y_dim <- ncdim_def(name = "y", units = "km", vals = y1, longname = "distance from radar in meridional direction")
t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                   vals= time_seconds, longname = "Time of the scan", unlim = TRUE)

#also make radar lat-lon variables
lat2d <- ncvar_def(name = "lat2d", units = "degrees_north", dim = list(x_dim, y_dim), prec = "float")
lon2d <- ncvar_def(name = "lon2d", units = "degrees_east", dim = list(x_dim, y_dim), prec = "float")


#This is the list of important "float" variables in the file to read and write

var_long_name <- get_var_longnames(infile, var_id)

fvar_names <- c("DBZc", "lat2d", "lon2d")
#close input file hear


# Now create netcdf variable
DBZc <- ncvar_def(name=var_id, units = "dBZ", dim = list(x_dim, y_dim, t_dim), prec = "float", compression = 7, 
                  chunksizes = c(length(x1), length(y1), 1), missval=-999.0)

#------------------------------- WRITING TO THE FILE --------------------------#
#open output netcdf file
ofile <- nc_create(filename = outfPath, vars =list(DBZc, lat2d, lon2d))
ncatt_put(nc=ofile, varid=DBZc, attname ="coordinates", attval ="lon2d lat2d")

nc_close(infile)

#write radar lat-lon to the file
ncvar_put(nc=ofile, varid = "lat2d", vals = lat1, start = c(1, 1), count = c(-1, -1))
ncvar_put(nc=ofile, varid = "lon2d", vals = lon1, start = c(1, 1), count = c(-1, -1))





print(paste("writing data in", basename(outfPath)))
pb = txtProgressBar(min =1, max = length(flist), initial = 1, style = 3) #progress bar
for (f in 1:length(flist)){
    setTxtProgressBar(pb, f)
    infile <- nc_open(flist[f])
    
    #get the data from input netcdf and put all float variables in output file
    
    data <- ncvar_get(infile, varid = var_id, start = c(1, 1, lev, 1), count = c(-1, -1, 1, 1))
    data <- replace(data, data<18, 0.0)
    data <- replace(data, is.na(data), 0.0)
    #data<-wt_cumsum(data, 5)
    
    ncvar_put(nc = ofile, varid = var_id, vals = data, start = c(1, 1, f),
              count = c(length(x1), length(y1), 1))
    nc_close(infile)
}
cat("\n")

#close file
nc_close(ofile)


print(proc.time()-start_t)
