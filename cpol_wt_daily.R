#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads daily CPOL netCDF files and computes wavelet
# transform to extract cumuliform clouds from the radar data. the output is
# saved in daily netCDF file.
#
#@todo :
#==========================================================================================
library(ncdf4)
library(stringr)
start_t <- proc.time()


#Other Paramters for WT
scale_max <- 3
min_rainFrac <- 0.02
dbz_threshold <- 0.0

#---------------------------------------------------------- FUNCTION DEFINITIONS
#Function computes WT of the 2d array up to given max_scale.
#Negative WT are removed. Not tested for non-square data.
wavelet <- function (data, max_scale){
    #make an output wt array
    dim_data <- dim(data)
    nlat <- dim_data[1]
    nlon <- dim_data[2]
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

#checks if rainFrac is sufficient and returns cumulative sum of all WT components
wt_cumsum<-function(data, max_scale){
    #check that rainfrac is more than threshold
    rainFrac <- length(data[data > dbz_threshold])/length(data)
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

setwd("/home/bhupendra/data/darwin_radar/test/CPOL/outNC/")
flist <- Sys.glob(path = "cpol*.nc")
inFileName <- flist[1]


ncfile <- nc_open(inFileName)

#read x, y, z, and time for reference
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
z <- ncvar_get(ncfile, varid = "z")
nx <- length(x)
ny <- length(y)
nz <- length(z)
nscan<-length(time)

time <- ncvar_get(ncfile, varid = "time")
time_posix <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")

#read bringi reflectivity and make smaller values zero
dbz_bringi <- ncvar_get(ncfile, varid = "zh")
dbz_bringi <- replace(dbz_bringi, dbz_bringi<0, 0.0)
dbz_bringi[is.na(dbz_bringi)] <- 0.0







#------------------------- MAKING DIMS & DATA VARIABLES -----------------------#
outFileName <- str_replace(inFileName, ".nc", "_wt3.nc")

#make output dim variables for x, y, z
x_dim <- ncdim_def(name = "x", units = "km", vals = x, longname = "distance from radar in zonal direction")
y_dim <- ncdim_def(name = "y", units = "km", vals = y, longname = "distance from radar in meridional direction")
z_dim <- ncdim_def(name = "z", units = "km", vals = z, longname = "altitude")
t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                   vals= time, longname = "Time of the scan", unlim = TRUE)

wt_var <- ncvar_def(name = "wt3", units = "dBZ", dim = list(x_dim, y_dim, z_dim, t_dim),
                    prec = "float", compression = 7, chunksizes = c(nx, ny, 1, 1), missval=-999.0)

outNC <- nc_create(filename = outFileName, vars =list(wt_var))

for(scan in 1:nscan){
    for(level in 1:nz){
        wt <- wt_cumsum(dbz_bringi[, , level, scan], max_scale = scale_max)
        ncvar_put(nc = outNC, varid = "wt3", vals = wt, start = c(1, 1, level, scan),
                  count = c(nx, ny, 1, 1))
    }

}

nc_close(outNC)


print(proc.time()-start_t)
