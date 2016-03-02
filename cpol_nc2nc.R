#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads CPOL netCDF files (Darwin radar) for each scan
# and writes them back as daily netCDF-4 files with unlimited time axis.
# The data is compressed.
#@todo :    1. check all paarameters for correctness,
#           2. add units to variables
#==========================================================================================

library(ncdf4)
library(plyr)
library(stringr)
start_t <- proc.time()

#---------------------------------------------------------- FUNCTION DEFINITIONS
#reads time from a list  of netcdf files
ncread_time_bulk <- function(filelist){
    time_all <- laply(filelist, ncread_time)
    invisible(time_all)
}

#reads time from a single netcdf file
ncread_time <- function(filename){
    ncfile <- nc_open(filename)
    time <- ncvar_get(ncfile, varid = "time")
    return(time)
}
#------------------------------------------------------------------------------#
indir <- "/home/bhupendra/data/darwin_radar/test/CPOL/"
#make output directory
outdir <- paste(indir, "outNC/", sep="")
dir.create(outdir)

#get all input file names
fpat <- "*.nc"
flist_all <- Sys.glob(paste(indir, fpat, sep=""))
#flist <- flist[1:2]
print(paste(length(flist_all), "file(s) in the folder."))

while(length(flist_all)>0) {

    firstFile <- flist_all[1] # this is first file
    fname_split <- unlist(strsplit(firstFile, "_")) #split the fileName
    date_str <- fname_split[length(fname_split)-1] #this is the date
    selectFiles <- str_detect(flist_all, date_str)  #find this pattern in fnames
    flist <- flist_all[selectFiles] # select these files for processing
    flist_all <- flist_all[!selectFiles] #remaining files are stored for next iteration
    print(paste(length(flist), "file(s) on the day", date_str))


    #open netcdf file for reading dims
    infile <- nc_open(filename = flist[1])

    #read x y z dims
    x1<- ncvar_get(nc = infile, varid = "x")
    x_max1 <- ncvar_get(nc= infile, varid = "x_max")
    x_min1 <- ncvar_get(nc= infile, varid = "x_min")
    nx1 <- ncvar_get(nc= infile, varid = "nx")

    y1<- ncvar_get(nc = infile, varid = "y")
    y_max1 <- ncvar_get(nc= infile, varid = "y_max")
    y_min1 <- ncvar_get(nc= infile, varid = "y_min")
    ny1 <- ncvar_get(nc= infile, varid = "ny")

    z1<- ncvar_get(nc = infile, varid = "z")
    z_max1 <- ncvar_get(nc= infile, varid = "z_max")
    z_min1 <- ncvar_get(nc= infile, varid = "z_min")
    nz1 <- ncvar_get(nc= infile, varid = "nz")

    #read time from all teh files
    time_seconds<- 86400 * ncread_time_bulk(flist)

    data_date1 <- ncvar_get(nc = infile, varid = "data_date")
    data_time1 <- ncvar_get(nc = infile, varid = "data_time")

    #also read other domain info
    radar_lat1 <- ncvar_get(nc = infile, varid = "radar_latitude")
    radar_lon1 <- ncvar_get(nc = infile, varid = "radar_longitude")


    #------------------------- MAKING DIMS & DATA VARIABLES -----------------------#
    #make output dim variables for x, y, z
    x_dim <- ncdim_def(name = "x", units = "km", vals = x1, longname = "distance from radar")
    y_dim <- ncdim_def(name = "y", units = "km", vals = y1, longname = "distance from radar")
    z_dim <- ncdim_def(name = "z", units = "km", vals = z1, longname = "altitude")
    t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                       vals= time_seconds, longname = "Time of the scan", unlim = TRUE)


    #This is the list of important "float" variables in the file to read and write
    fvar_names <- c("zh", "vr", "ve", "sw", "zd", "ps", "rs", "kd", "rr", "do", "nw")
    #read their long names
    fvar_longNames <- c(NULL)
    for(var in fvar_names){
        att<- unlist(ncatt_get(infile, varid=var, attname = "long_name"))
        fvar_longNames <- c(fvar_longNames, as.vector(att[2]))
    }

    # then combin and rename them to use in mlply
    fvar_df<-cbind(fvar_names, fvar_longNames)
    colnames(fvar_df) <- c("name", "longname")

    #Repeat the same for "integer" variables
    ivar_names <- c("hc", "cs","cm")
    #read their long names
    ivar_longNames <- c(NULL)
    for(var in ivar_names){
        att<- unlist(ncatt_get(infile, varid=var, attname = "long_name"))
        ivar_longNames <- c(ivar_longNames, as.vector(att[2]))
    }
    #close input file hear
    nc_close(infile)


    # combin and rename
    ivar_df<-cbind(ivar_names, ivar_longNames)
    colnames(ivar_df) <- c("name", "longname")


    #create all float netcdf variables at once using mlply
    out_fvar_list <- mlply(.data = fvar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim, t_dim), prec = "float",
                           compression = 7, chunksizes = c(length(x1), length(y1), 1, 1), missval=-999.0)


    # and also create integer variables the same way
    out_ivar_list <- mlply(.data = ivar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim, t_dim), prec = "integer",
                           compression = 7, chunksizes = c(length(x1), length(y1), 1, 1), missval=-99)

    #------------------------------- WRITING TO THE FILE --------------------------#
    #open output netcdf file
    outfName <- str_replace(flist[1], pattern = "_0000.nc", ".nc")
    outfPath <- paste(outdir, basename(outfName), sep="")
    ofile <- nc_create(filename = outfPath, vars =append(out_fvar_list, out_ivar_list))

    print(paste("writing data in", basename(outfPath)))
    pb = txtProgressBar(min =1, max = length(flist), initial = 1, style = 3) #progress bar
    for (f in 1:length(flist)){
        setTxtProgressBar(pb, f)
        infile <- nc_open(flist[f])
        #get the data from input netcdf and put all float variables in output file
        for(var in fvar_names){
            data <- ncvar_get(infile, varid = var)
            #data <- replace(data, data==NA, -999.0)
            ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1, f),
                      count = c(length(x1), length(y1), length(z1), 1))
        }

        #get the data from input netcdf and put all integer variables in output file.
        for(var in ivar_names){
            data <- ncvar_get(infile, varid = var)
            data <- replace(data, data==NA, -99)
            ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1, f),
                      count = c(length(x1), length(y1), length(z1), 1))
        }
        nc_close(infile)
    }

    print("")
    #close files
    nc_close(ofile)

} # while loop ends

print(proc.time()-start_t)
