
library(ncdf4)
library(plyr)

#get input file names
indir <- "~/Dropbox/Bhupendra_shared/Scripts/RScripts/darwin_project/cpol/"
fpat <- "*.nc"
flist <- Sys.glob(paste(indir, fpat, sep=""))
print(paste("found", length(flist), "file(s)."))

#make output directory
outdir <- paste(indir, "outNC/", sep="")
dir.create(outdir)

#open netcdf file for reading
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

#read time and transform to POSIXct time unit "seconds since 1970-01-01"
time1<- ncvar_get(nc = infile, varid = "time")
time1 <- as.POSIXct(time1*86400, origin = "1970-01-01", tz="UTC")

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

# combin and rename
ivar_df<-cbind(ivar_names, ivar_longNames)
colnames(ivar_df) <- c("name", "longname")


#create all float netcdf variables at once using mlply
out_fvar_list <- mlply(.data = fvar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim), prec = "float",
                      compression = 7, chunksizes = c(length(x1), length(y1), 1), missval=-999.0)


# and also create integer variables the same way
out_ivar_list <- mlply(.data = ivar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim), prec = "integer",
                      compression = 7, chunksizes = c(length(x1), length(y1), 1), missval=-99)

#------------------------------- WRITING TO THE FILE --------------------------#
#open output netcdf file
outfName <- paste(outdir, basename(flist[1]), sep="")
ofile <- nc_create(filename = outfName, vars =append(out_fvar_list, out_ivar_list))

#get the data from input netcdf and put all float variables in output file
for(var in fvar_names){
  data <- ncvar_get(infile, varid = var)
  #data <- replace(data, data==NA, -999.0)
  ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1),
            count = c(length(x1), length(y1), length(z1)))
}

#get the data from input netcdf and put all integer variables in output file
for(var in ivar_names){
    data <- ncvar_get(infile, varid = var)
    data <- replace(data, data==NA, -99)
    ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1),
              count = c(length(x1), length(y1), length(z1)))
}

#close files
nc_close(ofile)
nc_close(infile)
