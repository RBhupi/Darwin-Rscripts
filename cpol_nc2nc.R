#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads CPOL netCDF files (Darwin radar) for each scan
# and writes them back as daily netCDF-4 files with unlimited time axis.
# The data is compressed.
# Warning: file names are expected to follow the pattern some_prefix_string_DATE_TIME.nc
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

#@brief gives list of files belong to the same day as the first file in list_allFiles.
#Method: split the first file's name by "_" and extract 'dates' string,
#then search for the same 'dates' in the list and return them back.
get_1dayFiles <- function (list_allFiles) {
  firstFile <- list_allFiles[1]
  fname_split <- unlist(strsplit(firstFile, "_"))
  date_str <- fname_split[length(fname_split)-1]
  selectFiles <- str_detect(list_allFiles, date_str)
  flist_select <- list_allFiles[selectFiles]
  print(paste(length(flist_select), "file(s) found on the day", date_str))
  return(flist_select)
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

mk_outFilePath <- function (fname, outputDir) {
    fname_split <- unlist(strsplit(fname, "_"))
    pat_replace <- paste("_",fname_split[length(fname_split)], sep="")
    outfName <- str_replace(fname, pattern = pat_replace, ".nc")
    outfPath <- paste(outputDir, basename(outfName), sep="")
    return(outfPath)
}
#------------------------------------------------------------------------------#
indir <- "/home/bhupendra/data/darwin_radar/test/CPOL/"
#make output directory
outdir <- paste(indir, "outNC/", sep="")
dir.create(outdir)

#get all input file names
fpat <- "*.nc"
flist_all <- Sys.glob(paste(indir, fpat, sep=""))
print(paste(length(flist_all), "file(s) in the folder."))

while(length(flist_all)>0) {

  flist <- get_1dayFiles(flist_all)
  flist_all <- flist_all[flist_all!=flist] #remaining files are stored for next iteration

  #open a netcdf file for reading dims
  infile <- nc_open(filename = flist[1])

  #read x y z dims
  x1<- ncvar_get(nc = infile, varid = "x")
  y1<- ncvar_get(nc = infile, varid = "y")
  z1<- ncvar_get(nc = infile, varid = "z")

  #read time from all the files
  time_seconds<- 86400 * laply(flist, ncread_time)

  #also read radar latlon
  radar_lat1 <- ncvar_get(nc = infile, varid = "radar_latitude")
  radar_lon1 <- ncvar_get(nc = infile, varid = "radar_longitude")

  #------------------------- MAKING DIMS & DATA VARIABLES -----------------------#
  #make output dim variables for x, y, z
  x_dim <- ncdim_def(name = "x", units = "km", vals = x1, longname = "distance from radar in zonal direction")
  y_dim <- ncdim_def(name = "y", units = "km", vals = y1, longname = "distance from radar in meridional direction")
  z_dim <- ncdim_def(name = "z", units = "km", vals = z1, longname = "altitude")
  t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                     vals= time_seconds, longname = "Time of the scan", unlim = TRUE)

  #also make radar lat-lon variables
  radar_lat <- ncvar_def(name = "radar_latitude", units = "degrees_north", dim = list(), prec = "float")
  radar_lon <- ncvar_def(name = "radar_longitude", units = "degrees_east", dim = list(), prec = "float")
  radar_latlon_var <- list(radar_lat, radar_lon)

  #This is the list of important "float" variables in the file to read and write
  fvar_names <- c("zh", "vr", "ve", "sw", "zd", "ps", "rs", "kd", "rr", "do", "nw")
  fvar_df <- get_var_longnames(infile, fvar_names)

  #also get long names for "integer" variables
  ivar_names <- c("hc", "cs","cm")
  ivar_df <- get_var_longnames(infile, ivar_names)

  #close input file hear
  nc_close(infile)

  # Now create all float netcdf variables at once using mlply
  out_fvar_list <- mlply(.data = fvar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim, t_dim),
                         prec = "float", compression = 7, chunksizes = c(length(x1), length(y1), 1, 1), missval=-999.0)

  # and also create integer variables the same way and save them in a single list
  out_ivar_list <- mlply(.data = ivar_df, .fun = ncvar_def, units = "", dim = list(x_dim, y_dim, z_dim, t_dim),
                         prec = "integer", compression = 7, chunksizes = c(length(x1), length(y1), 1, 1), missval=-99)

  out_var_list <- append(out_fvar_list, out_ivar_list)

  #------------------------------- WRITING TO THE FILE --------------------------#
  #open output netcdf file
  outfPath <- mk_outFilePath(flist[1], outdir)
  ofile <- nc_create(filename = outfPath, vars =append(out_var_list, radar_latlon_var))

  #write radar lat-lon to the file
  ncvar_put(nc=ofile, varid = "radar_latitude", vals = radar_lat1)
  ncvar_put(nc=ofile, varid = "radar_longitude", vals = radar_lon1)

  print(paste("writing data in", basename(outfPath)))
  pb = txtProgressBar(min =1, max = length(flist), initial = 1, style = 3) #progress bar
  for (f in 1:length(flist)){
    setTxtProgressBar(pb, f)
    infile <- nc_open(flist[f])

    #get the data from input netcdf and put all float variables in output file
    for(var in fvar_names){
      data <- ncvar_get(infile, varid = var)
      ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1, f),
                count = c(length(x1), length(y1), length(z1), 1))
    }

    #get the data from input netcdf and put all integer variables in output file.
    for(var in ivar_names){
      data <- ncvar_get(infile, varid = var)
      ncvar_put(nc = ofile, varid = var, vals = data, start = c(1, 1, 1, f),
                count = c(length(x1), length(y1), length(z1), 1))
    }
    nc_close(infile)
  }
  cat("\n")

  #close file
  nc_close(ofile)
} # while loop ends

print(proc.time()-start_t)
