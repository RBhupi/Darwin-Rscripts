#' ---
#' Converts track files to 2d_scan files.
#' This program was written for Benjamin to convert netCDF tracks of convection
#' echoes to 2d scans with echo ids.
#' author: "Bhupendra Raut"
#' date: "December 24, 2016"
#' ---
library(ncdf4)
library(stringr)
library(EBImage)  # for bwlabel,
library(spatstat)  #Using for is.empty()

#' returns ids  of the echoes that are alive for the given time
select_track_ids<-function(track_time, scan_time){
    ids <- which(track_time==scan_time, arr.ind = TRUE) #index are ids
    return(ids)
}


#retruns negatively labeled radar scan with convection pixels only
get_RadarScan<-function(ncf_scans, scan_num){
    scan_data <- ncvar_get(ncf_scans, varid="steiner_class",
                           start = c(1, 1, scan_num), count=c(-1, -1, 1))

    #keep only convective pixels
    #scan_data <- replace(scan_data, scan_data==1, 0) #remove non-convective

    scan_data <- replace(scan_data, scan_data!= 2, 0.0) #remove non-convective
    scan_data <- replace(scan_data, is.na(scan_data), 0.0)     #remove NAs

    return(bwlabel_neg(scan_data))
}


#' Lables scan data with object lables using negative integer lables.
#' Puts missing values where appropriate.
bwlabel_neg <- function(scan_data){
    labled_data <- bwlabel(scan_data) #lable continuouse objects

    #we make all the temp lables negative so that they wont crash with echo_lables.
    labled_data<- labled_data *-1

    #we need to add missing values back, as they are relabled by bwlabel()
    labled_data <- replace(labled_data, labled_data==labled_data[1, 1], missing_data)
    labled_data <- replace(labled_data, labled_data==labled_data[61, 61], missing_data)
    return(labled_data)
}

#' Changes base epoch of. Default To_epoch is "1970-01-01.
change_baseEpoch <- function(time_seconds, From_epoch, To_epoch=as.Date("1970-01-01")){
    epoch_diff <- as.integer(From_epoch-To_epoch)
    epoch_diff_seconds <- epoch_diff * 86400 #seconds in a day
    time_newEpoch <- time_seconds + epoch_diff_seconds
    return(time_newEpoch)
}

open_2dncFile <- function(out_fname, infile_nc){
    #read units
    tunit_str <- ncatt_get(ncf_scans, varid ="time", attname = "units")
    xyunit_str <- ncatt_get(ncf_scans, varid ="x", attname = "units")

    #read XY dims
    x1 <- ncvar_get(ncf_scans, varid="x")
    y1<- ncvar_get(ncf_scans, varid="y")

    x_dim <- ncdim_def(name = "x", units = xyunit_str$value, vals = x1,
                       longname = "distance from radar in zonal direction")

    y_dim <- ncdim_def(name = "y", units = xyunit_str$value, vals = y1,
                       longname = "distance from radar in meridional direction")

    t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                       vals= time, longname = "Time of the scan", unlim = TRUE)

    var_scan <- ncvar_def("scan", units = "echo_id", longname = paste("labeled echo id, untracked", untracked_label),
                         dim = list(x_dim, y_dim, t_dim), missval = missing_data, prec = "integer",
                         compression = 9, shuffle = TRUE)

    outNC <- nc_create(filename = out_fname, vars = var_scan)
    return(outNC)
}


write_scan <- function(ncfile, scan_data, scan_num){
    dcount <- dim(scan_data)[1] #awsomeR
    ncvar_put(ncfile, varid = "scan", vals = scan_data, start =c(1, 1, scan_num),
              count = c(dcount, dcount, 1))
}


#' Takes in labeled image and label objects smaller than min_size with id= -1.
label_smallEchoes <- function(label_image, min_size) {
    size_table <- table(label_image[label_image!=0 & label_image!=missing_data])#remove zero/missing values
    onePix_objects <- as.numeric(names(which(size_table < min_size)))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, untracked_label)
    }
    invisible(label_image)
}

#' Labels each object in scan_2d with its given track_ids. Need to give time of the scan.
lable_scan_ids <- function(scan_2d, track_ids, time_of_scan){
    nids <- length(track_ids[, 2])

    for(track in 1:nids){
        echo_xy<-get_echo_xy_and_area(track_ids[track, ], time_of_scan)
        scan_2d <- label_id_inScan(scan_2d, echo_xy, track_ids[track, 2])
    }
    return(scan_2d)
}


#' Returns x-y position of the center of the given echo at given time.
get_echo_xy_and_area <- function(track_id, time){
    x<-ncvar_get(ncf_tracks, varid = "x", start = as.vector(track_id), count=c(1, 1))
    y<-ncvar_get(ncf_tracks, varid = "y", start = as.vector(track_id), count=c(1, 1))
    area <-ncvar_get(ncf_tracks, varid = "area", start = as.vector(track_id), count=c(1, 1))
    return(c(x, y, area))
}


#' This function labels the given echo with the id
label_id_inScan<-function(scan_bwlb, xy, id_label){
    temp_label <- scan_bwlb[xy[1], xy[2]]

    #if center is an empty pixel check surrounding pixels.
    if(temp_label == 0 || temp_label ==missing_data){
        temp_label <- scan_bwlb[(xy[1]-1):(xy[1]+1), (xy[2]-1):(xy[2]+1)]

        #find unique non-zero lables
        temp_label <- unique(temp_label[temp_label!=0 & temp_label!= missing_data])

        #now check the area of each temporary labled echo to match with echo_area.
        for(an_echo in temp_label){
            area_of_echo <- length(which(scan_bwlb==an_echo))
            if(area_of_echo==xy[3]){
                temp_label=an_echo
                scan_bwlb <- replace(scan_bwlb, scan_bwlb==temp_label, id_label)
                break
            }
            # if no echo was found in the scan for the given track-time
            # we simply return the image without labeling any pixels. This will create
            # a problem so I need to make that echo untracked, if it exists somewhere
            # in the scan. These instances are rare (1 in 500 or so), hence they are hard to find and fix.
        }
    }else{
        scan_bwlb <- replace(scan_bwlb, scan_bwlb==temp_label, id_label)
    }

    return(scan_bwlb)
}



################################################################################
untracked_label <- -1
missing_data <- -999

#All File names
fname_scans <- "~/projects/darwin/data/2d/cpol_2D_0607.nc"
#fname_tracks <- "~/Desktop/test_tracks.nc"
fname_tracks <- "~/projects/darwin/data/tracks/cpol_2D_0607_tracks_V17-03.nc"
fname_2dOut <- str_replace(fname_tracks, "_tracks_", "_2d-tracks_")

#---------------------read radar scan file-------------------------#
ncf_scans <- nc_open(fname_scans)
time <- ncvar_get(ncf_scans, varid="time")
time_units <- ncatt_get(ncf_scans, varid = "time", attname = "units")
time <- change_baseEpoch(time, From_epoch =as.Date("2004-01-01"))
nscans <- length(time)

#read XY dims
x1 <- ncvar_get(ncf_scans, varid="x")
y1<- ncvar_get(ncf_scans, varid="y")
missing_val <- ncatt_get(ncf_scans, varid = "steiner_class", attname ="_FillValue")

#----------------------read all track times-----------------------------#
ncf_tracks <- nc_open(fname_tracks, readunlim = FALSE)
track_times <- ncvar_get(ncf_tracks, varid = "record_time", start = c(1, 1), count = c(-1, -1))
min_size <- ncatt_get(ncf_tracks, varid = 0, attname = "min_echoSize_toTrack")

#open output file
outnc <- open_2dncFile(fname_2dOut, ncf_scans)

pb = txtProgressBar(min =1, max = nscans, initial = 1, style = 3) #progress bar

for(scan in seq(nscans)){
    setTxtProgressBar(pb, scan) #advance progress bar

    track_ids <- select_track_ids(track_times, time[scan])
    radar_scan <- get_RadarScan(ncf_scans, scan)
    radar_scan <- label_smallEchoes(radar_scan, min_size$value) #labels 1-pix objects

    #if no tracked echoes, label untracked echoes -1
    if(is.empty(track_ids)){
        write_scan(outnc, radar_scan, scan)
        next
    }

    labled_scan <- lable_scan_ids(radar_scan, track_ids, time[scan])
    write_scan(outnc, labled_scan, scan)
}

nc_close(outnc)


