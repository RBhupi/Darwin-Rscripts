library(ncdf4)
library(plot3D)
library(animation)
library(RColorBrewer)

#' Returns vertical classification Cg=1, Cb=2, Co=3
get_vertical_class <- function(incFile, nscans) {

    steiner <- ncvar_get(incFile, varid= "steiner_class", start = c(1, 1, 1), count = c(-1, -1, nscans))


    conv_height <- ncvar_get(incFile, varid= "zero_dbz_cont", start = c(1, 1, 1), count = c(-1, -1, nscans))
    conv_height <- replace(conv_height, steiner!=2, 0.0) #remove non-convective


    #min max scan levels for classification
    #names  <- c("Cg", "Cb", "Co")
    min_level <- c(5, 15, 31)
    max_level <- c(14, 30, 40)

    for(i in 1:length(min_level)){
        conv_height <- replace(conv_height, conv_height>=min_level[i] &
                                   conv_height<= max_level[i], i)
    }
    return(conv_height) #classified
}

#' returns scan number for first rainy scan
get_firstRainyScan<-function(vClass){
    dims <- dim(vClass)
    for (i in 1:dims[3]){
        if(any(vClass[, , i]>0, na.rm = TRUE))
            return(i)
    }
    return(NA)
}

#' Changes base epoch of. Default To_epoch is "1970-01-01.
change_baseEpoch <- function(time_seconds, From_epoch, To_epoch=as.Date("1970-01-01")){
    epoch_diff <- as.integer(From_epoch-To_epoch)
    epoch_diff_seconds <- epoch_diff * 86400 #seconds in a day
    time_newEpoch <- time_seconds + epoch_diff_seconds
    return(time_newEpoch)
}

#' returns track ids from the file for the given start/genesis time
get_track_ids<-function(nc_tracks, scan_time){
    track_time <- ncvar_get(nc_tracks, varid = "record_time", start = c(1, 1), count = c(-1, -1))
    ids <- which(track_time==scan_time, arr.ind = TRUE)
    return(ids)
}


#' plots track for given ids
plot_track<-function(incf_tracks, track_ids){
    nids <- length(track_ids[, 2])
    dur <- ncvar_get(incf_tracks, varid = "duration")

    for(track in 1:nids){
        xdist <- ncvar_get(incf_tracks, varid="x_dist", start = c(track_ids[track, 1], track_ids[track, 2]), count=c(1, 1))
        ydist <- ncvar_get(incf_tracks, varid="y_dist", start = c(track_ids[track, 1], track_ids[track, 2]), count=c(1, 1))
        if(dur[track_ids[track, 2]]>1)
            text(xdist, ydist, labels = toString(track_ids[track, 2]), cex = 1.0)
    }
}



#read tracks
ifile_tracks <- "~/Desktop/test_trial.nc"
incf_tracks  <- nc_open(ifile_tracks)

tracks <-  ncvar_get(incf_tracks, varid ="echo_id")
ntracks <- length(tracks)


#read radar data
setwd("~/data/darwin_radar/2d/")
ifile_radar <- "./cpol_2D_0506.nc" #a file for a season
ntime <- 1000
incf_radar <- nc_open(ifile_radar)
x <- ncvar_get(incf_radar, varid="x")
y <- ncvar_get(incf_radar, varid="y")
time <- ncvar_get(incf_radar, varid = "time")
time <- change_baseEpoch(time, From_epoch =as.Date("2004-01-01"))
time_posix <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")


#read the data and get height data from it
vClass <- get_vertical_class(incf_radar, ntime)
scan1 <- get_firstRainyScan(vClass)


colors <- rev(brewer.pal(name = "Set1", n=3))
colors <- c("white", colors)

empty_counter <- 0

saveGIF({
for(scan in 1:ntime) {
    track_ids <- get_track_ids(incf_tracks, time[scan])
    if(empty_counter>2 && is.empty(track_ids)) next
    else if (!is.empty(track_ids))
        empty_counter <-0

    image2D(vClass[, , scan], x=x, y=y, col = colors, breaks = c(-99, 0, 1, 2, 3), NAcol = "grey",
            xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]")
    title(main=strftime(time_posix[scan], tz = "UTC", usetz = TRUE))
    if(is.empty(track_ids)){
        empty_counter <- empty_counter + 1
        next
    }

    plot_track(incf_tracks, track_ids)
}
}, movie.name = "tracks_trial_0506.gif", interval = 0.5)

