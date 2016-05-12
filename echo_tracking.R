#!---------------------------------------------------------------------------
#' @title Tracking Convective Echoes in CPOL Radar Data.
#' @author Bhupendra Raut (www.baraut.info)
#' @description This R script reads netCDF files containing echo heights and rain classification
#' identify objects and estimates flow vectors, using phase correlation with FFT (Leese et al 1971).
#' Then for each object the search region is predicted for next frame and the objects in that
#' region (in frame2) are assign to the object in frame 1.
#' @reference :Leese, John A., Charles S. Novak, and Bruce B. Clark.
#'           "An automated technique for obtaining cloud motion from geosynchronous
#'           satellite data using cross correlation."
#'           Journal of applied meteorology 10.1 (1971): 118-132.

#' @todo
#' 1. Check overlapping pixels to identify better match.
#' 2.
#'
#==========================================================================================
# Start the clock!
start_time <- proc.time()

library(ncdf4)
library(EBImage)  # for bwlabel,
library(spatstat) #for smoothing
library(stringr)
library(clue)  #solve_LSAP()

#----------------------------------------------------------------------fucntions
#' Plots image with objects labels
plot_objects_label<-function(labeled_image, xvalues, yvalues){
    image2D(replace(labeled_image, labeled_image==0, NA), x=xvalues, y=yvalues)
    grid()

    for(object in 1:max(labeled_image)) {
        #center indices of the object assuming it is a rectangle
        obj_id1 <- object

        obj_index <- which(labeled_image==obj_id1, arr.ind = TRUE)

        r1 <- min(obj_index[, 1], na.rm = TRUE)
        r2 <- max(obj_index[, 1], na.rm = TRUE)
        c1 <- min(obj_index[, 2], na.rm = TRUE)
        c2 <- max(obj_index[, 2], na.rm = TRUE)

        obj_centerIndex <- c((r1+r2)/2, (c1+c2)/2)
        text(x=xvalues[obj_centerIndex[1]], y=yvalues[obj_centerIndex[2]],
             toString(object), cex=0.7)
    }
}

#' Returns a single radar scan with echo objects lebeled. smaller objects are removed.
get_filteredFrame <- function(ncfile, scan_num, min_size) {
    echo_height <- get_convHeight(ncfile, scan_num)
    labeled_echo <- bwlabel(echo_height)          #label objects
    frame <-clear_smallEchoes(labeled_echo, min_size)
}

#' Returns a single radar scan with echo objects lebeled. smaller objects are removed.
get_classFrame <- function(ncfile, scan_num) {
    echo_height <- get_convHeight(ncfile, scan_num)
    get_vertical_class(echo_height)
}


#' Reads height and classification data and replaces non-convective amd missing pixels with zero.
get_convHeight <- function(ncfile, scan) {
    dbz_height <- ncvar_get(ncfile, varid = "zero_dbz_cont", start = c(1, 1, scan), count = c(-1, -1, 1))
    steiner <- ncvar_get(ncfile, varid = "steiner_class", start = c(1, 1, scan), count = c(-1, -1, 1))

    dbz_height <-replace(dbz_height, steiner != 2, 0.0)      #set non-convective pixels to zeros
    dbz_height <- replace(dbz_height, is.na(dbz_height), 0.0)     #remove NAs
}

#' Returns labeled image after removing objects smaller than min_size
clear_smallEchoes <- function(label_image, min_size) {
    size_table <- table(label_image[label_image>0]) # remove zero values
    onePix_objects <- as.vector(which(size_table < min_size))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, 0.0)
    }

    label_image <- bwlabel(label_image)
    invisible(label_image)
}

#' Returns vertical classification Cg=1, Cb=2, Co=3
get_vertical_class <- function(conv_height) {
    #min max scan levels for classification
    min_level <- c(5, 15, 31)
    max_level <- c(14, 30, 40)

    for(i in 1:length(min_level)){
        conv_height <- replace(conv_height, conv_height>=min_level[i] &
                                   conv_height<= max_level[i], i)
    }
    return(conv_height) #classified
}


#' Changes base epoch of. Default To_epoch is "1970-01-01.
change_baseEpoch <- function(time_seconds, From_epoch, To_epoch=as.Date("1970-01-01")){
    epoch_diff <- as.integer(From_epoch-To_epoch)
    epoch_diff_seconds <- epoch_diff * 86400 #seconds in a day
    time_newEpoch <- time_seconds + epoch_diff_seconds
    return(time_newEpoch)
}


#' Given two images, the function identifies the matching
#' objects and pair them appropriatly.
get_matchPairs <- function(image1, image2) {
    nObjects1 <- max(image1) #objects in first image
    nObjects2 <- max(image2) #objects in second image

    if(nObjects1==0){ #error if first image is empty
        stop("No echoes found in the first scan.")
    } else if (nObjects2==0){ #all objects will be zero if second image is empty
        zero_pairs <- rep(0, nObjects1)
        return(zero_pairs)
    }

    obj_match <- locate_allObjects(image1, image2)
    pairs <- match_pairs(obj_match) #1-to-1
    return(as.vector(pairs))
}

#' Matches objects into pairs and removes bad matching.
match_pairs <- function(obj_match) {
    pairs <- solve_LSAP(obj_match)
    pairs <- as.vector(pairs)
    ## remove bad matching
    for(pair in 1:length(pairs)){
        if(obj_match[pair, pairs[pair]] >15){
            pairs[pair] <- 0
        }
    }
    return(pairs)
}


#' Matches all the obejects in image1 to the objects in image 2
locate_allObjects <- function(image1, image2) {
    nObjects1 <- max(image1) #objects in first image
    nObjects2 <- max(image2) #objects in second image


    if(nObjects2==0 || nObjects1==0){
        stop("No echoes to track!!!")
    }

    obj_match <- matrix(large_num, nrow = nObjects1, ncol = max(nObjects1, nObjects2))

    ## here we match each object in image1 to all the near-by objects in image2.
    for(obj_id1 in 1:nObjects1) {
        obj1_extent <- get_objExtent(image1, obj_id1) #location and radius
        shift <- get_std_flowVector(obj1_extent, image1, image2, flow_margin, stdFlow_mag)
        #print(paste("fft shift", toString(shift)))

        search_box <- predict_searchExtent(obj1_extent, shift, search_margin)
        search_box <- check_searchBox(search_box, dim(image2)) #search within the image
        obj_found <- find_objects(search_box, image2)  # gives possible candidates
        discrepancy <- get_discrepancy_all(obj_found, image2, search_box, obj1_extent)

        obj_match <- save_objMatch(obj_id1, obj_found, discrepancy, obj_match)

        #print(paste(obj_id1, "==>", toString(obj_found)))
    }

    invisible(obj_match)
}



#' Takes in a labeled image and finds the radius and the center of the given object.
get_objExtent <- function(labeled_image, obj_label) {
    #center indices of the object assuming it is a rectangle
    obj_index <- which(labeled_image==obj_label, arr.ind = TRUE)

    rlength <- max(obj_index[, 1]) - min(obj_index[, 1]) + 1
    clength <- max(obj_index[, 2]) - min(obj_index[, 2]) + 1

    obj_radius<- max(c(rlength, clength))/2 #maximum possible object radius
    obj_center <- c(min(obj_index[, 1])+obj_radius, min(obj_index[, 2]) + obj_radius)
    obj_area <- length(obj_index[, 1])  #size in pixels

    obj_extent<-list(obj_center=obj_center, obj_radius=obj_radius, obj_area=obj_area)
    return(obj_extent)
}


#' Takes in object info (radius and center) and two images to estimate ambient flow.
#' margin is the additional region arround the object used to comput the flow vectors.
get_objAmbientFlow <- function(obj_extent, img1, img2, margin) {
    #coordinates of the flowregion
    r1 <- obj_extent$obj_center[1] - obj_extent$obj_radius - margin
    r2 <- obj_extent$obj_center[1] + obj_extent$obj_radius + margin
    c1 <- obj_extent$obj_center[2] - obj_extent$obj_radius - margin
    c2 <- obj_extent$obj_center[2] + obj_extent$obj_radius + margin

    dims <- dim(img1)
    if(r1<=0 || c1 <=0 || r2>dims[1] || c2 > dims[2]){
        return(c(0, 0))                         #if echo is at the image boundary
    }

    flow_region1 <- img1[r1:r2, c1:c2]
    flow_region2 <- img2[r1:r2, c1:c2]

    return(fft_flowVectors(flow_region1, flow_region2))
}

#' Alternative to get_objAmbientFlow.
#' Flow vectors magnitude is clipped to given magnitude
get_std_flowVector<-function(obj_extent, img1, img2, margin, magnitude){
    shift <- get_objAmbientFlow(obj_extent, img1, img2, margin)
    shift <- replace(shift, shift > magnitude, magnitude)
    shift <- replace(shift, shift < magnitude*-1, magnitude*-1)
    return(shift)
}


#' Estimates flow vectors in two images using cross covariance of the images.
fft_flowVectors <- function (im1, im2) {
    if(max(im1)==0 || max(im2)==0){
        return(c(0, 0))                         #if no object found
    }

    #im1 <- replace(im1, im1==0, runif(1))       # add noise to image background
    #im2 <- replace(im2, im1==0, runif(1))       # This may help when objects are bigger

    crossCov <- fft_crossCov(im1, im2)
    cov_smooth <- blur(as.im(crossCov))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}


#' Computes cross-covariance using FFT, returns shifted covariance image
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
    return(fft_shift(crossCov))
}

#' Rearranges the crossCov matrix so that 'zero' frequency or DC component
#'  is in the middle of the matrix.
#'  This function is adopted from following discussion on stackOverflow
#'  http://stackoverflow.com/questions/30630632/performing-a-phase-correlation-with-fft-in-r
fft_shift <- function(fft_mat) {
    if(class(fft_mat)=='matrix') {
        rd2 <- floor(nrow(fft_mat)/2)
        cd2 <- floor(ncol(fft_mat)/2)

        ## Identify the first, second, third, and fourth quadrants
        q1 <- fft_mat[1:rd2,1:cd2]
        q2 <- fft_mat[1:rd2,(cd2+1):ncol(fft_mat)]
        q3 <- fft_mat[(rd2+1):nrow(fft_mat),(cd2+1):ncol(fft_mat)]
        q4 <- fft_mat[(rd2+1):nrow(fft_mat),1:cd2]

        ## rearrange the quadrants
        centered.t <- rbind(q4,q1)
        centered.b <- rbind(q3,q2)
        centered <- cbind(centered.b,centered.t)

        invisible(Re(centered))
    } else {
        stop("input to fft_shift() should be a matrix")
    }
}


#' Predicts search extent for the object in image2 given shift
predict_searchExtent <- function(obj1_extent, shift, search_radius){
    shifted_center <- obj1_extent$obj_center + shift

    x1 <- shifted_center[1] -search_radius
    x2 <- shifted_center[1] +search_radius
    y1 <- shifted_center[2] -search_radius
    y2 <- shifted_center[2] +search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}




#' Returns NA if search box  outside the image or very small.
check_searchBox <- function(search_box, img_dims){

    if(search_box$x1 <= 0){
        search_box$x1 <- 1
    }
    if(search_box$y1 <= 0){
        search_box$y1 <- 1
    }
    if(search_box$x2 > img_dims[1]){
        search_box$x2 <- img_dims[1]
    }
    if(search_box$y2 > img_dims[2]){
        search_box$y2 <- img_dims[2]
    }

    #search box should be large enough
    if(search_box$x2-search_box$x1 < 5 || search_box$y2-search_box$y1 < 5 ){
        return(NA)
    } else {
        return(search_box)
    }
}

#' Given the search box and image2, returns objects in the region
find_objects <- function(search_box, image2) {
    #if search box is NA then object left the image
    if(is.na(search_box[1])){
        obj_found <- NA
    } else {
        search_area <- image2[search_box$x1:search_box$x2, search_box$y1:search_box$y2]
        obj_found <- unique(as.vector(search_area))
    }
    return(obj_found)
}



#' Returns discrepancies of all the objects found within the search box or NA if
#' no object is present.
get_discrepancy_all <- function(obj_found, image2, search_box, obj1_extent) {
    if(is.na(obj_found[1]) || max(obj_found)==0) {
        obj_id2 <- 0
        dist_pred <- NA
        dist_actual <- NA
        discrepancy <- NA
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0

        if(length(obj_found)==1){ # if this is the only object
            discrepancy <- get_discrepancy(obj_found, image2, search_box, obj1_extent)
            if(discrepancy < 10) discrepancy <- 0 #lower the discrepancy if not too large

        } else { # when more than one objects
            discrepancy <- get_discrepancy(obj_found, image2, search_box, obj1_extent)
        }
    }
    return(discrepancy)
}

#' Saves discrepancy values in obj_match to obj_match array for appropriate objects
save_objMatch <- function(obj_id1, obj_found, discrepancy, obj_match) {
    if(discrepancy >15 || is.na(discrepancy)){
        obj_match[obj_id1, obj_found] <- large_num
    } else {
        obj_match[obj_id1, obj_found] <- discrepancy
    }
    return(obj_match)
}


#' Computes discrepancy for a single object. Check how it is computed.
#' This parameter has most effect on the acccuracy of tracks.
get_discrepancy <- function(obj_found, image2, search_box, obj1_extent) {
    dist_pred <- c(NULL)
    dist_actual <- c(NULL)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(image2, target_obj)
        euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist)

        euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_actual <- append(dist_actual, euc_dist)
        size_changed <- get_ratio(target_extent$obj_area, obj1_extent$obj_area) #change in size

        discrepancy <- dist_pred + size_changed + dist_actual

    }
    return(discrepancy)
}

#' Returns  Euclidean distance between two vectors or matrices
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}

#' Returns ratio (>=1) of bigger number to smaller number when given two number.
get_ratio<-function(x, y){
    if(x>=y)
        return(x/y)
    else
        return(y/x)
}



#' Creates output netcdf file for radar echo tracjecories.
create_outNC <- function(ofile, max_obs) {
    if(file.exists(ofile)){
        print(paste("removing existing file", basename(ofile)))
        file.remove(ofile)
    }


    dim_echo <- ncdim_def("conv_echo", vals=1, units = "", unlim = TRUE,
                          longname = "unique id of convective echo", create_dimvar = TRUE)

    dim_obs <- ncdim_def("obs", vals = seq(max_obs), units="",
                         longname = "observation of a convective echo", create_dimvar = TRUE)

    dim_time <- ncdim_def("time", vals=1, units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan", unlim = TRUE, create_dimvar = TRUE)

    dim_stat <- ncdim_def("stat", vals = seq(3), units="", longname = "1=lived, 2=died, 3=born")

    ## Define Variables
    var_survival <- ncvar_def("survival", units = "", longname = "survival stats for each scan",
                              dim=list(dim_stat, dim_time), missval = -999, prec="integer")

    var_time <- ncvar_def("obs_time", units = "seconds since 1970-01-01 00:00:00 UTC",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer")

    var_xdist <- ncvar_def("x_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float")

    var_ydist <- ncvar_def("y_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float")

    var_x <- ncvar_def("x", units = "", longname = "index along x-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999.0, prec = "integer")

    var_y <- ncvar_def("y", units = "", longname = "index along y-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999.0, prec = "integer")

    var_npix <- ncvar_def("area", units = "pixels", longname = "area of the echo in pixels",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer")

    var_ncg <- ncvar_def("Cg", units = "pixels", longname = "num of Cu Congestus pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer")

    var_ncb <- ncvar_def("Cb", units = "pixels", longname = "num of Cumulonimbus pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer")

    var_nco <- ncvar_def("Co", units = "pixels", longname = "num of Cu overshooting pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer")

    var_list <- list(var_time, var_survival, var_xdist, var_ydist, var_x, var_y,
                     var_npix, var_ncg, var_ncb, var_nco)



    outNC <- nc_create(filename = ofile, vars = var_list)

    #for CF standards
    ncatt_put(outNC, varid = "conv_echo", attname = "cf_role", attval = "trajectory_id")
    ncatt_put(outNC, varid = 0, attname = "featureType", attval = "trajectory")

    description <- paste("The CPOL radar echoes of convective types were separated using Steiner classification scheme and tracked.")

    ncatt_put(outNC, varid = 0, attname = "_description",
              attval = description, prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_creator",
              attval = "Bhupendra Raut", prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_url",
              attval = "www.baraut.info", prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_date_created",
              attval = date(), prec = "text")

    invisible(outNC)
}

#' Writes properties and uids for all objects.
write_update<-function(outNC, current_objects, obj_props, obs_time){
    nobj <- length(current_objects$id1)

    for(object in seq(nobj)){
        nc_start <- c(current_objects$obs_num[object], current_objects$uid[object])
        nc_count <- c(1, 1)

        ncvar_put(outNC, varid = "obs_time", obs_time, start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x", obj_props$x[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y", obj_props$y[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x_dist", obj_props$xdist[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y_dist", obj_props$ydist[object], start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "area", obj_props$area[object],  start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "Cg", obj_props$Cg[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "Cb", obj_props$Cb[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "Co", obj_props$Co[object], start = nc_start, count = nc_count)
    }

}


#' Write survival stats to the file for each scan
write_survival <- function(outNC, survival_stat, time, scan){
    if(!is.atomic(survival_stat)){
        survival_stat <- unlist(survival_stat, use.names = FALSE)
    }

    ncvar_put(outNC, varid = "survival", vals = survival_stat, start = c(1, scan), count = c(3, 1))
    ncvar_put(outNC, varid = "time", vals = time, start = scan-1, count=1)
}




#' Returns a dataframe for objects with ids in frame1 and frame2 and uids (same as ids for first frame).
init_uids <- function(first_frame, pairs){
    nobj <- max(first_frame) #number of objects in frame1
    objects_mat <- matrix(data = NA, ncol = 4, nrow = nobj)

    objects_mat[, 1] <- seq(nobj) #id1
    objects_mat[, 2] <- next_uid(count = nobj) #unique ids
    objects_mat[, 3] <- as.vector(pairs) #as they are in frame2
    objects_mat[, 4] <-rep(1, nobj) #observation number for the echo
    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num")
    return(current_objects)
}


#' Removes dead objects, updates living objects and assign new uids to new born objects.
#' Also, updates number of observations for each echo.
update_current_objects <- function(frame1, pairs, current_objects){
    nobj <- max(frame1)
    objects_mat <- matrix(data = NA, ncol = 4, nrow = nobj)

    objects_mat[, 1] <- seq(nobj) #this is id1

    for (obj in seq(nobj)){
        if(obj %in% current_objects$id2){
            objects_mat[obj, 2] <- current_objects$uid[current_objects$id2==obj]
            objects_mat[obj, 4] <- current_objects$obs_num[current_objects$id2==obj]+1
        } else {
            objects_mat[obj, 2] <- next_uid()
            objects_mat[obj, 4] <- 1 #first observation of the echo
        }
    }

    objects_mat[, 3] <- as.vector(pairs) #as they are in frame2

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num")
    invisible(current_objects)
}


#' Returns sequence of next unique ids and increament the uid_counter.
next_uid<-function(count=1){
    this_uid <- uid_counter + 1:count
    uid_counter <<- uid_counter + count
    return(this_uid)
}

#' Return object's size, location and classification info, xyDist should be a list
get_objectProp <- function(image1, class1, xyDist){
    objprop <- c(NULL)
    nobj <- max(image1)

    for(obj in seq(nobj)){
        obj_index <- which(image1==obj, arr.ind = TRUE)
        objprop$id1 <- append (objprop$id1, obj)  #id in frame1
        objprop$x <- append(objprop$x, floor(median(obj_index[, 2]))) #center column
        objprop$y <- append(objprop$y, floor(median(obj_index[, 1]))) #center row
        objprop$area <- append(objprop$area, length(obj_index[, 1]))

        obj_class <- class1[image1==obj] #class of convection for the object
        #store number of pixels with classification Cg, Cb, Co etc.
        objprop$Cg <- append(objprop$Cg, length(obj_class[obj_class==1]))
        objprop$Cb <- append(objprop$Cb, length(obj_class[obj_class==2]))
        objprop$Co <- append(objprop$Co, length(obj_class[obj_class==3]))
    }

    objprop <- attach_xyDist(objprop, xyDist$x, xyDist$y)
    invisible(objprop)
}


#' Attaches y and x distance from radar in km to object location indices
attach_xyDist<-function(obj_props, xdist, ydist){
    obj_props$xdist <- xdist[obj_props$x]
    obj_props$ydist <- ydist[obj_props$y]
    invisible(obj_props)
}


#' Returns a list with number of objects lived, died and born in this step.
survival_stats <- function(pairs, num_obj2) {
    pairs_vec <- as.vector(pairs)
    obj_lived <- length(pairs_vec[pairs_vec>0])
    obj_died <- length(pairs_vec)-obj_lived
    obj_born <- num_obj2 - obj_lived
    return(list(lived=obj_lived, died=obj_died, born=obj_born))
}
#==============================================================================#

#------------------- Settings for tracking method etc. ------------------------#
search_margin <- 5      #pixels
flow_margin <- 10       #pixels
stdFlow_mag <- 3        #fft_flow will not be faster than this
large_num <- 100000     #a very large number
max_obs<- 100           #longest track that is likely to be recorded
uid_counter <- 0        #(a global variable) start unique id with zero.
min_size <- 2           #objects smaller than this will be filter
#==============================================================================#

#----------------------------------------------------------------Calling Program
setwd("~/data/darwin_radar/2d/")
infile_name <- "./cpol_2D_0506.nc" #a file for a season
outfile_name <- str_replace(infile_name, ".nc", "_tracks.nc")
#outfile_name <- "~/Desktop/test.nc"
print(paste("Opening output file", basename(outfile_name)))
outNC <- create_outNC(outfile_name, max_obs)

#read x, y and time from the file
ncfile <- nc_open(infile_name)
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")

time <- ncvar_get(ncfile, varid="time")
time <- change_baseEpoch(time, From_epoch = as.Date("2004-01-01"))



nscans <- length(time)
newRain <- TRUE         #is this new rainy scan after dry period?

print(paste("Total scans in this file", nscans))
pb = txtProgressBar(min =2, max = nscans, initial = 2, style = 3) #progress bar

# We read the first frame and call it second frame so that in the loop,
# this frame will be copied to frame1 and next frame will be frame2.
frame2 <- get_filteredFrame(ncfile, 1, min_size)
class2 <- get_classFrame(ncfile, 1) #classifictaion

for(scan in 2:nscans){
    setTxtProgressBar(pb, scan) #progress bar

    frame1 <- frame2
    class1 <- class2

    frame2 <- get_filteredFrame(ncfile, scan, min_size)
    class2 <- get_classFrame(ncfile, scan)

    #skip if no echoes in frame 1
    if(max(frame1)==0){         #if no echoes in frame1
        newRain = TRUE          #next rain will be newRain
        write_survival(outNC, survival_stat = rep(0, 3),
                       time = time[scan-1], scan = scan)
        next
    }

    pairs <- get_matchPairs(frame1, frame2)
    obj_props <- get_objectProp(frame1, class1, list(x=x, y=y)) #of frame1

    if(newRain){                #if this is newRain scan, init ids
        current_objects <- init_uids(frame1, pairs) #initiate ids and return
        newRain <- FALSE
    } else {                    #else update old ids
        current_objects <- update_current_objects(frame1, pairs, current_objects)
    }
    write_update(outNC, current_objects, obj_props, time[scan-1]) #for frame1

    #Survival is from frame1 to frame2
    num_obj2 <- max(frame2)
    obj_survival <- survival_stats(pairs, num_obj2)
    write_survival(outNC, survival_stat = obj_survival,
                   time = time[scan-1], scan = scan)
}
cat("\n") #new line required for progress bar

print("closing files")
nc_close(ncfile)
#write unlimited dim and close
ncvar_put(outNC, varid = "conv_echo", vals = seq(uid_counter), start = 1, count = uid_counter)
nc_close(outNC)

# Stop the clock and print the time elapsed
time_elapsed <- (proc.time() - start_time)
print(paste("time elapsed", round(time_elapsed[3]/60), "minutes"))
