#' ---
#' title: "Tracking Convection Echoes in CPOL Radar Data"
#' author: "Bhupendra Raut"
#' date: "August 10, 2016"
#' ---

#'
#' This R script reads netCDF files containing maximum echo heights for 0 dBZ
#' and rain classification from steiner method. I obtained these files from Benjamin Moebis.
#' We first remove non-convective pixels and identify the contiguous convective regions.
#' The flow vectors are estimated for the region of convection using phase correlation method (Leese et al 1971).
#' Using these flow vector as the first guess,  the search region is predicted for each object.
#' In the next frame, all objects in that region (in frame2) are assigned to the object in frame 1.
#' The disparity/cost function is computed for each combination of pairs and assignment is done using
#' Hungarian method. For the object pairs with large disparity, we assumed that the old object is dead in frame1
#' and the new object was born in frame2.
#'
#' The unique identity numbers are produced for new objects and was kept constant through out the life of that object.
#' No merging or spliting ia taken care of in this version of the code.
#'
#'
#'
#' ToDo
#' 1. Check overlapping pixels to identify better match.
#' 2. Use vertical profile of clouds for better assignment
#' 3. change size_change factor
#+ echo=FALSE
#==========================================================================================
# Start the clock!
start_time <- proc.time()

#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
#' Following R packages are required.
library(ncdf4)    #Read/Write netcdf-4 files
library(EBImage)  # for bwlabel,
library(spatstat) #for smoothing in fft functions
library(stringr)  #string manipulations
library(clue)     #solve_LSAP() assignment problem

#+ echo=FALSE
#----------------------------------------------------------------------fucntions

#' Plots image with objects labels. This is used to test images when processed 1-by-1.
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

#' Returns a single radar scan from the netcdf file.
#' Smaller objects are removed and the rest are lebeled.
get_filteredFrame <- function(ncfile, scan_num, min_size) {
    echo_height <- get_convHeight(ncfile, scan_num)
    labeled_echo <- bwlabel(echo_height)          #label objects
    frame <-clear_smallEchoes(labeled_echo, min_size)
    invisible(frame)
}

#' Returns a single radar scan of classification for given scan_num.
#' Convective objects lebeled with vertical class 1=Cg, 2=Cb, 3=Co
get_classFrame <- function(ncfile, scan_num) {
    echo_height <- get_convHeight(ncfile, scan_num)
    get_vertical_class(echo_height) #returns this
}


#' Reads height and classification data and replaces non-convective and missing pixels with zero.
get_convHeight <- function(ncfile, scan) {
    dbz_height <- ncvar_get(ncfile, varid = "zero_dbz_cont",
                            start = c(1, 1, scan), count = c(-1, -1, 1))

    steiner <- ncvar_get(ncfile, varid = "steiner_class",
                         start = c(1, 1, scan), count = c(-1, -1, 1))

    dbz_height <-replace(dbz_height, steiner != 2, 0.0)      #set non-convective pixels to zeros
    dbz_height <- replace(dbz_height, is.na(dbz_height), 0.0)     #remove NAs
}

#' Takes in labeled image removes objects smaller than min_size and returns re-labeled image.
clear_smallEchoes <- function(label_image, min_size) {
    size_table <- table(label_image[label_image>0]) # remove zero values
    onePix_objects <- as.vector(which(size_table < min_size))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, 0.0)
    }

    label_image <- bwlabel(label_image)
    invisible(label_image)
}

#' Given the convective height image, it returns vertical classification (Cg=1, Cb=2, Co=3)
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


#' Changes base epoch of 'time in seconds'. Default To_epoch is "1970-01-01.
#' From_epoch should be in the format as.Date("2004-01-01").
change_baseEpoch <- function(time_seconds, From_epoch, To_epoch=as.Date("1970-01-01")){
    epoch_diff <- as.integer(From_epoch-To_epoch)
    epoch_diff_seconds <- epoch_diff * 86400 #seconds in a day
    time_newEpoch <- time_seconds + epoch_diff_seconds
    return(time_newEpoch)
}


#' Given two images, the function identifies the matching
#' objects and pair them appropriatly. See disparity function.
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
#' The bad matching is when disparity is more than the set value.
match_pairs <- function(obj_match) {
    pairs <- solve_LSAP(obj_match)
    pairs <- as.vector(pairs)
    ## remove bad matching
    for(pair in 1:length(pairs)){
        if(obj_match[pair, pairs[pair]] > max_desparity){
            pairs[pair] <- 0
        }
    }
    return(pairs)
}


#' Matches all the obejects in image1 to the objects in image 2.
#' This is the main function to be called on two sets of radar images, for tracking.
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

        if(exists("current_objects"))
            shift <- correct_shift(shift, current_objects, obj_id1)

        search_box <- predict_searchExtent(obj1_extent, shift, search_margin)
        search_box <- check_searchBox(search_box, dim(image2)) #search within the image
        obj_found <- find_objects(search_box, image2)  # gives possible candidates
        disparity <- get_disparity_all(obj_found, image2, search_box, obj1_extent)

        obj_match <- save_objMatch(obj_id1, obj_found, disparity, obj_match)
    }

    invisible(obj_match)
}


#' takes in flow vector based shift and current_object dataframe which has last
#' headings, and check if they are resonably close if not rejects or modify shift and return.
#' Note:  frame2 of last timestep is now frame1, but current_objects still has it as frame2.
#' So id2 in the last frame2 are actually ids related to frame1 now.
correct_shift<-function(this_shift, current_objects, object_id1){
    last_heads <- c(current_objects$xhead[current_objects$id2==object_id1],
    current_objects$yhead[current_objects$id2==object_id1])

    #for small shifts and empty last shifts
    if(is.na(last_heads) || all(last_heads<=1 && last_heads>=-1, na.rm = TRUE))
        return(this_shift)
    else if(any(abs(this_shift-last_heads)>4))  #if they are too different
        return(last_heads)                      #then trust last_heads
    else return((this_shift+last_heads)/2) #else retun the average of both
}


#' Takes in a labeled image and finds the radius, area and the center of the given object.
get_objExtent <- function(labeled_image, obj_label) {
    #center indices of the object assuming it is a rectangle
    obj_index <- which(labeled_image==obj_label, arr.ind = TRUE)

    xlength <- max(obj_index[, 1]) - min(obj_index[, 1]) + 1
    ylength <- max(obj_index[, 2]) - min(obj_index[, 2]) + 1

    obj_radius<- max(c(xlength, ylength))/2 #maximum possible object radius
    obj_center <- c(min(obj_index[, 1])+obj_radius, min(obj_index[, 2]) + obj_radius)
    obj_area <- length(obj_index[, 1])  #size in pixels

    obj_extent<-list(obj_center=obj_center, obj_radius=obj_radius,
                     obj_area=obj_area, obj_index=obj_index)
    return(obj_extent)
}


#' Returns object_extent with number of pixel of each class for the given object.
get_objClass_extent <- function(label_image, class_image, obj_label){
    objExtent <- get_objExtent(label_image, obj_label)
    objClass <- get_object_vertProfile(label_image, class_image, obj_label)
    objExtent$Cg <- objClass$Cg
    objExtent$Cb <- objClass$Cb
    objExtent$Co <- objClass$Co
    return(objExtent)
}



#' Takes in object info (radius and center) and two images to estimate ambient flow.
#' Margin is the additional region arround the object used to comput the flow vectors.
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
    magnitude_negative <- magnitude * -1
    shift <- replace(shift, shift < magnitude_negative, magnitude_negative)
    return(shift)
}


#' Estimates flow vectors in two images using cross covariance of the images.
#'
#' Leese, John A., Charles S. Novak, and Bruce B. Clark.
#'           "An automated technique for obtaining cloud motion from geosynchronous
#'           satellite data using cross correlation."
#'           Journal of applied meteorology 10.1 (1971): 118-132.
fft_flowVectors <- function (im1, im2) {
    if(max(im1)==0 || max(im2)==0){
        return(c(0, 0))                         #if no object found
    }

    #im1 <- replace(im1, im1==0, runif(1))       # add noise to image background
    #im2 <- replace(im2, im1==0, runif(1))       # when objects are big and smooth

    crossCov <- fft_crossCov(im1, im2)
    cov_smooth <- blur(as.im(crossCov))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}


#' Computes cross-covariance using FFT method, returns shifted covariance image
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
    return(fft_shift(crossCov))
}


#' Rearranges the crossCov matrix so that 'zero' frequency or DC component
#'  is in the middle of the matrix. Taken from stackoverflow Que. 30630632
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


#' Predicts search extent/region for the object in image2 given the image shift.
predict_searchExtent <- function(obj1_extent, shift, search_radius){
    shifted_center <- obj1_extent$obj_center + shift

    x1 <- shifted_center[1] -search_radius
    x2 <- shifted_center[1] +search_radius
    y1 <- shifted_center[2] -search_radius
    y2 <- shifted_center[2] +search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}


#' Returns NA if search box  outside the image or search box is very small.
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
    if(search_box$x2-search_box$x1 < 4 || search_box$y2-search_box$y1 < 4 ){
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


#' Returns disparities of all the objects found within the search box or NA if
#' no object is present.
get_disparity_all <- function(obj_found, image2, search_box, obj1_extent) {


    if(is.na(obj_found[1]) || max(obj_found)==0) {
        obj_id2 <- 0
        dist_pred <- NA
        dist_actual <- NA
        disparity <- NA
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0

        if(length(obj_found)==1){ # if this is the only object
            disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
            if(disparity <= 2) disparity <- 0 #lower the disparity if not too large

        } else { # when more than one objects to match
            disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
        }
    }
    return(disparity)
}


#' If disparity is large then it saves a large number for the value to reduce
#' the chances of this pairing to zero, else it save the value in the obj_match array.
save_objMatch <- function(obj_id1, obj_found, disparity, obj_match) {
    if(disparity > max_desparity || is.na(disparity)){
        obj_match[obj_id1, obj_found] <- large_num
    } else {
        obj_match[obj_id1, obj_found] <- disparity
    }
    return(obj_match)
}


#' Computes disparity for a single object. Check how it is computed for detail.
#' This parameter has most effect on the acccuracy of the tracks.
#'
get_disparity <- function(obj_found, image2, search_box, obj1_extent) {
    dist_pred <- c(NULL)
    dist_actual <- c(NULL)
    change <- c(NULL)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(image2, target_obj)

        euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist)

        euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_actual <- append(dist_actual, euc_dist)
        size_changed <- get_sizeChange(target_extent$obj_area, obj1_extent$obj_area) #change in size
        change <- append(change, size_changed)

        #overlap <- obj_overlap(target_extent$obj_index, obj1_extent$obj_index)
    }

    #This is crucial parameter that affect the results
    disparity <- dist_pred + change # + dist_actual
    return(disparity)
}


#' Returns  Euclidean distance between two vectors or matrices.
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}


#' Returns change in size of the eacho as ratio of bigger number to smaller
#' number when given two number, minus 1.
get_sizeChange<-function(x, y){
    if(x < 5 && y <5 ) # if too small, return zero
        return(0)
    else if(x>=y)
        return(x/y - 1)
    else
        return(y/x - 1)
}




#' Creates output netcdf file for radar echo tracjecories. This is the longest function.
create_outNC <- function(ofile, max_obs) {
    if(file.exists(ofile)){
        print(paste("removing existing file", basename(ofile)))
        file.remove(ofile)
    }
    deflat <- 9

    dim_echo <- ncdim_def("echo_id", vals=1, units = "", unlim = TRUE,
                          longname = "unique id of convective echo", create_dimvar = TRUE)

    dim_obs <- ncdim_def("records", vals = seq(max_obs), units="",
                         longname = "observation records")

    dim_time <- ncdim_def("time", vals=1, units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan", unlim = TRUE, create_dimvar = TRUE)

    dim_stat <- ncdim_def("stat", vals = seq(3), units="", longname = "survival stats; lived, died, born")

    ## Define Variables
    var_survival <- ncvar_def("survival", units = "", longname = "survival stats for each scan",
                              dim=list(dim_stat, dim_time), missval = -999, prec="integer",
                              compression = deflat, shuffle = TRUE)

    var_dur <- ncvar_def("duration", units = "", longname = "duration of echo in time-steps",
                         dim=dim_echo, missval = -999, prec="integer",
                         compression = deflat, shuffle = TRUE)

    var_parent <- ncvar_def("parent", units="", longname = "id of the parent echo",
                            dim=dim_echo, missval = 0, prec = "integer")

    var_time <- ncvar_def("record_time", units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan for each record",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                          compression = deflat, shuffle = TRUE)

    var_xdist <- ncvar_def("x_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float",
                           compression = deflat)

    var_ydist <- ncvar_def("y_dist", units = "Km", longname = "distance from Radar",
                           dim = list(dim_obs, dim_echo), missval = -999.0, prec = "float",
                           compression = deflat)


    var_x <- ncvar_def("x", units = "", longname = "index along x-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                       compression = deflat, shuffle = TRUE)

    var_y <- ncvar_def("y", units = "", longname = "index along y-coordinate",
                       dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                       compression = deflat, shuffle = TRUE)

    var_npix <- ncvar_def("area", units = "pixels", longname = "area of the echo in pixels",
                          dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                          compression = deflat, shuffle = TRUE)

    var_ncg <- ncvar_def("Cg", units = "pixels", longname = "num of Cu Congestus pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                         compression = deflat, shuffle = TRUE)

    var_ncb <- ncvar_def("Cb", units = "pixels", longname = "num of Cumulonimbus pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                         compression = deflat, shuffle = TRUE)

    var_nco <- ncvar_def("Co", units = "pixels", longname = "num of Cu overshooting pixels",
                         dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
                         compression = deflat, shuffle = TRUE)

    var_list <- list(var_time, var_survival, var_dur, var_parent, var_xdist, var_ydist,
                     var_x, var_y, var_npix, var_ncg, var_ncb, var_nco)


    outNC <- nc_create(filename = ofile, vars = var_list)

    #for CF standards
    ncatt_put(outNC, varid = "echo_id", attname = "cf_role", attval = "trajectory_id")
    ncatt_put(outNC, varid = 0, attname = "featureType", attval = "trajectory")

    description <- paste("The CPOL radar echoes of convective types were separated using Steiner classification scheme and tracked. Added 'parent' variable.")

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


#' Writes properties and uids for all objects into output netcdf file.
write_update<-function(outNC, current_objects, obj_props, obs_time){
    nobj <- length(current_objects$id1)

    for(object in seq(nobj)){
        nc_start <- c(current_objects$obs_num[object], current_objects$uid[object])
        nc_count <- c(1, 1)

        ncvar_put(outNC, varid = "record_time", obs_time, start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x", obj_props$x[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y", obj_props$y[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x_dist", obj_props$xdist[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y_dist", obj_props$ydist[object], start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "area", obj_props$area[object],  start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "Cg", obj_props$Cg[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "Cb", obj_props$Cb[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "Co", obj_props$Co[object], start = nc_start, count = nc_count)
    }

    write_duration(outNC, current_objects)
}


#' Writes number of observations for dead objects.
#' The time written here is one step earlier (needs correction).
write_duration <- function(outNC, current_objects){
    nobj <- length(current_objects$id1)

    for (obj in seq(nobj)){
        if(current_objects$id2[obj]==0){
            ncvar_put(outNC, varid = "duration", current_objects$obs_num[obj],
                      start=current_objects$uid[obj], count=1)
            ncvar_put(outNC, varid = "parent", current_objects$parent[obj],
                      start=current_objects$uid[obj], count=1)
        }
    }
}


#' Write survival stats (number of lived, dead and born objects) to the file for each scan.
write_survival <- function(outNC, survival_stat, time, scan){
    if(!is.atomic(survival_stat)){
        survival_stat <- unlist(survival_stat, use.names = FALSE)
    }

    ncvar_put(outNC, varid = "survival", vals = survival_stat, start = c(1, scan), count = c(3, 1))
    ncvar_put(outNC, varid = "time", vals = time, start = scan-1, count=1)
}


#' Returns a dataframe for objects with unique ids and their corresponding ids in frame1 and frame2.
#' This function is called when new rainy scan is seen after the period of no rain or the first time.
init_uids <- function(first_frame, second_frame, pairs){
    nobj <- max(first_frame) #number of objects in frame1
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj)               #id1
    objects_mat[, 2] <- next_uid(count = nobj) #unique ids
    objects_mat[, 3] <- as.vector(pairs)#as they are in frame2
    objects_mat[, 4] <-rep(1, nobj)     #observation number for the echo
    objects_mat[, 5] <-rep(0, nobj)



    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "parent")
    current_objects <- attach_xyheads(first_frame, second_frame, current_objects)
    return(current_objects)
}


#' Attaches last xyheads to current objects for future use.
attach_xyheads <- function(frame1, frame2, current_objects) {

    nobj <- length(current_objects$uid)
    xhead <- yhead <- NULL

    for (obj in seq(nobj)) {
        if(current_objects$id1[obj]>0 && current_objects$id2[obj]>0){
            center1<- get_objectCenter(current_objects$id1[obj], frame1)
            center2<- get_objectCenter(current_objects$id2[obj], frame2)
            xhead <- append(xhead, center1[1]-center2[1])
            yhead <- append(yhead, center1[2]-center2[2])
        }else{ #if object is dead write NA
            xhead<-append(xhead, NA)
            yhead<-append(yhead, NA)
        }
    }
        return(cbind(current_objects, xhead, yhead)) #attach values
}


#' Returns index of center pixel of the given object id from a labeled image.
get_objectCenter<-function(obj_id, labeled_image){
    obj_index <- which(labeled_image==obj_id, arr.ind = TRUE)
    center_x <- floor(median(obj_index[, 1])) #center column
    center_y <- floor(median(obj_index[, 2])) #center row
    return(c(center_x, center_y))
}


#'
#' Removes dead objects, updates living objects and assign new uids to new born objects.
#' Also, updates number of valid observations for each echo.
#' This function is called when rain continues from the last frame.
#' This is a complicated function to understand.
#'
#' See how the pairs vector looks like for a real case. The pairs
#' shows mapping of the current frame1 and frame2. This shows that frame2 has 4 objects.
#' The objects [1, 2, 3, 4] in current frame2 are mapped with objects [0, 1, 2, 3]
#' in current frame1. Thus, object 1 in frame2 is new born. Others can be traced back to
#' frame1.
#'
#' pairs>>
#'
#' 0, 1, 2, 3
#'
#' Now check old_objects and remember that at this instant, id2 (in the old_objects)
#' correspond to the objects in current frame1 which was frame2 in the earlier
#' time-step, and that they are the same frame.
#'
#' old_objects>>
#'
#' id1, uid, id2, obs_num, xhead, yhead
#'
#'  1, 1, 1,   2,   1,  0
#'
#'  2,  12,  3, 2, 0,  -1
#'
#' So the object 1 and 3 in current frame1 (earlier it was frame2 with id2) existed
#' before and has "uid" (11 and 12). We will copy their "uid" to our object_matrix
#' and increament the observation number (obs_num).
#' For object 2 and 4 in current frame2 which do not exist in frame1,
#' we will ask for new uids. This information will be written  in
#' current_objects and return for writting in to the output file.
update_current_objects <- function(frame1, frame2, pairs, old_objects){
    nobj <- max(frame1)
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj) # this is id1 at current step

    for (obj in seq(nobj)){
        if(obj %in% old_objects$id2){ # but same was id2 in the last step
            # so they should get same uid as last time
            objects_mat[obj, 2] <- old_objects$uid[old_objects$id2==obj]
            objects_mat[obj, 4] <- old_objects$obs_num[old_objects$id2==obj] + 1
            objects_mat[obj, 5] <- old_objects$parent[old_objects$id2==obj]
        } else {
            objects_mat[obj, 2] <- next_uid()
            objects_mat[obj, 4] <- 1 #first observation of the echo
            objects_mat[obj, 5] <- get_parent_uid(obj, frame1, old_objects)
        }
    }

    objects_mat[, 3] <- as.vector(pairs) #match as they are in frame2

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "parent")
    current_objects <- attach_xyheads(frame1, frame2, current_objects)
    invisible(current_objects)
}


#' returns unique id of the parent (or zero) for given object in frame1.
#' Also remember that old object id2 is actual id1 in frame1, as we still have
#' to update the object_ids.
get_parent_uid<-function(obj, frame1, old_objects){
    parent_id <- find_parent(obj, frame1)
    if (parent_id==0) return(0)

    parent_index <- which(old_objects$id2==parent_id)

    # If it is first observation of the object then it will not be recorded in
    # old_objects$id2, ans it will not be suitable as the parent.
    if(!(parent_id %in% old_objects$id2)) return(0)

    parent_uid <- old_objects$uid[parent_index]
    return(parent_uid)
}


#' This function checks near by objects in the frame for the given new-born object.
find_parent <- function(id1_newObj, frame1){
    if(max(frame1)==1) return(0) # If there is only one object, then dont look for parent

    #get length and indices of the given object pixels
    object_ind <- which(frame1==id1_newObj, arr.ind = TRUE)
    object_size <- length(object_ind[,1])

    #Do this for all other objects in the frame
    neighbour_ind <- which(frame1>0 & frame1!=id1_newObj, arr.ind = T)
    neighbour_size <- length(neighbour_ind[, 1])

    #make empty vectors
    neighbour_dist <- NULL
    neighbour_id <- NULL
    size_ratio <- NULL
    size_diff <- NULL

    # We are chekcing for all object pixels and finding the nearest pixel.
    for(pix in seq(object_size)){
        for(neighbour in seq(neighbour_size)){
            euc_dist <- euclidean_dist(as.vector(object_ind[pix, ]),
                                       as.vector(neighbour_ind[neighbour, ]))
            neighbour_dist <- append(neighbour_dist, euc_dist)

            pix_id <- as.vector(neighbour_ind[neighbour, ])
            neighbour_id <- append(neighbour_id, frame1[pix_id[1], pix_id[2]])
        }
    }

    nearest_object_id <- neighbour_id[which(neighbour_dist<4)]
    the_nearest_object <- neighbour_id[which(neighbour_dist==min(neighbour_dist))]

    if (is.empty(nearest_object_id)) #if no close neighbour return 0
        return(0)

    # This is to take care of multiple objects in the neighbouring region.
    neigh_objects <- unique(nearest_object_id)
    for(object in neigh_objects){
        nearest_object_size <- length(frame1[frame1==object])
        size_ratio <- append(size_ratio, nearest_object_size/object_size)
        size_diff <- append(size_diff, nearest_object_size - object_size)
    }

    # id of the object which has max size_ratio
    big_ratio_obj <- neigh_objects[which(size_ratio==max(size_ratio))]
    big_diff_obj  <- neigh_objects[which(size_diff==max(size_diff))]

    #if both are same call it the parent
    if(big_ratio_obj==big_diff_obj)
        return(big_diff_obj[1])
    else
        return(big_diff_obj[1])
    # NOTE: 1. At this time we are calling big_diff_obj as parent in all the situations.
    # This looks like good a first guess. But if needed we can make it more
    # complex and use ratio and size_diff as cost function.
    # 2. We are not considering the possibility of multiple potential parents
    # beyond this point.
}



#' Returns sequence of next unique ids and increament the uid_counter.
next_uid<-function(count=1){
    this_uid <- uid_counter + 1:count
    uid_counter <<- uid_counter + count
    return(this_uid)
}


#' Return all the object's size, location and classification info,
#' xyDist should be a list of x_dist and y_dist in km.
get_objectProp <- function(image1, class1, xyDist){
    objprop <- c(NULL)
    nobj <- max(image1)

    for(obj in seq(nobj)){
        obj_index <- which(image1==obj, arr.ind = TRUE)
        objprop$id1 <- append (objprop$id1, obj)  #id in frame1
        objprop$x <- append(objprop$x, floor(median(obj_index[, 1]))) #center column
        objprop$y <- append(objprop$y, floor(median(obj_index[, 2]))) #center row
        objprop$area <- append(objprop$area, length(obj_index[, 1]))

        obj_class <- get_object_vertProfile(image1, class1, obj_label = obj) #class of convection for the object

        #store number of pixels with classification Cg, Cb, Co etc.
        objprop$Cg <- append(objprop$Cg, obj_class$Cg)
        objprop$Cb <- append(objprop$Cb, obj_class$Cb)
        objprop$Co <- append(objprop$Co, obj_class$Co)
    }
    objprop <- attach_xyDist(objprop, xyDist$x, xyDist$y)
    invisible(objprop)
}


#' Returns number of Cg, Cb and Co type pixels in the given object.
get_object_vertProfile <- function(label_image, class_image, obj_label){
    obj_class <- class_image[label_image==obj_label]
    nCg <- length(obj_class[obj_class==1])
    nCb <- length(obj_class[obj_class==2])
    nCo <- length(obj_class[obj_class==3])
    return(list(Cg = nCg, Cb = nCb, Co = nCo))
}


#' Attaches y and x distance from radar in km to object location indices
attach_xyDist<-function(obj_props, xdist, ydist){
    obj_props$xdist <- xdist[obj_props$x]
    obj_props$ydist <- ydist[obj_props$y]
    invisible(obj_props)
}


#' Returns a list with number of objects lived, died and born between this step and the next one.
survival_stats <- function(pairs, num_obj2) {
    pairs_vec <- as.vector(pairs)
    obj_lived <- length(pairs_vec[pairs_vec>0])
    obj_died <- length(pairs_vec)-obj_lived
    obj_born <- num_obj2 - obj_lived
    return(list(lived=obj_lived, died=obj_died, born=obj_born))
}
#+ echo=FALSE
#==============================================================================#

#+ echo=TRUE
#'----------------------- Settings for tracking method ------------------------#
search_margin <- 4          #pixels
flow_margin <- 3            #pixels
stdFlow_mag <- 5            #fft_flow will not be faster than this
min_signif_movement <- 2    #not used at this time
large_num <- 1000           #a large number
max_obs<- 60                #longest recoreded track (eles show error).
min_size <- 2               #objects smaller than this will be filter
max_desparity <- 15         # two objects with more desparity than this value, are not same.
#==============================================================================#

#'Here we call the above functions to track convective echoes in the subsequent images.
#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE

setwd("~/data/darwin_radar/2d/")
file_list <- Sys.glob(paths = "./cpol_2D_????.nc")

for(infile_name in file_list){
    outfile_name <- str_replace(infile_name, ".nc", "_tracks_V16_9.nc")
    #outfile_name <- paste("../tracks/", basename(outfile_name), sep="")
    print(paste("Opening output file", outfile_name))
    outNC <- create_outNC(outfile_name, max_obs)
    uid_counter <- 0            #(a global variable) to start uid with 1, count from 0.

    #read x, y and time from the file
    ncfile <- nc_open(infile_name)
    x <- ncvar_get(ncfile, varid = "x")
    y <- ncvar_get(ncfile, varid = "y")


    time <- ncvar_get(ncfile, varid="time")
    time <- change_baseEpoch(time, From_epoch = as.Date("2004-01-01"))


    start_scan <- 1
    end_scan <- length(time)

    newRain <- TRUE         #is this new rainy scan after dry period?

    print(paste("Total scans in this file", end_scan-start_scan+1))
    pb = txtProgressBar(min =start_scan, max = end_scan, initial = 1, style = 3) #progress bar

    #' To continuousely track the echoes in the images,
    #' we read the first frame and save it in to frame2 then in the loop,
    #' this frame will be copied to frame1 and next frame will be frame2.
    #'
    #' frame1 <-- frame2 <-- new frame  (repeat)
    #+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
    frame2 <- get_filteredFrame(ncfile, start_scan, min_size)
    class2 <- get_classFrame(ncfile, start_scan) #classifictaion

    for(scan in (start_scan+1):end_scan){
        setTxtProgressBar(pb, scan) #progress bar

        frame1 <- frame2
        class1 <- class2
        frame2 <- get_filteredFrame(ncfile, scan, min_size)
        class2 <- get_classFrame(ncfile, scan)

        if(scan==end_scan){ #if this is the last scan make it zero.
            frame2 <- replace(frame2, frame2>0, 0)
        }


        #skip if no echoes in frame 1
        if(max(frame1, na.rm = TRUE)==0){         #if no echoes in frame1
            newRain = TRUE          #next rain will be newRain
            write_survival(outNC, survival_stat = rep(0, 3),
                           time = time[scan-1], scan = scan)

            if(exists("current_objects"))
                rm(current_objects)
            next
        }

        pairs <- get_matchPairs(frame1, frame2)
        obj_props <- get_objectProp(frame1, class1, list(x=x, y=y)) #of frame1

        if(newRain){                #if this is newRain scan, init ids
            current_objects <- init_uids(frame1, frame2, pairs) #initiate ids and return
            newRain <- FALSE
        } else {                    #else update old ids
            current_objects <- update_current_objects(frame1, frame2, pairs, current_objects)
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
    ncvar_put(outNC, varid = "echo_id", vals = seq(uid_counter), start = 1, count = uid_counter)
    nc_close(outNC)

    # Stop the clock and print the time elapsed
    time_elapsed <- (proc.time() - start_time)
    print(paste("time elapsed", round(time_elapsed[3]/60), "minutes"))

}#file loop end


#+ echo=FALSE, eval=TRUE, warning=FALSE, error=FALSE, message=FALSE, fig.width=9, fig.height=9, fig.cap="Function call graph"
#library(RColorBrewer)
#library(mvbutils)
#foodweb(border = TRUE, boxcolor = "skyblue", textcolor = "black", cex = 0.8, lwd=2)

#'Run following command to generate documentation
#'
#' rmarkdown::render("~/projects/darwin/src/echo_tracking.R", output_format = "pdf_document")
