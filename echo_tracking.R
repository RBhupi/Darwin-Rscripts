#' ---
#' title: "Tracking Convection Echoes in CPOL Radar Data"
#' author: "Bhupendra Raut"
#' date: "December 01, 2016"
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
#' John A Leese, Charles S Novak, and Bruce B Clark. An automated technique for obtaining
#' cloud motion from geosynchronous satellite data using cross correlation.
#' Journal of applied meteorology, 10(1):118–132, 1971.
#'
#' 1. The unique identity numbers are produced for new objects and was kept
#' constant through out the life of that object.
#' 2. All the distances used for tracking are in pixels, hence this code can be
#' customised for newer problems of similar kind.
#' 3. Splitting/Merging is considered in this version.
#' Uid of the origin/product echos written in arrays named `origin' and `merged'.
#'
#'
#'Issues: 1. NetCDF files missing value issue persists. http://stackoverflow.com/q/37202454/1227454
#'      2. Check correctness of the output "survival_stat" for all the conditions.
#'
#' ToDo
#' 1. Check overlapping pixels to identify better match.
#' 2. Use vertical profile of clouds for better assignment (KF)
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
#----------------------------------------------------------------------functions

#' Plots image with objects labels. 
#' 
#' This is used in development stage to test images when processed 1-by-1.
plot_objects_label <- function(labeled_image, xvalues, yvalues){
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

#' Returns a single radar scan from the input netcdf file.
#' 
#' Smaller objects are removed and the rest are lebeled.
#' @param ncfile input netcdf file object from ncdf4 package.
#' @param var_id string name of the varibale in the file.
#' @param scan_num index of frame to be read from the file.
#' @param min_size in pixels objects smaller than this will be removed. Default 2.
#' @export
get_filteredFrame <- function(ncfile, var_id, scan_num, min_size=2) {
    frame <- ncvar_get(nc=ncfile, varid = var_id, 
                        start = c(1, 1, scan_num), count = c(-1, -1, 1))
    labeled_echo <- bwlabel(frame)          #label objects
    frame <-clear_smallEchoes(labeled_echo, min_size)
    invisible(frame)
}


#' Removed objects smaller than the given size.
#' 
#' Takes in labeled image removes objects smaller than min_size and returns re-labeled image.
#' Pixels attached diagonally are not considered as continuouse part of the object. 
clear_smallEchoes <- function(label_image, min_size) {
    size_table <- table(label_image[label_image>0]) # remove zero values
    onePix_objects <- as.numeric(names(which(size_table < min_size)))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, 0.0)
    }

    label_image <- bwlabel(label_image)
    invisible(label_image)
}



#'returns pairs of object ids from both frames.
#'
#' Given two images, the function identifies the matching 
#' objects and pair them appropriatly using Hungarian method and desparity. 
#' 
#' @export
#' @seealso  code{disparity function}.
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
#' 
#' The bad matching is when disparity is more than the set value.
match_pairs <- function(obj_match) {
    obj_match[obj_match<0] <- 0
    pairs <- solve_LSAP(obj_match)
    pairs <- as.vector(pairs)
    # remove bad matching
    for(pair in 1:length(pairs)){
        if(obj_match[pair, pairs[pair]] > max_desparity){
            pairs[pair] <- 0
        }
    }
    return(pairs)
}


#' Matches all the obejects in image1 to the objects in image2.
#' 
#' This is the main function to be called on two sets of radar images, for tracking.
locate_allObjects <- function(image1, image2) {
    nObjects1 <- max(image1) #objects in first image
    nObjects2 <- max(image2) #objects in second image


    if(nObjects2==0 || nObjects1==0){
        stop("No echoes to track!!!")
    }

    obj_match <- matrix(large_num, nrow = nObjects1, ncol = max(nObjects1, nObjects2))

    # here we match each object in image1 to all the near-by objects in image2.
    for(obj_id1 in 1:nObjects1) {
        obj1_extent <- get_objExtent(image1, obj_id1) #location, ind and radius
        shift <- get_std_flowVector(obj1_extent, image1, image2, flow_margin, maxFlow_mag)

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


#' FFT shifts are corrected using last headings.
#' 
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


#' Returns object extent properties.
#' 
#' Takes in a labeled image and finds the radius, area and the center of the given object.
get_objExtent <- function(labeled_image, obj_label) {
    #center indices of the object assuming it is a rectangle
    obj_index <- which(labeled_image==obj_label, arr.ind = TRUE)

    xlength <- max(obj_index[, 1]) - min(obj_index[, 1]) + 1
    ylength <- max(obj_index[, 2]) - min(obj_index[, 2]) + 1

    obj_radius<- max(c(xlength, ylength))/2 #maximum possible object radius
    #obj_center <- c(min(obj_index[, 1]) + obj_radius, min(obj_index[, 2]) + obj_radius)

    #definition of object center based on median, This needs some thought
    obj_center <- c(median(obj_index[, 1]), median(obj_index[, 2]))

    obj_area <- length(obj_index[, 1])  #size in pixels

    obj_extent<-list(obj_center=obj_center, obj_radius=obj_radius,
                     obj_area=obj_area, obj_index=obj_index)

    return(obj_extent)
}




#' Computes flow in the vicinity of the object.
#' 
#' Takes in object info (radius and center) and two images to estimate ambient flow using FFT phase correlation.
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
#' 
#' Flow vectors magnitude is clipped to given magnitude to controll erratic output from FFT flow.
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


#' Rearranges the crossCov matrix 
#' 
#' Rearranges the crossCov matrix so that 'zero' frequency or DC component
#'  is in the middle of the matrix. Taken from stackoverflow Que. 30630632
fft_shift <- function(fft_mat) {
    if(class(fft_mat)=='matrix') {
        rd2 <- floor(nrow(fft_mat)/2)
        cd2 <- floor(ncol(fft_mat)/2)

        # Identify the first, second, third, and fourth quadrants
        q1 <- fft_mat[1:rd2,1:cd2]
        q2 <- fft_mat[1:rd2,(cd2+1):ncol(fft_mat)]
        q3 <- fft_mat[(rd2+1):nrow(fft_mat),(cd2+1):ncol(fft_mat)]
        q4 <- fft_mat[(rd2+1):nrow(fft_mat),1:cd2]

        # rearrange the quadrants
        centered.t <- rbind(q4,q1)
        centered.b <- rbind(q3,q2)
        centered <- cbind(centered.b,centered.t)

        invisible(Re(centered))
    } else {
        stop("input to fft_shift() should be a matrix")
    }
}


#' predict search region.
#' 
#' Predicts search extent/region for the object in image2 given the image shift.
predict_searchExtent <- function(obj1_extent, shift, search_radius){
    shifted_center <- obj1_extent$obj_center + shift

    x1 <- shifted_center[1] -search_radius
    x2 <- shifted_center[1] +search_radius
    y1 <- shifted_center[2] -search_radius
    y2 <- shifted_center[2] +search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}


#' checks that the search box is in the domain.
#' 
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


#' Returns vector of objects ids in the given reion.
#' 
#' Given the search box and image2, returns objects in the region.
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


#' Returns Disparity of all the objects in the region.
#' 
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
            if(disparity <= 3) disparity <- 0 #lower the disparity if not too large

        } else { # when more than one objects to match
            disparity <- get_disparity(obj_found, image2, search_box, obj1_extent)
        }
    }
    return(disparity)
}


#' Corrects and saves disparity for the object matching stage.
#' 
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


#' Actually computes desparity for a single object.
#' 
#' Check how it is computed for detail.
#' This parameter has most effect on the acccuracy of the tracks.
get_disparity <- function(obj_found, image2, search_box, obj1_extent) {
    dist_pred <- c(NULL)
    dist_initial <- c(NULL)
    change <- c(NULL)
    overlap <- c(NULL)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(image2, target_obj)
        overlap_area <- check_bigOverlap(obj1_extent, target_extent)
        overlap <- append(overlap, overlap_area)

        euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist)

        euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_initial <- append(dist_initial, euc_dist)
        size_changed <- get_sizeChange(target_extent$obj_area, obj1_extent$obj_area) #change in size
        change <- append(change, size_changed)

    }
    
    #This is crucial parameter that affect the results
    disparity <- dist_pred + change + dist_initial - overlap


    return(disparity)
}

#' checks overlapping regoin for big objects in both the frames.
#' 
#' Checks overlapping area in pixels, size of the object and return if overlapping is considerable.
check_bigOverlap <- function(obj_extend, target_extend){
    duplicates <- duplicated(rbind(obj_extend$obj_index, target_extend$obj_index))
    overlap_area <- length(duplicates[duplicates==TRUE])
    if(obj_extend$obj_area > big_obj_size & overlap_area >= obj_extend$obj_area/2){
        return(overlap_area)
    } else 
        return(0)
}

#' standard Euclidean distance.
#' 
#' Returns  Euclidean distance between two vectors or matrices.
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}


#' returns size change between the frames.
#' 
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




#' Creates output netcdf file for radar echo trajectories. 
#'
#' This is the longest function.
#' @param ofile name and path of the output file.
#' @param max_obs longest possible track record. Deafult maximum 65 observation per track.
#' @export
create_outNC <- function(ofile, max_obs) {
    if(file.exists(ofile)){
        print(paste("removing existing file", basename(ofile)))
        file.remove(ofile)
    }
    deflat <- 9

    dim_echo <- ncdim_def("echo_id", vals=1, units = "", unlim = TRUE,
                          longname = "unique id of convection echo", create_dimvar = TRUE)

    dim_obs <- ncdim_def("records", vals = seq(max_obs), units="",
                         longname = "observation records")

    dim_time <- ncdim_def("time", vals=1, units = "seconds since 1970-01-01 00:00:00 UTC",
                          longname = "time of the scan", unlim = TRUE, create_dimvar = TRUE)

    dim_stat <- ncdim_def("stat", vals = seq(4), units="", longname = "object survival vector; lived, died, born, total")

    # Define Variables
    var_survival <- ncvar_def("survival", units = "", longname = "survival from the last scan",
                              dim=list(dim_stat, dim_time), missval = -999, prec="integer",
                              compression = deflat, shuffle = TRUE)

    var_dur <- ncvar_def("duration", units = "", longname = "duration of echo in time-steps",
                         dim=dim_echo, missval = -999, prec="integer",
                         compression = deflat, shuffle = TRUE)

    var_origin <- ncvar_def("origin", units="", longname = "id from which the echo split up.",
                            dim=dim_echo, missval = 0, prec = "integer")

    var_merged <- ncvar_def("merged", units="", longname = "id in which the echo merged.",
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

    #var_ncg <- ncvar_def("Cu_cong", units = "pixels", longname = "num of Cu Congestus pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    #var_ncb <- ncvar_def("Cu_deep", units = "pixels", longname = "num of deep convection pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    #var_nco <- ncvar_def("Cu_over", units = "pixels", longname = "num of overshooting convection pixels",
    #                     dim = list(dim_obs, dim_echo), missval = -999, prec = "integer",
    #                     compression = deflat, shuffle = TRUE)

    var_list <- list(var_time, var_survival, var_dur, var_origin, var_merged, var_xdist, var_ydist,
                     var_x, var_y, var_npix) #, var_ncg, var_ncb, var_nco)


    outNC <- nc_create(filename = ofile, vars = var_list)

    write_settingParms_toNC(outNC)

    #for CF standards
    ncatt_put(outNC, varid = "echo_id", attname = "cf_role", attval = "trajectory_id")
    ncatt_put(outNC, varid = 0, attname = "featureType", attval = "trajectory")

    description <- paste("The CPOL (Darwin) radar echoes of convective types were separated using Steiner classification scheme and tracked.",
                         "Merging and splitting is added with echo ids.")

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


#' Writes all the setting parameters (as attributes) for the tracking. 
#' 
#' These parameters affect the sensitivity of the tracks, mergers and split definitions etc.
write_settingParms_toNC <- function(outNC){
    ncatt_put(outNC, varid = 0, attname = "search_margin", attval =search_margin, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "flow_margin", attval =flow_margin, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "maxFlow_magnitude", attval =maxFlow_mag, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "min_echoSize_toTrack", attval =min_size, prec = "short")
    ncatt_put(outNC, varid = 0, attname = "max_desparity", attval =max_desparity, prec = "short")
}



#' Writes properties and uids for all objects into output netcdf file.
#' 
#' @param outNC output netcdf file object from function \code{create_outNC()}
#' @param current_objects output of \code{update_current_objects}
#' @param obj_props output of \code{get_object_prop()}
#' @param obs_time time of first scan in POSIX format. units="seconds since 1970-01-01".
#' @export
write_update<-function(outNC, current_objects, obj_props, obs_time){
    nobj <- length(current_objects$id1) #num of objects in frame1

    for(object in seq(nobj)){
        nc_start <- c(current_objects$obs_num[object], current_objects$uid[object])
        nc_count <- c(1, 1)

        ncvar_put(outNC, varid = "record_time", obs_time, start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x", obj_props$x[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y", obj_props$y[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "x_dist", obj_props$xdist[object], start = nc_start, count = nc_count)
        ncvar_put(outNC, varid = "y_dist", obj_props$ydist[object], start = nc_start, count = nc_count)

        ncvar_put(outNC, varid = "area", obj_props$area[object],  start = nc_start, count = nc_count)
    }

    write_duration(outNC, current_objects)
}


#' Write duration of dead objects in output NC file.
#' 
#' Writes number of observations for dead objects. Duration is in time-steps.
write_duration <- function(outNC, current_objects){
    nobj <- length(current_objects$id1)

    for (obj in seq(nobj)){
        if(current_objects$id2[obj]==0){
            ncvar_put(outNC, varid = "duration", current_objects$obs_num[obj],
                      start=current_objects$uid[obj], count=1)
            ncvar_put(outNC, varid = "origin", current_objects$origin[obj],
                      start=current_objects$uid[obj], count=1)

            #check for merging
            merged_in <- check_merging(obj, current_objects, obj_props)
            ncvar_put(outNC, varid="merged", merged_in,
                      start=current_objects$uid[obj], count=1)
        }
    }
}


#' Checks possible merging of the dead objects.
#'
#'This function takes in two R-lists containing information about current objects
#' in the frame1 and their properties, such as center location and area. If the
#' I am using an arbitrary crieterion for merging. If euclidean distance between centers
#' of the two objects c_dist < or = to r=sqrt(area), then merging is considered.
#' Here, if we assume square objects, then the r is length of a sides of the square.
check_merging<-function(dead_obj_id1, current_objects, object_props){
    nobj_frame1 <- length(current_objects$id1)
    c_dist_all <- NULL
    checked_id1 <- NULL

    #' If all objects are dead in frame2 then no merging happened.
    if(all(current_objects$id2==0)) return(0)

    for(check_obj in seq(nobj_frame1)){
        # skip checking dead objects of frame2
        if(current_objects$id2[check_obj]!=0){
            dead_xy <- c(obj_props$x[dead_obj_id1], obj_props$y[dead_obj_id1])
            merge_xy <- c(obj_props$x[check_obj], obj_props$y[check_obj])
            c_dist <- euclidean_dist(merge_xy, dead_xy)
            if(c_dist < sqrt(object_props$area[check_obj])){
                c_dist_all <-append(c_dist_all, c_dist)
                checked_id1 <- append(checked_id1, check_obj)
            }

        }
    }
    #if this is null, then no merging, else nearest object is the product of merging.
    if(is.null(c_dist_all)) return(0)
    else {
        product_id1 <- checked_id1[which(c_dist_all==min(c_dist_all))]
        return(current_objects$uid[product_id1])
    }


}


#' Write survival stats 
#' 
#' write number of lived, dead and born objects to the file for each scan.
#' 
#' @export
write_survival <- function(outNC, survival_stat, time, scan){
    if(!is.atomic(survival_stat)){
        survival_stat <- unlist(survival_stat, use.names = FALSE)
    }

    ncvar_put(outNC, varid = "survival", vals = survival_stat, start = c(1, scan), count = c(4, 1))
    ncvar_put(outNC, varid = "time", vals = time, start = scan, count=1)
}


#' Returns a dataframe for objects with unique ids and their corresponding ids in frame1 and frame2.
#' 
#' This function is called when new rainy scan is seen after the period of no rain or the first time.
#' @param first_frame First image for tracking.
#' @param second_frame Second image for tracking.
#' @param pairs output of \code{get_match_pairs()}
#' @export
init_uids <- function(first_frame, second_frame, pairs){
    nobj <- max(first_frame) #number of objects in frame1
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj)               #id1
    objects_mat[, 2] <- next_uid(count = nobj) #unique ids
    objects_mat[, 3] <- as.vector(pairs) #as they are in frame2
    objects_mat[, 4] <-rep(1, nobj)     #observation number for the echo
    objects_mat[, 5] <-rep(0, nobj)

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "origin")
    current_objects <- attach_xyheads(first_frame, second_frame, current_objects)
    return(current_objects)
}


#' saves last x y movements of the objects.
#' 
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


#' return center indices of the object.
#' 
#' Returns indices of center pixel of the given object id from a labeled image. 
#' This may be done in better way for non-oval objects.
get_objectCenter<-function(obj_id, labeled_image){
    obj_index <- which(labeled_image==obj_id, arr.ind = TRUE)
    center_x <- floor(median(obj_index[, 1])) #center column
    center_y <- floor(median(obj_index[, 2])) #center row
    return(c(center_x, center_y))
}


#' Removes dead objects, updates living objects and assign new uids to new born objects.
#' 
#' Also, updates number of valid observations for each echo.
#' This function is called when rain continues from the last frame.
#' This is a complicated function to understand.
#'
#' @details See how the pairs vector looks like for a real case. The pairs
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
#' @export
update_current_objects <- function(frame1, frame2, pairs, old_objects){
    nobj <- max(frame1)
    objects_mat <- matrix(data = NA, ncol = 5, nrow = nobj)

    objects_mat[, 1] <- seq(nobj) # this is id1 at current step

    for (obj in seq(nobj)){
        if(obj %in% old_objects$id2){ # but same was id2 in the last step
            # so they should get same uid as last time
            objects_mat[obj, 2] <- old_objects$uid[old_objects$id2==obj]
            objects_mat[obj, 4] <- old_objects$obs_num[old_objects$id2==obj] + 1
            objects_mat[obj, 5] <- old_objects$origin[old_objects$id2==obj]
        } else {
            objects_mat[obj, 2] <- next_uid()
            objects_mat[obj, 4] <- 1 #first observation of the echo
            objects_mat[obj, 5] <- get_origin_uid(obj, frame1, old_objects)
        }
    }

    objects_mat[, 3] <- as.vector(pairs) #match as they are in frame2

    current_objects <- data.frame(objects_mat, row.names = NULL)
    colnames(current_objects) <- c("id1", "uid", "id2", "obs_num", "origin")
    current_objects <- attach_xyheads(frame1, frame2, current_objects)
    invisible(current_objects)
}


#' Find id of the parent of the new born object.
#' 
#' returns unique id of the origin (or zero) for given object in frame1.
#' Also remember that old object id2 is actual id1 in frame1, as we still have
#' to update the object_ids.
get_origin_uid<-function(obj, frame1, old_objects){
    origin_id <- find_origin(obj, frame1)
    if (origin_id==0) return(0)

    origin_index <- which(old_objects$id2==origin_id)

    # If it is first observation of the object then it will not be recorded in
    # old_objects$id2, ans it will not be suitable as the origin.
    if(!(origin_id %in% old_objects$id2)) return(0)

    origin_uid <- old_objects$uid[origin_index]
    return(origin_uid)
}


#' Checks for parent in the vicinity.
#' 
#' This function checks near by objects in the frame for the given new-born object.
#' origin is an object which existed before the new born objects,
#' has comparable or larger size and is close enough to the offspring.
find_origin <- function(id1_newObj, frame1){
    if(max(frame1)==1) return(0) # If there is only one object, then dont look for origin

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

    #if both are same call it the origin
    if(big_ratio_obj==big_diff_obj)
        return(big_diff_obj[1])
    else
        return(big_diff_obj[1])
    # NOTE: 1. At this time we are calling big_diff_obj as origin in all the situations.
    # This looks like a good first guess. But if needed we can make it more
    # complex and use ratio and size_diff as cost function.
    # 2. We are not considering the possibility of multiple potential origins
    # beyond this point.
}


#' 
#' 
#' Returns sequence of next unique ids and increament the uid_counter.
next_uid<-function(count=1){
    this_uid <- uid_counter + 1:count
    uid_counter <<- uid_counter + count
    return(this_uid)
}


#' Return all the object's size, location and classification info.
#' 
#' xyDist should be a list of x_dist and y_dist in km.
#' @export
get_objectProp <- function(image1, xyDist){
    objprop <- c(NULL)
    nobj <- max(image1)

    for(obj in seq(nobj)){
        obj_index <- which(image1==obj, arr.ind = TRUE)
        objprop$id1 <- append (objprop$id1, obj)  #id in frame1
        objprop$x <- append(objprop$x, floor(median(obj_index[, 1]))) #center column
        objprop$y <- append(objprop$y, floor(median(obj_index[, 2]))) #center row
        objprop$area <- append(objprop$area, length(obj_index[, 1]))
    }
    objprop <- attach_xyDist(objprop, xyDist$x, xyDist$y)
    invisible(objprop)
}




#' 
#' 
#' Attaches y and x distance from radar in km to object location indices
attach_xyDist<-function(obj_props, xdist, ydist){
    obj_props$xdist <- xdist[obj_props$x]
    obj_props$ydist <- ydist[obj_props$y]
    invisible(obj_props)
}


#' Returns a list with number of objects lived, died and born between the current and the previousstep.
#' 
#' @export
survival_stats <- function(pairs, num_obj2) {
    pairs_vec <- as.vector(pairs)
    obj_lived <- length(pairs_vec[pairs_vec>0])
    obj_died <- length(pairs_vec)-obj_lived
    obj_born <- num_obj2 - obj_lived
    return(list(lived=obj_lived, died=obj_died, born=obj_born, total=num_obj2))
}
#+ echo=FALSE
#==============================================================================#


#+ echo=TRUE
#'----------------------- Settings for tracking method ------------------------#
search_margin <- 4          #pixels
flow_margin <- 4            #pixels
maxFlow_mag <- 5            #fft_flow will not be faster than this
min_signif_movement <- 2    #not used at this time
large_num <- 1000           #a large number for Hungarian method
max_obs<- 100                #longest recoreded track (eles show error).
min_size <- 25               #objects smaller than this will be filter
max_desparity <- 20         # two objects with more desparity than this value, are not same.
big_obj_size <- 25
#==============================================================================#

var_name<-"DBZc"


#'Here we call the above functions to track convective echoes in the subsequent images.
#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE

setwd("/home/bhupendra/data/netcdf_solapur/2dNC")
file_list <- Sys.glob(paths = "./SolaRadar2d_lev3.nc")
#file_list <- Sys.glob(paths = "./cpol_2D_2004-11-03.nc")

for(infile_name in file_list){
    outfile_name <- str_replace(infile_name, ".nc", "_tracks_testDelete.nc")
    #outfile_name <- "~/Desktop/test_tracks.nc"
    print(paste("Opening output file", outfile_name))

    outNC <- create_outNC(outfile_name, max_obs)
    uid_counter <- 0            #(a global variable) to start uid with 1, count from 0.

    #read x, y and time from the file
    ncfile <- nc_open(infile_name)
    x <- ncvar_get(ncfile, varid = "x")
    y <- ncvar_get(ncfile, varid = "y")


    time <- ncvar_get(ncfile, varid="time")
    #time <- change_baseEpoch(time, From_epoch = as.Date("2004-01-01"))


    start_scan <- 79
    print(paste("start scan = ", start_scan))
    end_scan <- 100 #length(time)

    newRain <- TRUE         #is this new rainy scan after dry period?

    print(paste("Total scans in this file", end_scan-start_scan+1))
    pb = txtProgressBar(min =start_scan, max = end_scan, initial = 1, style = 3) #progress bar

    #' To continuousely track the echoes in the images,
    #' we read the first scan and save it in to array named frame2 then in the loop,
    #' this frame will be copied to the array named frame1 and next scan will be frame2.
    #'
    #' frame1 <-- frame2 <-- next scan  (repeat)
    #+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
    #frame2 <- get_filteredFrame(ncfile, start_scan, min_size)
    #class2 <- get_classFrame(ncfile, start_scan) #classifictaion
    frame2 <- get_filteredFrame(ncfile, var_name, start_scan, min_size)
    
    for(scan_ind in (start_scan+1):end_scan){
        setTxtProgressBar(pb, scan_ind) #progress bar

        # save earlier frame2 to frame1
        frame1 <- frame2
        #class1 <- class2

        # and read next scan to frame2
        #frame2 <- get_filteredFrame(ncfile, scan, min_size)
        #class2 <- get_classFrame(ncfile, scan)
        frame2 <- get_filteredFrame(ncfile, var_name, scan_ind, min_size)

        #if this is the last scan make it zero. This kills all the objects.
        if(scan_ind==end_scan){
            frame2 <- replace(frame2, frame2>0, 0)
        }


        #skip if no echoes in frame 1
        if(max(frame1, na.rm = TRUE)==0){       #if no echoes in frame1
            newRain = TRUE                      #next rain will be newRain
            if(exists("current_objects")){
                rm(current_objects)
            }

            #write zeros for empty frames
            write_survival(outNC, survival_stat = rep(0, 4),
                          time = time[scan_ind], scan = scan_ind)

            next
        }

        
        # track when echoes are present in frame1
        pairs <- get_matchPairs(frame1, frame2)
        obj_props <- get_objectProp(frame1, list(x=x, y=y)) #of frame1

        #when echoes are found in frame1 and it is newRain, init uids
        if(newRain){
            current_objects <- init_uids(frame1, frame2, pairs) #init ids and return list of objects

            #for newRain, all the objects are born in this frame1
            num_obj1 <- max(frame1)
            survival <- c(rep(0, 2), num_obj1, num_obj1)

            #except for the first scan where values should be missing
            if(scan_ind==start_scan+1)
                survival <- c(rep(-999, 3), num_obj1)

            write_survival(outNC, survival_stat = survival,
                           time = time[scan_ind-1], scan = scan_ind-1)

            newRain <- FALSE

        } else {            #else update old ids
            current_objects <- update_current_objects(frame1, frame2, pairs, current_objects)
        }

        write_update(outNC, current_objects, obj_props, time[scan_ind-1]) #for frame1

        #Survival for frame2
        num_obj2 <- max(frame2)
        obj_survival <- survival_stats(pairs, num_obj2)
        write_survival(outNC, survival_stat = obj_survival, time = time[scan_ind], scan = scan_ind)
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
#library(RCu_overlorBrewer)
#library(mvbutils)
#foodweb(border = TRUE, boxcolor = "skyblue", textcolor = "black", cex = 0.8, lwd=2)

#'Run following command to generate documentation
#'
#' rmarkdown::render("~/projects/darwin/src/echo_tracking.R", output_format = "pdf_document")

