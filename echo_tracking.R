#!---------------------------------------------------------------------------
#' @author Bhupendra Raut (www.baraut.info)
#' @brief This R script reads netCDF files containing echo heights and rain classification
#' identify objects and estimates flow vectors, using phase correlation with FFT (Leese et al 1971).
#' Then for each object the search region is predicted for next frame and the objects in that
#' region (in frame2) are assign to the object in frame 1.
#' @reference :Leese, John A., Charles S. Novak, and Bruce B. Clark.
#'           "An automated technique for obtaining cloud motion from geosynchronous
#'           satellite data using cross correlation."
#'           Journal of applied meteorology 10.1 (1971): 118-132.

#' @todo
#' 1. Check overlapping pixels to identify better match.
#'
#==========================================================================================
# Start the clock!
start_time <- proc.time()

library(ncdf4)
library(EBImage)  # for bwlabel,
library(plot3D)
library(spatstat) #for smoothing
library(stringr)
library(plyr)
library(clue)  #solve_LSAP()

#----------------------------------------------------------------------fucntions
## Plots image with objects labels
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
        text(x=xvalues[obj_centerIndex[1]], y=yvalues[obj_centerIndex[2]], toString(object), cex=0.7)
    }
}

## Returns labeled image after removing single pixel objects
clear_onePix_objects <- function(label_image) {
    size_table <- table(label_image[label_image>0]) # remove zero values
    onePix_objects <- as.vector(which(size_table==1))

    for(obj in onePix_objects){
        label_image <- replace(label_image, label_image == obj, 0.0)
    }

    label_image <- bwlabel(label_image)
    invisible(label_image)
}


## Takes in a labeled image and finds the redius and the center of the given object.
get_objExtent <- function(labeled_image, obj_label) {
    #center indices of the object assuming it is a rectangle
    obj_index <- which(labeled_image==obj_label, arr.ind = TRUE)

    rlength <- max(obj_index[, 1]) - min(obj_index[, 1]) + 1
    clength <- max(obj_index[, 2]) - min(obj_index[, 2]) + 1

    obj_radius<- max(c(rlength, clength))/2 #maximum possible object radius
    obj_center <- c(min(obj_index[, 1])+obj_radius, min(obj_index[, 2]) + obj_radius)
    return(list(obj_center=obj_center, obj_radius=obj_radius))
}


## Takes in object info (radius and center) and two images to estimate ambient flow.
# margin is the additional region arround the object used to comput the flow vectors.
get_objAmbientFlow <- function(obj_extent, img1, img2, margin) {
    #coordinates of the flowregion
    r1 <- obj_extent$obj_center[1] - obj_extent$obj_radius - margin
    r2 <- obj_extent$obj_center[1] + obj_extent$obj_radius + margin
    c1 <- obj_extent$obj_center[2] - obj_extent$obj_radius - margin
    c2 <- obj_extent$obj_center[2] + obj_extent$obj_radius + margin

    dims <- dim(img1)
    if(r1<=0 || c1 <=0 || r2>dims[1] || c2 > dims[2]){
        return(c(0, 0))
    }

    flow_region1 <- img1[r1:r2, c1:c2]
    flow_region2 <- img2[r1:r2, c1:c2]

    return(fft_flowVectors(flow_region1, flow_region2))
}


#computs cross-covariance using FFT
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
}

## Rearranges the crossCov matrix so that 'zero' frequency or DC component
#  is in the middle of the matrix.
#   This function is adopted from following discussion on stackOverflow
#   http://stackoverflow.com/questions/30630632/performing-a-phase-correlation-with-fft-in-r
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

## Estimates flow vectors in two images
fft_flowVectors <- function (im1, im2) {
    if(max(im1)==0 || max(im2)==0){
        return(c(0, 0))
    }

    im1 <- replace(im1, im1==0, runif(1)) # add noise to image background
    im2 <- replace(im2, im1==0, runif(1)) # This may help when objects are bigger

    crossCov <- fft_crossCov(im1, im2)
    cov_shifted <- fft_shift(crossCov)
    cov_smooth <- blur(as.im(cov_shifted))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}

get_std_flowVector<-function(obj_extent, img1, img2, margin, magnitude){
    shift <- get_objAmbientFlow(obj_extent, img1, img2, margin)
    shift <- replace(shift, shift > magnitude, magnitude)
    shift <- replace(shift, shift < magnitude*-1, magnitude*-1)
    return(shift)
}


##
predict_searchExtent <- function(obj1_extent, shift){
    shifted_center <- obj1_extent$obj_center + shift
    search_radius <-5

    x1 <- shifted_center[1] -search_radius
    x2 <- shifted_center[1] +search_radius
    y1 <- shifted_center[2] -search_radius
    y2 <- shifted_center[2] +search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}


## Retuns  Euclidean distance between two vectors or matrices
euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}

## search box shoudn't be outside the image or very small.
check_searchBox <- function(search_box, sample_img){
    dims <- dim(sample_img)
    if(search_box$x1 <= 0){
        search_box$x1 <- 1
    }
    if(search_box$y1 <= 0){
        search_box$y1 <- 1
    }
    if(search_box$x2 > dims[1]){
        search_box$x2 <- dims[1]
    }
    if(search_box$y2 > dims[2]){
        search_box$y2 <- dims[2]
    }

    #search box should be large enough
    if(search_box$x2-search_box$x1 < 5 || search_box$y2-search_box$y1 < 5 ){
        return(NA)
    } else {
        return(search_box)
    }
}

## Given the search box and image2, finds objects in the region
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


## Retuns ratio (>=1) of bigger number to smaller number when given two number.
get_ratio<-function(x, y){
    if(x>=y) return(x/y)
    else if (y>x) return(y/x)
}


#computes discrepancy for objects
get_discrepancy <- function(obj_found, temp2, search_box, obj1_extent) {
    dist_pred <- c(NULL)
    dist_actual <- c(NULL)
    for(target_obj in obj_found){
        target_extent <- get_objExtent(temp2, target_obj)
        euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred)
        dist_pred <- append(dist_pred, euc_dist)

        euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
        dist_actual <- append(dist_actual, euc_dist)
        size_changed <- get_ratio(target_extent$obj_radius, obj1_extent$obj_radius) #change in size

        discrepancy <- dist_pred + size_changed + dist_actual

    }
    return(discrepancy)
}



#----------------------------------------------------------------Calling Program
setwd("~/data/darwin_radar/2d/")
infile_name <- "./cpol_2D_0506.nc"

ncfile <- nc_open(infile_name)
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
dbz_height <- ncvar_get(ncfile, varid = "zero_dbz_cont", start = c(1, 1, 1), count = c(-1, -1, 100))
steiner <- ncvar_get(ncfile, varid = "steiner_class", start = c(1, 1, 1), count = c(-1, -1, 100))

dbz_height <-replace(dbz_height, steiner != 2, 0.0)      #set non-convective pixels to zeros
dbz_height <- replace(dbz_height, is.na(dbz_height), 0.0)     #remove NAs
labeled_echo <- bwlabel(dbz_height)               #identify and label objects


rm(dbz_height)
rm(steiner)



scan <-96
search_margin <- 5 #pixels
flow_margin <- 10 #pixels
stdFlow_mag <- 3
large_num <- 100000

temp1<-labeled_echo[, , scan]
temp2<-labeled_echo[, , scan+1]

temp1 <- clear_onePix_objects(temp1)
temp2 <- clear_onePix_objects(temp2)


nObjects1 <- max(temp1) #objects in first image
nObjects2 <- max(temp2) #objects in first image
if(nObjects2==0 || nObjects1==0){
    stop("No echoes to track!!!")
}

if(nObjects2 > nObjects1){
    obj_match <- matrix(large_num, nrow = nObjects1, ncol = nObjects2)
} else {
    obj_match <- matrix(large_num, nrow = nObjects1, ncol = nObjects1)
}


for(obj_id1 in 1:nObjects1) {

    obj1_extent <- get_objExtent(temp1, obj_id1)
    shift <- get_std_flowVector(obj1_extent, temp1, temp2, flow_margin, stdFlow_mag)
    print(paste("fft shift", toString(shift)))

    search_box <- predict_searchExtent(obj1_extent, shift)
    search_box <- check_searchBox(search_box, temp2)
    obj_found <- find_objects(search_box, temp2)

    if(is.na(obj_found[1]) || max(obj_found)==0) {
        obj_id2 <- 0
        dist_pred <- NA
        dist_actual <- NA
        discrepancy <- NA
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0

        if(length(obj_found)==1){ # if this is the only object
            discrepancy <- get_discrepancy(obj_found, temp2, search_box, obj1_extent)
            if(discrepancy < 10) discrepancy <- 0 #lower the discrepancy if not too large

        } else { # when more than one objects
            discrepancy <- get_discrepancy(obj_found, temp2, search_box, obj1_extent)
        }
    }


    if(discrepancy >15 || is.na(discrepancy)){
        obj_match[obj_id1, obj_found] <- large_num
    } else {
        obj_match[obj_id1, obj_found] <- discrepancy
    }
    print(paste(obj_id1, "==>", toString(obj_found)))
}

pairs <- solve_LSAP(obj_match)

## remove bad matching
for(pair in 1:length(pairs)){
    if(obj_match[pair, pairs[pair]] >15){
        pairs[pair] <- 0
    }
}




#plot
pdf(paste("object_label_", scan, ".pdf", sep=""), width=12, height=8)
par(mfrow=c(1,2))
plot_objects_label(temp1, x, y)
plot_objects_label(temp2, x, y)
dev.off()

# Stop the clock
print(proc.time() - start_time)
