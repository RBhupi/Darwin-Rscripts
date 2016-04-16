#@ToDo
# 1. Check overlapping pixels to identify better match.
# 2. check the actual area of the echo for both images.
#
#



library(ncdf4)
library(EBImage)
library(plot3D)
library(spatstat) #for smoothing
library(stringr)
library(plyr)

plot_objects_label<-function(labeled_image, xvalues, yvalues){
    image2D(replace(labeled_image, labeled_image==0, NA), x=xvalues, y=yvalues)
    grid()

    for(object in 1:max(labeled_image)) {
        #center indices of the object assuming it is a rectangle
        obj_id1 <- object

        obj_index <- which(labeled_image==obj_id1, arr.ind = TRUE)

        r1 <- min(obj_index[, 1])
        r2 <- max(obj_index[, 1])
        c1 <- min(obj_index[, 2])
        c2 <- max(obj_index[, 2])


        obj_centerIndex <- c((r1+r2)/2, (c1+c2)/2)

        text(x=xvalues[obj_centerIndex[1]], y=yvalues[obj_centerIndex[2]], toString(object), cex=0.7)
    }
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

## Rearranges the crossCov matrix so that 'zero' frequency or DC component is in the middle of the matrix.
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
    crossCov <- fft_crossCov(im1, im2)
    cov_shifted <- fft_shift(crossCov)
    cov_smooth <- blur(as.im(cov_shifted))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}



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



scan <- 97
search_margin <- 4 #pixels
temp1<-labeled_echo[, , scan]
temp2<-labeled_echo[, , scan+1]




for(obj_id1 in 1:max(temp1)) {
    obj1_extent <- get_objExtent(temp1, obj_id1)
    shift <- get_objAmbientFlow(obj1_extent, temp1, temp2, search_margin)

    search_box <- predict_searchExtent(obj1_extent, shift)
    search_area <- temp2[search_box$x1:search_box$x2, search_box$y1:search_box$y2]

    obj_found <- unique(as.vector(search_area))
    #print(paste(obj_id1, "==>", toString(obj_found)))


    if(max(obj_found)==0) {
        obj_id2 <- 0
    } else {
        obj_found <- obj_found[obj_found>0] #remove 0

        dist_pred <- c(NULL)
        dist_actual <- c(NULL)

        for(target_obj in obj_found){
            target_extent <- get_objExtent(temp2, target_obj)
            euc_dist<- euclidean_dist(target_extent$obj_center, search_box$center_pred) #find the distsnce from predicted
            dist_pred <- append(dist_pred, euc_dist)

            euc_dist<- euclidean_dist(target_extent$obj_center, obj1_extent$obj_center)
            dist_actual <- append(dist_actual, euc_dist)
        }

        obj_id2 <- obj_found[which(dist<3) && which(dist_actual<4)]
    }

    print(paste(obj_id1, "==>", toString(obj_id2)))
}


#plot
pdf("object_label.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot_objects_label(temp1, x, y)
plot_objects_label(temp2, x, y)
dev.off()
##
predict_searchExtent <- function(obj1_extent, shift){
    shifted_center <- obj1_extent$obj_center + shift

    #x1 <- shifted_center[1] - obj1_extent$obj_radius -2
    #x2 <- shifted_center[1] + obj1_extent$obj_radius +2
    #y1 <- shifted_center[2] - obj1_extent$obj_radius -2
    #y2 <- shifted_center[2] + obj1_extent$obj_radius +2

    search_radius <-5

    x1 <- shifted_center[1] -search_radius
    x2 <- shifted_center[1] +search_radius
    y1 <- shifted_center[2] -search_radius
    y2 <- shifted_center[2] +search_radius

    return(list(x1=x1, x2=x2, y1=y1, y2=y2, center_pred=shifted_center))
}

euclidean_dist <- function(vec1, vec2){
    sqrt(sum((vec1-vec2)^2))
}

