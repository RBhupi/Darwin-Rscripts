#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads netCDF files obtained from Benjamine and
# estimates flow vectors, using phase correlation with FFT (Leese et al 1971).
#
#@todo :
#==========================================================================================






library(EBImage)
library(ncdf4)
library(plot3D)
library(spatstat) #for smoothing


## Reads a single frame from netcdf file, replaces NA with 'zeros'
read_ncFrame <- function(ncfile, var_name, frame_num) {
    start_vec <- c(1, 1, frame_num)
    count_vec <- c(-1, -1, 1)

    data <- ncvar_get(ncfile, varid = var_name, start = start_vec, count = count_vec)
    data <- replace(data, is.na(data), 0.0)
    invisible(data)
}

#comput cross-covariance using FFT
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
    invisible(crossCov)
}

## Rearrange the crossCov matrix so that 'zero' frequency or DC component is in the middle of the matrix.
#   This function is adopted from following discussion on stackOverflow
#   http://stackoverflow.com/questions/30630632/performing-a-phase-correlation-with-fft-in-r
fft_shift <- function(x) {
    if(class(x)=='matrix') {
        rd2 <- floor(nrow(x)/2)
        cd2 <- floor(ncol(x)/2)

        ## Identify the first, second, third, and fourth quadrants
        q1 <- x[1:rd2,1:cd2]
        q2 <- x[1:rd2,(cd2+1):ncol(x)]
        q3 <- x[(rd2+1):nrow(x),(cd2+1):ncol(x)]
        q4 <- x[(rd2+1):nrow(x),1:cd2]

        ## rearrange the quadrants
        centered.t <- rbind(q4,q1)
        centered.b <- rbind(q3,q2)
        centered <- cbind(centered.b,centered.t)

        return(Re(centered))
    }
    if(class(x)!='matrix') {
        print('sorry, this class of input x is not supported yet')
    }
}

## Estimates flow vectors in two images
fft_flowVectors <- function (im1, im2) {
    crossCov <- fft_crossCov(im1, im2)
    cov_shifted <- fft_shift(crossCov)
    cov_smooth <- blur(as.im(cov_shifted))

    dims<-dim(im1)

    pshift <- which(cov_smooth$v==max(cov_smooth$v),arr.ind=TRUE)
    pshift <- pshift-(dims[1]/2)

    return(c(pshift[1], pshift[2]))
}


#set directory, Open file
setwd("~/data/darwin_radar/2d/")
infile_name <- "./cpol_2D_0405.nc"
ncfile <- nc_open(infile_name)

#read dimention variables
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
time <- ncvar_get(ncfile, varid = "time")

boxLength <- 11 #pixels
nbox_xdim<-11 #number of boxes in x dimention

#read first frame
img1 <- read_ncFrame(ncfile, var_name = "rain_rate", 1)
img1_tiles <- untile(img1, nim = c(boxLength, boxLength), lwd = 0) #untile it in to  boxes

#read next frame
img2 <- read_ncFrame(ncfile, var_name = "rain_rate", 2)
img2_tiles <- untile(img2, nim = c(boxLength, boxLength), lwd = 0) #untile it in to  boxes
img_dims <- dim(img2_tiles)

#save headings here




headings <- list()
for(tile in 1:img_dims[3]){
im1 <- img1_tiles[, , tile]
im2 <- img2_tiles[, , tile]

if(max(im1)==0 || max(im2)==0) {
    headings<-append(headings, list(c(NA, NA)))
    next
}
headings <- append(headings, list(fft_flowVectors(im1, im2)))
}

x_headings<-matrix(data = NA, nrow = nbox_xdim, ncol = nbox_xdim)
y_headings<-matrix(data = NA, nrow = nbox_xdim, ncol = nbox_xdim)
for (rdim in 1:nbox_xdim){
    for (cdim in 1:nbox_xdim) {
        x_headings[cdim, rdim] = headings[[(rdim-1)*nbox_xdim + cdim]][1]
        y_headings[cdim, rdim] = headings[[(rdim-1)*nbox_xdim + cdim]][2]

   }
}

image2D(img1, x=x, y=y)
image2D(img2, x=x, y=y)
x_vec <- x[seq(1, 121, by=12)]
y_vec <- y[seq(1, 121, by=12)]

image2D(x_headings, x=x_vec, y=y_vec)
image2D(y_headings, x=x_vec, y=y_vec)
