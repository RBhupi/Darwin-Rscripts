#!---------------------------------------------------------------------------
# @author Bhupendra Raut (www.baraut.info)
# @brief This R script reads netCDF files obtained from Benjamine and
# estimates flow vectors, using phase correlation with FFT (Leese et al 1971).
#
#Reference :Leese, John A., Charles S. Novak, and Bruce B. Clark.
#           "An automated technique for obtaining cloud motion from geosynchronous
#           satellite data using cross correlation."
#           Journal of applied meteorology 10.1 (1971): 118-132.

#@todo :
# 1. Not tested for the non-square images.
# 3. Not tested for variable tile sizes or numbers
#==========================================================================================
# Start the clock!
start_time <- proc.time()

library(ncdf4)
library(plot3D)
library(spatstat) #for smoothing
library(stringr)
library(EBImage)
library(plyr)

#----------------------------------------------------------function Definitions#

## Reads a single frame from netcdf file, replaces NAs with 'zeros'
read_ncFrame <- function(ncfile, var_name, frame_num) {
    start_vec <- c(1, 1, frame_num)
    count_vec <- c(-1, -1, 1)

    data <- ncvar_get(ncfile, varid = var_name, start = start_vec, count = count_vec)
    data <- replace(data, is.na(data), 0.0)
    invisible(data)
}


## Returns radar frame, with non-convective pixels set to zero.
get_convection_frame <- function(ncfile, var_name, frame_num){
    data <- read_ncFrame(ncfile, var_name = var_name, frame_num)
    steiner <- read_ncFrame(ncfile, var_name = "steiner_class", frame_num)
    data<-replace(data, steiner != 2, 0.0)
}

#computs cross-covariance using FFT
fft_crossCov <- function (img1, img2) {
    fft1_conj <- Conj(fft(img1)) #complex conjugate
    fft2 <- fft(img2)

    C <- (fft2*fft1_conj)/abs(fft2*fft1_conj) #crossCov in Freq domain

    crossCov <- fft(C, inv=TRUE)/length(C)
    crossCov <- Re(crossCov)
    invisible(crossCov)
}

## Rearranges the crossCov matrix so that 'zero' frequency or DC component is in the middle of the matrix.
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

## Takes in two frames and desire tile size to compute flow over tiles of the images.
#  Returns list of x and y headings for each tile.
get_imageFlow <- function(frame1, frame2, tileSize){
    frame1_tiles <- untile(frame1, nim = c(tileSize, tileSize), lwd = 0) #untile it in to  boxes
    frame2_tiles <- untile(frame2, nim = c(tileSize, tileSize), lwd = 0)

    headings <- list()

    num_xytiles <- dim(frame1)/tileSize # num_tiles in x-y direction

    frame_dims <- dim(frame2_tiles)

    for(tile in 1:frame_dims[3]){
        im1 <- frame1_tiles[, , tile]
        im2 <- frame2_tiles[, , tile]

        if(max(im1)==0 || max(im2)==0) {
            headings<-append(headings, list(c(NA, NA)))
            next
        }
        headings <- append(headings, list(fft_flowVectors(im1, im2)))
    }

    x_headings<-matrix(data = NA, nrow = num_xytiles[1], ncol = num_xytiles[2])
    y_headings<-matrix(data = NA, nrow = num_xytiles[1], ncol = num_xytiles[2])
    for (rdim in 1:num_xytiles[1]){
        for (cdim in 1:num_xytiles[2]) {
            x_headings[cdim, rdim] = headings[[(rdim-1)*num_xytiles[1] + cdim]][1]
            y_headings[cdim, rdim] = headings[[(rdim-1)*num_xytiles[2] + cdim]][2]
        }
    }
    invisible(list(x_headings, y_headings))
}


## Returns headings arrays as missing values for writing
get_missing_headings <- function(frame, tileSize){
    frame_tiles <- untile(frame, nim = c(tileSize, tileSize), lwd = 0) #untile it in to  boxes
    num_xytiles <- dim(frame)/tileSize # num_tiles in x-y direction

    x_headings<-matrix(data = NA, nrow = num_xytiles[1], ncol = num_xytiles[1])
    y_headings<-matrix(data = NA, nrow = num_xytiles[2], ncol = num_xytiles[2])

    invisible(list(x_headings, y_headings))
}

## A function to create output netcdf file for flow vectors using input
#   netcdf dimentions and attributes.
create_outNC <- function(ncfile) {
    #read dimention variables
    x <- ncvar_get(ncfile, varid = "x")
    y <- ncvar_get(ncfile, varid = "y")
    time <- ncvar_get(ncfile, varid = "time")

    # x y for flow vectors matrix
    x_vec <- x[seq(ceiling(boxLength/2), 121, by=11)]
    y_vec <- y[seq(ceiling(boxLength/2), 121, by=11)]

    #define dimentions for netcdf
    x_dim <- ncdim_def(name = "x", vals = x_vec, units = ncfile$dim$x$units,
                       longname = "Distance of the center of tiles from Radar")
    y_dim <- ncdim_def(name = "y", vals = y_vec, units = ncfile$dim$y$units,
                       longname = "Distance of the center of tiles from Radar")
    t_dim <- ncdim_def(name = "time", vals = time[-1], units = ncfile$dim$time$units, unlim = TRUE,
                       longname = "time of the second scan used to compute the vectors")

    #define variables
    uVec <- ncvar_def(name="U_Vec", dim = list(x_dim, y_dim, t_dim), units = "pixels per time step",
                      missval = -999.0, compression = 5, prec = "float", longname = "Convection flow vector along x-axis")
    vVec <- ncvar_def(name="V_Vec", dim = list(x_dim, y_dim, t_dim), units = "pixels per time step",
                      missval = -999.0, compression = 5, prec = "float", longname = "Convection flow vector along y-axis")

    #create output file
    outFile <- str_replace(ncfile$filename, pattern = ".nc", replacement = "_fftFlow.nc")
    if(file.exists(outFile)) {
        print(paste("replacing existing file.", outFile))
        file.remove(outFile)
    }
    outNC <- nc_create(outFile, vars = list(uVec, vVec))

    #add attributes
    description <- paste("The CPOL radar field of 121 x 121 pixels was devided in to 121 tiles of 11 x 11 pixels.",
                         "For each tile, the flow vectors were computed using a method described in Leese et. al. (1971).",
                         "Only convective echos as classified by Steiner et al (1995) were used to comput the flow vectors.")

    ncatt_put(outNC, varid = 0, attname = "_description", attval = description, prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_email", attval = "Bhupendra.Raut@monash.edu", prec = "text")
    ncatt_put(outNC, varid = 0, attname = "_date_created", attval = date(), prec = "text")

    invisible(outNC)
}

## main_program does everything for a single input file.
main_program <- function(inFile, var_name, tileSize) {

    ncfile <- nc_open(inFile)

    time_seconds <- ncvar_get(ncfile, varid = "time")
    ntimes <- ncfile$dim$time$len

    outNC <- create_outNC(ncfile)

    print(paste("Computing flow vectors for ", basename(inFile)))
    pb = txtProgressBar(min =2, max = ntimes, initial = 2, style = 3) #progress bar

    #read first frame and call it img2, for convinience
    img2 <- get_convection_frame(ncfile, var_name = ncvar_name, 1)

    missing_scans <- 0
    for(scan in 2:ntimes) { #from second frame
        setTxtProgressBar(pb, scan)

        img1 <- img2
        img2 <- get_convection_frame(ncfile, var_name = ncvar_name, scan)

        #missing scan, if scans are more than 12 minutes apart
        if(time_seconds[scan]-time_seconds[scan-1]>720){
            heads <- get_missing_headings(img1, tileSize)
            missing_scans <- missing_scans + 1
        } else {
            heads <- get_imageFlow(img1, img2, tileSize)
        }

        # Write to the output file
        ncvar_put(nc = outNC, varid = outNC$var$U_Vec$name, vals = heads[[1]],
                  start = c(1, 1, scan-1), count=c(dim(heads[[1]]), 1))
        ncvar_put(nc = outNC, varid = outNC$var$V_Vec$name, vals = heads[[2]],
                  start = c(1, 1, scan-1), count=c(dim(heads[[1]]), 1))
    }
    cat("\n") #new line
    print(paste("missing radar scans in this file = ", missing_scans))

    nc_close(ncfile)
    nc_close(outNC)
}

#----------------------------------------------------------------Calling Program
#initial settings
boxLength <- 11 #in pixels
ncvar_name <- "rain_rate"

#set directory  and get file names
setwd("/home/bhupendra/projects/darwin/data/2d")
flist <- Sys.glob(paths = "./cpol_2D_????.nc")

l_ply(flist, .fun = main_program, ncvar_name, boxLength)

# Stop the clock
print(proc.time() - start_time)

