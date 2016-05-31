library(ncdf4)
library(EBImage)


#' read radar scan and return image of convective pixels only.
read_scan_cHeight<- function(ncfile, nscan){
    steiner <- ncvar_get(ncfile, varid = "steiner_class", start = c(1, 1, nscan), count = c(-1, -1, 1))
    conv_height <- ncvar_get(ncfile, varid= "zero_dbz_cont", start = c(1, 1, nscan), count = c(-1, -1, 1))
    conv_height <- replace(conv_height, steiner!=2, 0.0)    #remove non-convective
    conv_height <- replace(conv_height, is.na(steiner), 0.0)    #remove NA
}


#' First table() call returns, How many pixels has same lable?
#' This is the size of each echo.
#' Second, table() call tells how many n-pixel objects are found?
count_objectPixels <- function(label_img, size_hist){
    size_table <- table(label_img[label_img>0])
    size_count <- table(size_table)
    index_str <- names(size_count)
    for(i in index_str){
        size_hist[as.numeric(i)] <- size_hist[as.numeric(i)]+size_count[i]
    }
    return(size_hist)
}

#read radar data
setwd("~/data/darwin_radar/2d/")
ifile_radar <- "./cpol_2D_0405.nc" #a file for a season
incf_radar <- nc_open(ifile_radar)

time <- ncvar_get(incf_radar, varid = "time")
ntime <- 1000#length(time)

size_hist <- rep(0, 1000)

for(scan in seq(ntime)){
conv_height <- read_scan_cHeight(incf_radar, scan)
height_lable <- bwlabel(conv_height)
if(max(height_lable)==0) next
size_hist <- count_objectPixels(height_lable, size_hist)
}

size_hist_pc <- (size_hist*100)/sum(size_hist)

max_size <- 50
plot(seq(max_size), size_hist_pc[seq(max_size)], log = "y", type="l",
     xlab="size in pixels", ylab="% of Conv. Echoes")
grid()






