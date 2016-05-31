#+ echo=FALSE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
library(ncdf4)
library(EBImage)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())


# read radar scan and return image of convective pixels only.
read_scan_cHeight<- function(ncfile, nscan){
    steiner <- ncvar_get(ncfile, varid = "steiner_class", start = c(1, 1, nscan), count = c(-1, -1, 1))
    conv_height <- ncvar_get(ncfile, varid= "zero_dbz_cont", start = c(1, 1, nscan), count = c(-1, -1, 1))
    conv_height <- replace(conv_height, steiner!=2, 0.0)    #remove non-convective
    conv_height <- replace(conv_height, is.na(steiner), 0.0)    #remove NA
}


# Given the convective height image, it returns vertical classification (Cg=1, Cb=2, Co=3)
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


# First table() call returns, How many pixels has same lable?
# This is the size of each echo.
# Second, table() call tells how many n-pixel objects are found?
count_objectPixels <- function(label_img, size_hist){
    size_table <- table(label_img[label_img>0])
    size_count <- table(size_table)
    index_str <- names(size_count)
    for(i in index_str){
        size_hist[as.numeric(i)] <- size_hist[as.numeric(i)]+size_count[i]
    }
    return(size_hist)
}

# For each object in the image, check each pixel and find its height and increament the counter.
count_objectClassPixels <- function(label_img, height_img, size_hist){
    nObj <- max(label_img)

    for(i in seq(nObj)){                                #for each object
        obj_index <- which(label_img==i, arr.ind = TRUE)
        obj_size <- length(obj_index[, 1])
        for(pixel in seq(obj_size)){                    #check each pixel
            pix_ind <- as.vector(obj_index[pixel, ])
            pix_height <- height_img[pix_ind[1], pix_ind[2]]

            size_hist[obj_size, pix_height] <- size_hist[obj_size, pix_height] +1
        }
    }
    return(size_hist)
}

#read radar data
setwd("~/projects/darwin/data/2d/")
ifile_radar <- "./cpol_2D_0405.nc" #a file for a season
incf_radar <- nc_open(ifile_radar)

time <- ncvar_get(incf_radar, varid = "time")
ntime <- length(time)

max_objSize <- 1000

size_hist <- matrix(data=0, nrow = max_objSize, ncol = 3)
echo_size <- rep(0, max_objSize)

for(scan in seq(ntime)){
conv_height <- read_scan_cHeight(incf_radar, scan)
height_lable <- bwlabel(conv_height)
if(max(height_lable)==0) next
conv_class <- get_vertical_class(conv_height)
size_hist <- count_objectClassPixels(height_lable, conv_class, size_hist)
echo_size <- count_objectPixels(height_lable, echo_size)
}

sizes <- seq(1:max_objSize)


pdf("~/projects/darwin/graphs/height_objectSize.pdf")

#' Size of convective echoes in pixels (X-axis) Vs. % frequency of occurrence of echoes (red) and % of total area by echo sizes (blue).
#+ echo=FALSE, eval=TRUE, warning=FALSE, error=FALSE, message=FALSE
echo_size_pc <- echo_size*100/sum(echo_size)
echo_area <- echo_size*sizes
echo_area_pc <- echo_area*100/sum(echo_area)

echo_size_df <- data.frame(echo_size_pc, echo_area_pc, sizes)
colnames(echo_size_df) <- c("freq", "area", "size")

echo_size_df <- melt(echo_size_df, id.vars = "size")

ggplot(data=echo_size_df, aes(x=size, y=value, color=variable))+geom_line(stat = "identity")+
    coord_cartesian(xlim=c(0, 150)) +scale_y_log10() +annotation_logticks(sides="l")


#'Size of convective echoes in pixels (X-axis) Vs. Percentage of pixels of Type Cg, Cb, Co in these echoes.
#+ echo=FALSE, eval=TRUE, warning=FALSE, error=FALSE, message=FALSE
size_hist_pc <- size_hist*100/sum(size_hist)
size_hist_all <- apply(size_hist_pc, MARGIN = 1, FUN = sum)

size_hist_df<-data.frame(size_hist_pc, size_hist_all, sizes)
colnames(size_hist_df)<-c("Cg", "Cb", "Co", "All", "size")

size_hist_df <- melt(size_hist_df, id.vars = "size", variable.name = "Type", value.name = "Pixels")

ggplot(size_hist_df, aes(x = size, y = Pixels, color=Type)) +
    geom_line(stat='identity')+coord_cartesian(xlim=c(0, 150), ylim=c(0.001, 100)) +
    scale_y_log10() + annotation_logticks(sides="l")

#-------------------------------------------------plotting percentage in ggplot2
#' Plot shows percentage of Cg, Cb, and Co, pixels that constituted the echoes of
#' various sizes.
#+ echo=FALSE, eval=TRUE, warning=FALSE, error=FALSE, message=FALSE
percentage <- function(data){
    return(data*100/sum(data, na.rm = TRUE))
}

size_hist_pc<-apply(size_hist, MARGIN =1, FUN = percentage)
size_hist_pc <- t(size_hist_pc)


size_hist_df<-data.frame(size_hist_pc, sizes)
colnames(size_hist_df)<-c("Cg", "Cb", "Co", "size")

size_hist_df <- melt(size_hist_df, id.vars = "size", variable.name = "Type", value.name = "Frequency")

ggplot(size_hist_df, aes(x = size, y = Frequency, fill=Type)) +
    geom_bar(stat='identity')+coord_cartesian(xlim=c(0, 150))

dev.off()




