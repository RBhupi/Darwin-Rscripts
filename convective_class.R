library(ncdf4)
library(plot3D)
library(animation)
library(stringr)
#-------------------------------------------------------------------------------
setwd("/home/bhupendra/data/darwin_radar/2d")
flist<-Sys.glob(path="*.nc")

ncfile <- nc_open(filename = flist[2])

#read x, y, z, and time for reference.
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
z <- ncvar_get(ncfile, varid = "z")
time <- ncvar_get(ncfile, varid = "time")
time_posix <- as.POSIXct(time, origin = "2004-01-01", tz="UTC")

nimage <- 144
start_img <- 15600

#read data
zero_dbz_cont<- ncvar_get(ncfile, varid = "zero_dbz_cont", start = c(1, 1, start_img), count=c(-1, -1, nimage))
steiner_class <- ncvar_get(ncfile, vari= "steiner_class", start = c(1, 1, start_img), count=c(-1, -1, nimage))
timeOfScan<-as.POSIXct(ncvar_get(nc = ncfile, varid = "time", start = c(start_img), count=c(nimage)),
                       origin = "2004-01-01", tz="UTC")

#replace all non-convective pixels by 0.0
conv_height <- replace(zero_dbz_cont, steiner_class!=2, 0.0)
rm(list=c("zero_dbz_cont", "steiner_class"))

hst <- hist(conv_height[conv_height>0], breaks = (seq(5, 40, 1)), plot=F)

pdf(str_replace(string = flist[1], pattern = ".nc", ".pdf"))
plot(hst$mids/2, hst$counts, type = "b", main="steiner convection echo-top height distribution ", cex.main = 0.5,
     xlab="Zero dBZ Height [Km]", ylab="Pixel Count")
grid()
dev.off()

#min max heigths for classification
min_level <- c(5, 15, 31)
max_level <- c(14, 30, 40)

conv_class <- conv_height
for(i in 1:length(min_level)){
conv_class <- replace(conv_class, conv_class>=min_level[i] & conv_class<= max_level[i], i)
}

#make GIF
saveGIF({
    for(i in 1:80){
        image2D(z=conv_class[, , i], x=x, y=y, breaks =1:5)
        title(main="classified convection")
        text(x =-150, y=140, labels = timeOfScan[i], pos = 4)
        grid()
    }
}, movie.name = str_replace(string = flist[1], pattern = ".nc", ".gif"), interval = 0.4)

