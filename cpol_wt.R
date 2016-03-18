# Use wavelet transform to extract cumuliform clouds from the radar data.


library(ncdf4)
library(plot3D)
library(misc3d)


#Function computes WT of the 2d array up to given max_scale.
#Negative WT are removed. Not tested for non-square data.
wavelet <- function (data, max_scale) {
    dim_data <- dim(data)
    nlon <- dim_data[1]
    nlat <- dim_data[2]

    sf=c(0.0625, 0.25, 0.375) #scaling coef
    wt <- array(dim = c(max_scale, nlon, nlat))
    temp1<-array(dim = dim(data))
    temp2<-array(dim = dim(data))

    #start Wavelet loop
    for(scale in 1:max_scale){
        x1 <- 2^(scale-1)
        x2 <- 2 * x1

        #Row-wise (longitude) smoothing
        for (i in 1:nlon){

            #find the indices for prev and next points on the line
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)

            #If these indices are outside the image use next values
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }

            if(next1 > nlon | next2 > nlon){
                next1 <- prev1
                next2 <- prev2
            }


            for (j in 1:(nlat)) {
                l2  <-  data[j, prev2]
                l1  <-  data[j, prev1]
                r1  <-  data[j, next1]
                r2  <-  data[j, next2]
                temp1[j, i]  <-  sf[1] * (l2+r2) + sf[2] * (l1 + r1) + sf[3] * data[j, i]
            }
        }


        #column-wise (latitude) smoothing
        for(i in 1:nlat){

            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)

            #If these indices are outside the image use next values
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }

            if(next1 > nlat){
                next1  <-  (2 * nlat) - next1
            }

            if(next2 > nlat){
                next2  <-  (2 * nlat) - next2
            }

            for(j in 1:nlon){
                t2  <-  temp1[prev2, j]
                t1  <-  temp1[prev1, j]
                b1  <-  temp1[next1, j]
                b2  <-  temp1[next2, j]
                temp2[i, j]  <-  sf[1] * (t2+b2) + sf[2] * (t1 + b1) + sf[3] * temp1[i, j]
            }
        }

        wt[scale, , ] <- data-temp2
        data <- temp2


    }
    wt<-replace(wt, wt<0, 0)
    invisible(wt)
}

#return cumulative sum of all WT components
wt_cumsum<-function(data, max_scale){
    wt <- wavelet(data, max_scale)
    wt_sum <- apply(wt, MARGIN = c(2, 3), FUN = sum)
    invisible(wt_sum)
}


colors20 <- c("#FFE4F3","#FFD9FE","#FFD0FF","#FFCBFF","#F0C8FF","#CEC8FF",
            "#A5C9FF","#73C9F8","#2BC8E3","#00C5C9","#00BFAB","#00B78A",
            "#19AD65","#52A239","#709400","#868600","#957700","#9F6600",
            "#A65408","#A9403E")

bg_color <- "gray80"

setwd("/home/bhupendra/data/darwin_radar/test/CPOL/outNC/")
inFileName <- "cpol_3d_fields_classes_20050402.nc"
ncfile <- nc_open(inFileName)

#read x, y, z, and time for reference
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
z <- ncvar_get(ncfile, varid = "z")
nx <- length(x)
ny <- length(y)


time <- ncvar_get(ncfile, varid = "time")
time_posix <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")

#read steiner and merhala cloassification and bringi reflectivity and hdrometeor classification
dbz_bringi <- ncvar_get(ncfile, varid = "zh")
dbz_bringi <- replace(dbz_bringi, dbz_bringi<0, NA)

class_steiner <- ncvar_get(ncfile, varid = "cs")
class_merhala <- ncvar_get(ncfile, varid="cm")
class_hydromet <- ncvar_get(ncfile, varid = "hc")

scan <- 65
scale_max <- 3

#make pdf device
pdf(paste("dBZ_WT3_compare_levels", scan, ".pdf", sep=""), width=9, height=6)
par(mfrow=c(2, 3))

#select data scan
for(level in 1:40){
#level <- 4
data <- dbz_bringi[, , level, scan]
data[is.na(data)] <- 0

#Take wavelet transform

wt <- wt_cumsum(data, max_scale = scale_max)

wt_cleaned <- replace(wt, wt < mean(wt)+2*sd(wt), 0)



subtitle <- paste(z[level], "Km", "; ", strftime(time_posix[scan], tz = "UTC", usetz = TRUE), sep="")
#plot Reflectivity in First panel
colors<-colors20
image2D(z = data, x = x, y = y, col=colors, breaks = seq(0, 60, by=3),
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color)
title(main="Bringi Reflectivity [dBZ]")
text(labels = subtitle, x = min(x), y=(max(y) - 0.1*max(y)), pos = 4)



#plot WT3 raw
colors<-colors20
image2D(z = wt, x = x, y = y, col=colors, breaks = seq(0, 60, by=3),
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color)
size <- 2^(scale_max-1) * 2.5
title(main=paste("Sum of Raw WT (size", size, "Km)"))


#plot WT3 cleaned
colors<-colors20
image2D(z = wt_cleaned, x = x, y = y, col=colors, breaks = seq(0, 60, by=3),
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color)
title(main="WT with noise reduction", cex=0.6)



#plot Classification
colors <- c("white", "#CB7D8D","#67A160","#7494C8")
image2D(z = class_merhala[, , level, scan], x = x, y = y,
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color, col = colors)
title(main="Merhala Classification")

colors <- c("white", "#CB7D8D","#67A160")
image2D(z = class_steiner[, , level, scan], x = x, y = y,
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color, col = colors)
title(main="Steiner Classification")

colors<- c("white", "#E16A86","#D07C42","#B08D00","#7F9C00","#00A742","#00AD84",
           "#00AAB7","#009EDA","#8986E6","#C86DD7")

image2D(z = class_hydromet[, , level, scan], x = x, y = y,
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = bg_color, col = colors)
title(main="Hydrometeor Classification")


}
dev.off()
