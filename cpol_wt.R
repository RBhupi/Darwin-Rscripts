# Use wavelet transform to extract cumuliform clouds from the radar data.




#Function computes WT of the 2d array up to given max_scale.
#Negative WT are removed. Not tested for non-square data.
wavelet <- function (data, max_scale) {
    dim_data <- dim(data)
    nx <- dim_data[1]
    ny <- dim_data[2]

    sf=c(0.0625, 0.25, 0.375) #scaling coef
    wt <- array(dim = c(max_scale, nx, ny))
    temp1<-array(dim = dim(data))
    temp2<-array(dim = dim(data))

    #start Wavelet loop
    for(scale in 1:max_scale){
        x1=2^(scale-1)
        x2=2 * x1

        #Row-wise (longitude) smoothing
        for (i in 1:nlon){

            #find the indices for prev and next points on the line
            prev2=abs(i-x2)
            prev1=abs(i-x1)
            next1=(i+x1)
            next2=(i+x2)

            #If these indices are outside the image use next values
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }

            if(next1 > nlon){
                next1 = (2 * nlon) - next1
            }

            if(next2 > nlon){
                next2 = (2 * nlon) - next2
            }

            for (j in 1:(nlat)) {
                l2 = data[j, prev2]
                l1 = data[j, prev1]
                r1 = data[j, next1]
                r2 = data[j, next2]
                temp1[j, i] = sf[1] * (l2+r2) + sf[2] * (l1 + r1) + sf[3] * data[j, i]
            }
        }


        #column-wise (latitude) smoothing
        for(i in 1:nlat){

            prev2=abs(i-x2)
            prev1=abs(i-x1)
            next1=(i+x1)
            next2=(i+x2)

            #If these indices are outside the image use next values
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }

            if(next1 > nlat){
                next1 = (2 * nlat) - next1
            }

            if(next2 > nlat){
                next2 = (2 * nlat) - next2
            }

            for(j in 1:nlon){
                t2 = temp1[prev2, j]
                t1 = temp1[prev1, j]
                b1 = temp1[next1, j]
                b2 = temp1[next2, j]
                temp2[i, j] = sf[1] * (t2+b2) + sf[2] * (t1 + b1) + sf[3] * temp1[i, j]
            }
        }

        wt[scale, , ]=data-temp2
        data=temp2


    }
    wt<-replace(wt, wt<0, 0)
    invisible(wt)
}




colors <- c("#FFE4F3","#FFDEF8","#FFD9FE","#FFD4FF","#FFD1FF","#FFCEFF",
            "#FFCBFF","#FFC9FF","#F3C8FF","#E4C8FF","#D3C8FF","#C0C8FF",
            "#ABC8FF","#95C9FF","#7CC9FB","#61C9F2","#3EC8E8","#00C7DC",
            "#00C6D0","#00C4C2","#00C1B3","#00BEA4","#00BA94","#00B583",
            "#00B171","#29AB5E","#45A649","#58A031","#679904","#749300",
            "#7E8C00","#878400","#8F7D00","#967500","#9B6D00","#A06500",
            "#A35D00","#A6540C","#A84A2B","#A9403E")



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

#select data scan
level <- 6
scan <- 1
data <- dbz_bringi[, , level, scan]
data[is.na(data)] <- 0

#Take wavelet transform
scale_max=5
wt <- wavelet(data, max_scale = scale_max)

#make pdf device
pdf("dBZ_WT.pdf", width=6, height=9)
par(mfrow=c(3, 2))
subtitle <- paste(z[level], "Km", "; ", strftime(time_posix[scan], tz = "UTC", usetz = TRUE), sep="")

#plot Reflectivity in First panel
image2D(z = data, x = x, y = y, col=colors,
        xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey")
title(main="Bringi Reflectivity [dBZ]")
text(labels = subtitle, x = min(x), y=(max(y) - 0.1*max(y)), pos = 4)

#and now all WT scales
for(scale in 1:scale_max) {
    image2D(z = wt[scale, , ], x = x, y = y, col=colors,
            xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey")
    title(main=paste("WT scale", scale))
}

dev.off()
