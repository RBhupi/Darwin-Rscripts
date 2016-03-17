library(ncdf4)
library(plot3D)
library(animation)

colors <- c("#FFDBEA","#FFD4C3","#FCCFA0","#DACB87","#B2C67C","#86C080",
            "#55B88D","#15AC9B","#009EA6","#278DAD","#5578AC")

#read reflectivity and other data from radar data
setwd("/home/bhupendra/data/darwin_radar/test/CPOL/outNC/")
inFileName <- "cpol_3d_fields_classes_20050402.nc"
ncfile <- nc_open(inFileName)

#read x, y, z, and time for reference
x <- ncvar_get(ncfile, varid = "x")
y <- ncvar_get(ncfile, varid = "y")
z <- ncvar_get(ncfile, varid = "z")
time <- ncvar_get(ncfile, varid = "time")
time_posix <- as.POSIXct(time, origin = "1970-01-01", tz="UTC")

#read steiner and merhala cloassification and bringi reflectivity and hdrometeor classification
steiner_rainClass <- ncvar_get(ncfile, varid = "cs")
merhala_rainClass <- ncvar_get(ncfile, varid = "cm")
dbz_bringi <- ncvar_get(ncfile, varid = "zh")
dbz_bringi <- replace(dbz_bringi, dbz_bringi<0, NA)
hydro_class <- ncvar_get(ncfile, varid = "hc")


colors <- c("#FFE4F3","#FFDEF8","#FFD9FE","#FFD4FF","#FFD1FF","#FFCEFF",
            "#FFCBFF","#FFC9FF","#F3C8FF","#E4C8FF","#D3C8FF","#C0C8FF",
            "#ABC8FF","#95C9FF","#7CC9FB","#61C9F2","#3EC8E8","#00C7DC",
            "#00C6D0","#00C4C2","#00C1B3","#00BEA4","#00BA94","#00B583",
            "#00B171","#29AB5E","#45A649","#58A031","#679904","#749300",
            "#7E8C00","#878400","#8F7D00","#967500","#9B6D00","#A06500",
            "#A35D00","#A6540C","#A84A2B","#A9403E")

#make GIF
#plot Reflectivity
saveGIF({
for(i in 1:144){
    level <- 6
    scan <- i
    subtitle <- paste(z[level], "Km", "; ", strftime(time_posix[scan], tz = "UTC", usetz = TRUE), sep="")
    image2D(z = dbz_bringi[, , level, scan], x = x, y = y, breaks = seq(0, 60, by=1.5), col = colors,
            xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey")
    title(main="Bringi Reflectivity")
    text(labels = subtitle, x = min(x), y=(max(y) - 0.1*max(y)), pos = 4)
}
}, movie.name = "dbz_motion.gif", interval = 0.2)



pdf("plot2D.pdf")
for(i in seq(2, 15, by = 1)){
level <- i
scan <- 1
subtitle <- paste(z[level], "Km", " ", strftime(time_posix[scan], tz = "UTC", usetz = TRUE), sep="")

par(mfrow = c(2, 2))
#plot Reflectivity
image2D(z = dbz_bringi[, , level, scan], x = x, y = y, xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey")
title(main="Bringi Reflectivity")
text(labels = subtitle, x = min(x), y=(max(y)-0.1* max(y)), pos = 4, cex=0.8)

#plot Classification
colors <- c("white", "#CB7D8D","#67A160","#7494C8")
image2D(z = merhala_rainClass[, , level, scan], x = x, y = y, xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey", col = colors)
title(main="Merhala Classification")

colors <- c("white", "#CB7D8D","#67A160")

image2D(z = steiner_rainClass[, , level, scan], x = x, y = y, xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey", col = colors)
title(main="Steiner Classification")

image2D(z = hydro_class[, , level, scan], x = x, y = y, xlab="Distance from Radar [Km]", ylab="Distance from Radar [Km]", NAcol = "grey", col = colors)
title(main="Hydrometeor Classification")

}
dev.off()







