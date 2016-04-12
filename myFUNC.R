#! This is my collection of easy to use functions in 'R'

#====================================to insert the string at a specific location
split_str_by_index <- function(target, index) {
    index <- sort(index)
    substr(rep(target, length(index) + 1),
           start = c(1, index),
           stop = c(index -1, nchar(target)))
}

#Taken from https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
interleave <- function(v1,v2)
{
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}

insert_str <- function(target, insert, index) {
    insert <- insert[order(index)]
    index <- sort(index)
    paste(interleave(split_str_by_index(target, index), insert), collapse="")
}
#-------------------------------------------------------------------------------




# ============ START ======NCDF4 functions=========== START ============#
#get start and count vectors for given lat-lon values. It also adds the time dim for default reading.
get.start_count<-function(ncfile, lons, lats){
    lon_ind <- start_count(ncfile, lons, 1)
    lat_ind <- start_count(ncfile, lats, 2)
    out<-list(c(lon_ind[1], lat_ind[1], 1), c(lon_ind[2], lat_ind[2], -1))
    return(out)
}

#this function finds start and count indices for lat and lon
start_count<-function(ncfile, reqVals, dim){
    dim_name=ncfile$dim[[dim]]$name             #dim name in the file
    inVals<-ncvar_get(ncfile, varid = dim_name) #dim values in the file
    ind1<-which(abs(inVals-reqVals[1])==min(abs(inVals-reqVals[1])))
    ind2<-which(abs(inVals-reqVals[2])==min(abs(inVals-reqVals[2])))
    start<-min(ind1, ind2)
    end<-max(ind1, ind2)
    count<-end-start
    return(c(start, count))
}
# ============ END ============NCDF4 functions========== END ============#

#to find file paths for given expression
fPath<-function(path, pat){
    path<-list.files(path=path, pattern = pat, recursive = T, full.names = T)
    if(length(path)!=1){
        print(path)
        stop(paste("ERROR:", length(path), "file(s) found.") )
    }
    return(path)
}

dbz2rr<-function(dbz){
    zra=240
    zrb=1.5
    ((10^(dbz/10))/zra)^(1/zrb)
}

#to comput wind direction from the U and V components 
windDir <- function(u, v) {
  (180 / pi) * atan(u/v) + ifelse(v>0,180,ifelse(u>0,360,0))
}


