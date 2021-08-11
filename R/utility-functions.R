####################################################################################
# Utility functions for working with raster files in R
# Author: Clement Bataille, 2020
####################################################################################

## function returns raster of posterior probabilities for univariate normal data
calcCellProb <- function(x,isoscape,std){
  m <- getValues(isoscape)
  s <- getValues(std) # use this if you have a raster of standard errors/deviations
  m <- dnorm(x,mean=m,sd=s)
  cell.dens <- setValues(isoscape,m)
  return(cell.dens)
}

## function returns raster of posterior probabilities for bivariate normal data
## x is the unknown tissue of interest, will have two values, one for each isotope
## m is a 2-D vector, all the values in the raster for each isotope
## v is the same as m, but for variances
## r is a single number - the covariance. Can be vector if estimated as non-stationary
## ras is a raster that will serve as a template for the final product
calcCellProb2D <- function(x,m,v,r,ras) {
  pd <- 1/(2*pi*sqrt(v[,1])*sqrt(v[,2])*sqrt(1-r^2))*exp(-(1/(2*(1-r^2)))*
                                                           ((x[1]-m[,1])^2/v[,1]+(x[2]-m[,2])^2/v[,2]-(2*r*(x[1]-m[,1])*
                                                                                                         (x[2]-m[,2]))/(sqrt(v[,1])*sqrt(v[,2]))))
  pdras <- setValues(ras,pd)
  return(pdras)
}

## function returns raster of posterior probability distribution
calcPostProb <- function(x){
  pp <- x/cellStats(x,sum)
  return(pp)
}

## function returns raster identifying upper quantiles across raster values
qtlRaster <- function(ras,q=0.95) {
  qras <- ras >= quantile(ras,probs=q)
  return(qras)
}

## function returns raster of log likelihoods for each cell
llRaster <- function(x,m,v,r,ras) {
  n <- nrow(x)
  k <- ncol(x)
  for (i in 1:nrow(m)) {
    mu <- m[i,]
    cv <- matrix(v[i,1],r,v[i,2],r,nrow=2)
    for (s in 1:n) {
      xs <- x[s,]
      kernSum  <- crossprod(x,solve(cv))%*%x
    }
    ll[i] <- - (n * k / 2) * (log(2*pi)) - (n/2) * log(abs(cv)) - kernSum/2
  }
  llras <- setValues(ras,ll)
  return(llras)
}

## function to get normalized cell probabilites
calcNormProb <- function(x){
  np <- x/cellStats(x,max)
  return(np)
}

## function to convert degrees (DD) to radians
deg2rad <- function(deg) return(deg*pi/180)

## function to get geodesic distance between two points using Spherical Law of Cosines
## xy specified by radian latitude/longitude 
calcDist <- function(long1, lat1, long2, lat2) {
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  R <- 6371 # earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # distance in km
}
