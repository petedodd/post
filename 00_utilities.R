## shared utilities

## libraries
library(ggplot2); theme_set(theme_bw())
library(scales)
library(data.table)
library(imputeTS)


Here <- function(x) gsub("([^\\/\\.]+\\/\\.\\.\\/)","",here(x))


absspace <- function(x,...) {             #works
   format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

## looking at big numbers
see <- function(x,ns=3)formatC(signif(x,ns),big.mark = ",",format='fg') #for reading big numbers
## unc for  X*Y
xfun <- function(x,y,x.sd,y.sd) sqrt(x^2 * y.sd^2 + x.sd^2 * y.sd^2 + y^2 * x.sd^2)
## unc for X*Y*Z
xfun3 <- function(x,y,z,x.sd,y.sd,z.sd) xfun(x, y * z, x.sd, xfun(y,z,y.sd,z.sd) )

## adding uncorrelated errors
Ssum <- function(x,...) sqrt(sum(x^2,...))

## var sum x * y
## errors with corr(Xi,Yj) = 0; i !=j: corr(Xi,Xj) = 1 & corr(Yi,Yj)=1
sumxy <- function(x,y,x.sd,y.sd){
  ans2 <- sum(x * y.sd)^2 + sum(y * x.sd)^2 + sum(y.sd * x.sd)^2
  sqrt(ans2)
}

## var sum x * y
## errors with corr(Xi,Yj) = 0; i !=j: corr(Xi,Xj) = 0 & corr(Yi,Yj)=1
Ssumxy <- function(x,y,x.sd,y.sd){
  ans2 <- sum(x * y.sd)^2 + sum(y^2 * x.sd^2) + sum(y.sd^2 * x.sd^2)
  sqrt(ans2)
}

## var sum x * y
## errors with corr(Xi,Yj) = 0; i !=j: corr(Xi,Xj) = 0 & corr(Yi,Yj)=0
SSumxy <- function(x,y,x.sd,y.sd){
  Ssum(xfun(x,y,x.sd,y.sd))
}



## ---- 3 term items

## var sum x * y * z
## errors with corr(Xi,Yj)=corr(Yi,Zj)=0; i!=j: corr(Xi,Xj)=corr(Yi,Yj)=corr(Zi,Zj)=0
SSSumxyz <- function(x,y,z,x.sd,y.sd,z.sd){
  Ssum(xfun3(x,y,z,x.sd,y.sd,z.sd))
}

## var sum x * y * z
## errors with corr(Xi,Yj)=corr(Yi,Zj)=0; i!=j: corr(Xi,Xj)=0, corr(Yi,Yj)=corr(Zi,Zj)=1
Sssumxyz <- function(x,y,z,x.sd,y.sd,z.sd){
  ans2a <- sum(x * y.sd * z.sd)^2 + sum(x * y * z.sd)^2 + sum(x * z * y.sd)^2
  ans2b <- sum(x.sd^2 * (y^2 * z^2 + y.sd^2 * z^2 + y^2 * z.sd^2 + y.sd^2 * z.sd^2))
  ans2 <- ans2a + ans2b
  sqrt(ans2)
}

## var sum x * y * z
## errors with corr(Xi,Yj)=corr(Yi,Zj)=0; i!=j: corr(Xi,Xj)=corr(Yi,Yj)=corr(Zi,Zj)=1
sssumxyz <- function(x,y,z,x.sd,y.sd,z.sd){
  ans2a <- sum(x * y.sd * z.sd)^2 + sum(x.sd * y * z.sd)^2 + sum(x.sd * z * y.sd)^2
  ans2b <- sum(x * y * z.sd)^2 + sum(x.sd * y * z)^2 + sum(x * z * y.sd)^2
  ans2 <- ans2a + ans2b + sum(x.sd * y.sd * z.sd)^2
  sqrt(ans2)
}


## ## ## ====== checks ======
## xfun(1,1,.1,.1)
## sd(rnorm(1e7,1,.1)*rnorm(1e7,1,.1))
## xfun3(1,1,2,.1,.1,.3)
## sd(rnorm(1e7,1,.1)*rnorm(1e7,1,.1)*rnorm(1e7,2,.3))

## ## block-wise correlations
## z <- rep(1,1e1)

## ## cn is rep
## ## x,y uncor but perfect cor in each unit
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, x=z+rnorm(1)/10, y=z+rnorm(1)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(x*y)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]                          #sqrt(2)
## sumxy(z,z,0.1*z,0.1*z)                 #correct

## ## x,y uncor but perfect cor in each unit for y only
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, x=z+rnorm(10)/10, y=z+rnorm(1)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(x*y)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]                          #1.06
## Ssumxy(z,z,0.1*z,0.1*z)                #1.05 correct

## ## x,y uncor completely
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, x=z+rnorm(10)/10, y=z+rnorm(10)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(x*y)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]                          #0.45
## Ssum(xfun(z,z,0.1*z,0.1*z))            #0.45 correct
## SSumxy(z,z,0.1*z,0.1*z)            #0.45 correct

## ## -------- 3 part tests -----

## ## cn is rep
## ## x uncor, y,z corr but perfect cor in each unit
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, X=z+rnorm(10)/10,
##                                       Y=z+2*rnorm(1)/10,
##                                       Z=2*z+rnorm(1)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(X*Y*Z)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]
## Sssumxyz(z,z,2*z,0.1*z,0.2*z,0.1*z)       #correct



## ## cn is rep
## ## x, y,z corr but perfect cor in each unit
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, X=z+rnorm(1)/10,
##                                       Y=z+2*rnorm(1)/10,
##                                       Z=2*z+rnorm(1)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(X*Y*Z)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]
## sssumxyz(z,z,2*z,0.1*z,0.2*z,0.1*z)       #correct

## ## x,y,z uncor
## zd <- list()
## for(i in 1:1e4) zd[[i]] <- data.table(cn=i, X=z+rnorm(10)/10,
##                                       Y=z+2*rnorm(10)/10,
##                                       Z=2*z+rnorm(10)/10 )
## zd <- rbindlist(zd)
## tmp <- zd[,.(sxy=sum(X*Y*Z)),by=cn]
## tmp[,mean(sxy)]
## tmp[,sd(sxy)]
## SSSumxyz(z,z,2*z,0.1*z,0.2*z,0.1*z)       #correct



## ====== demography turn N into 1 year age group interpolator for right years
acts <- c('04','514','1524','2534','3544','4554','5564','65')
racts <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')
ev5 <- seq(from = 1,to=99,by=5)

## flags
plt <- TRUE
overwrite <- TRUE                      #NOTE doesn't apply to age maps in this file


## age maps
fn <- here('../tmpdata/amap.Rdata')
if(file.exists(fn)){
  load(here('../tmpdata/amap.Rdata'))
  load(here('../tmpdata/lamap.Rdata'))
  load(here('../tmpdata/isokey.Rdata'))
  hivcountries <- scan(here('texto/hivcountries.txt'),what='character')
} else {
  load(here('../indata/N_simple.Rdata'))           #5 year age groups
  ## make an age map between 5year ages and TB ages
  amap <- data.table(AgeGrp=N[,unique(AgeGrp)])
  rracts <- racts[c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8,8,8,8,8,8)]
  amap[,acats:=rracts]
  amap
  save(amap,file=here('../tmpdata/amap.Rdata'))
  ## with HIV data
  hivcountries <- scan(here('texto/hivcountries.txt'),what='character')
  load(file=here('../tmpdata/HH.Rdata'))
  (anmz <- as.character(HH[,unique(age_name)]))
  anmzl <- c(anmz,rep(rev(anmz)[1],4))
  amap[,age_name:=anmzl]
  lamap <- amap[rep(1:20,each=5)]
  lamap[,age:=0:99]
  lamap
  save(lamap,file=here('../tmpdata/lamap.Rdata'))
  ## make isokey
  tmp <- fread(here('../indata/TB_burden_countries_2020-02-24.csv'))
  isokey <- unique(tmp[,.(iso3,g_whoregion)])
  save(isokey,file=here('../tmpdata/isokey.Rdata'))
}

