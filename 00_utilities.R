## shared utilities

## libraries
library(ggplot2); theme_set(theme_bw())
library(scales)
library(data.table)
library(imputeTS)


absspace <- function(x,...) {             #works
   format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

see <- function(x,ns=3)formatC(signif(x,ns),big.mark = ",",format='fg') #for reading big numbers
xfun <- function(x,y,x.sd,y.sd) sqrt(x^2*y.sd^2+x.sd^2*y.sd^2+y^2*x.sd^2)
Ssum <- function(x,...) sqrt(sum(x^2,...))

## ## check
## xfun(1,1,.1,.1)
## sd(rnorm(1e7,1,.1)*rnorm(1e7,1,.1))


## ====== demography turn N into 1 year age group interpolator for right years
acts <- c('04','514','1524','2534','3544','4554','5564','65')
racts <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')
ev5 <- seq(from = 1,to=99,by=5)

## flags
plt <- FALSE
overwrite <- FALSE                      #NOTE doesn't apply to age maps in this file


## age maps
fn <- here('../tmpdata/amap.Rdata')
if(file.exists(fn)){
  load(here('../tmpdata/amap.Rdata'))
  load(here('../tmpdata/lamap.Rdata'))
  load(here('../tmpdata/isokey.Rdata'))
  hivcountries <- scan(here('../tmpdata/hivcountries.txt'),what='character')
} else {
  load(here('../indata/N_simple.Rdata'))           #5 year age groups
  ## make an age map between 5year ages and TB ages
  amap <- data.table(AgeGrp=N[,unique(AgeGrp)])
  rracts <- racts[c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8,8,8,8,8,8)]
  amap[,acats:=rracts]
  amap
  save(amap,file=here('../tmpdata/amap.Rdata'))
  ## with HIV data
  hivcountries <- scan(here('../tmpdata/hivcountries.txt'),what='character')
  load(file=here('../tmpdata/HH.Rdata'))
  (anmz <- as.character(HH[,unique(age_name)]))
  anmzl <- c(anmz,rep(rev(anmz)[1],4))
  amap[,age_name:=anmzl]
  lamap <- amap[rep(1:20,each=5)]
  lamap[,age:=0:99]
  lamap
  save(lamap,here('../tmpdata/lamap.Rdata'))
  ## make isokey
  tmp <- fread(here('../indata/TB_burden_countries_2020-02-24.csv'))
  isokey <- unique(tmp[,.(iso3,g_whoregion)])
  save(isokey,file=here('../tmpdata/isokey.Rdata'))
}

## from https://rdrr.io/github/hadley/bigvis/src/R/weighted-stats.r
#' Compute quantiles of weighted data.
#'
#' @details
#' Currently only implements the type 7 algorithm, as described in
#' \code{\link{quantile}}. Based on \code{\link{quantile}} written by R-core.
#'
#' @inheritParams weighted.var
#' @param probs numeric vector of probabilities between 0 and 1
#' @param na.rm If \code{TRUE} will automatically remove missing values
#'   in \code{x} or \code{w}.
#' @family weighted statistics
#' @export
#' @examples
#' x <- runif(200)
#' w <- rpois(200, 5) + 1
#' weighted.quantile(x, w)
weighted.quantile <- function (x, w, probs = seq(0, 1, 0.25), na.rm = FALSE) {
  stopifnot(length(x) == length(w))
  na <- is.na(x) | is.na(w)
  if (any(na)) {
    if (na.rm) {
      x <- x[!na]
      w <- w[!na]
    } else {
      stop("Missing values not allowed when na.rm is FALSE", call. = FALSE)
    }
  }

  # Ensure x and w in ascending order of x
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  # Find closest x just below and above index
  n <- sum(w)
  index <- 1 + (n - 1) * probs
  j <- floor(index)

  wts <- cumsum(w)
  lo <- x[lowerBound(j, wts)]     # X_j
  hi <- x[lowerBound(j + 1, wts)]

  g <- index - j
  ifelse(lo == hi, lo, (1 - g) * lo + g * hi)
}
