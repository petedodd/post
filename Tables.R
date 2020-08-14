## this will be the tables
rm(list=ls())
library(here)
library(data.table)

whozg <- c('AFR','AMR','EMR','EUR','SEA','WPR')
whozt <- c('Africa','The Americas','Eastern Mediterranean','Europe','South-East Asia',
           'Western Pacific')
wrk <- data.table(g_whoregion=whozg,name=whozt)

ftb <- Vectorize(function(x){
  ## formatter according to GTB rounding rules
  stopifnot (!is.character(x))
  if (!is.na(x)){
    smallpos <- x>0 & x<0.01
    dg <- ifelse(abs(x)>0.01 & abs(x)<100, 2, 3)
    x <- signif(x, 3)
    x <- format(x, digits=dg, nsmall=0, big.mark=" ", justify='right', drop0trailing = T)
    if (smallpos) x <- '<0.01'
  } else x <- '-'
  return(x)}, 'x')


Ssum <- function(x,...) sqrt(sum(x^2,...))

## Function for rounding anything left with a decimal point and removing '<'
## (for number formatting for tables)
zapdotlt <- function(x){
  x <- gsub("<","",x)
  got <- grep("\\.",x)
  x[got] <- as.character(round(as.numeric(x[got])))
  x
}

## Function for formatting big numbers >=1e6, which are breaking with ftb
fmtbig <- function(x,sf=3){
  x <- signif(x, sf)
  x <- format(x, digits=3,
              nsmall=0, big.mark=" ", justify='right',
              drop0trailing = TRUE,scientific=FALSE)
  ## x <- stringi::stri_enc_toutf8(x)      #ensure UTF8
  x
}
fmtbig(123451000)
fmtbig(123451000,sf=4)


## ============ Table 1
load(here('../figdat/t1r1.Rdata'))         #Total new tuberculosis cases 1980-2019
load(here('../figdat/t1r2.Rdata'))         #New tuberculosis cases treated 1980-2019
load(here('../figdat/t1r3.Rdata'))         #New tuberculosis cases not treated 1980-2019
## load(here('figdat/t1r4.Rdata')) # = 5 + 6 #Total tuberculosis survivors alive in 2020
load(here('../figdat/t1r5.Rdata'))         #Treated tuberculosis survivors alive in 2020
load(here('../figdat/t1r6.Rdata'))       #Untreated tuberculosis survivors alive in 2020
## load(here('figdat/t1r7.Rdata')) # = 8 + 9 #Life-years post-tuberculosis
load(here('../figdat/t1r8.Rdata'))         #Life-years post-tuberculosis lived by treated
load(here('../figdat/t1r9.Rdata'))       #Life-years post-tuberculosis lived by untreated

t1r4 <- merge(t1r5,t1r6,by='g_whoregion')
t1r7 <- merge(t1r8,t1r9,by='g_whoregion')
t1r4[,value:=value.x+value.y]; t1r7[,value:=value.x+value.y]
t1r4[,value.sd:=sqrt(value.sd.x^2+value.sd.y^2)];
t1r7[,value.sd:=sqrt(value.sd.x^2+value.sd.y^2)]
t1r4[,quantity:='tot2020']; t1r7[,quantity:='LY2020'];
t1r4[,c('quantity.x','quantity.y'):=NULL];t1r7[,c('quantity.x','quantity.y'):=NULL];
t1r4[,c('value.x','value.y'):=NULL];t1r4[,c('value.sd.x','value.sd.y'):=NULL];
t1r7[,c('value.x','value.y'):=NULL];t1r7[,c('value.sd.x','value.sd.y'):=NULL];

## add row numbers
t1r1[,rn:=1]; t1r2[,rn:=2]; t1r3[,rn:=3]
t1r4[,rn:=4]; t1r5[,rn:=5]; t1r6[,rn:=6]
t1r7[,rn:=7]; t1r8[,rn:=8]; t1r9[,rn:=9]

## 2 missing SD
## t1r2[,value.sd:=0]

t1 <- rbindlist(list(t1r1,t1r2,t1r3,t1r4,t1r5,t1r6,t1r7,t1r8,t1r9),use.names = TRUE)

namekey <- data.table(rn=1:9,
                      nm=c('Total new tuberculosis cases 1980-2019',
                           'New tuberculosis cases treated 1980-2019',
                           'New tuberculosis cases not treated 1980-2019',
                           'Total tuberculosis survivors alive in 2020',
                           'Treated tuberculosis survivors alive in 2020',
                           'Untreated tuberculosis survivors alive in 2020',
                           'Life-years lived by tuberculosis survivors 1980-2020',
                         'Life-years lived by treated tuberculosis survivors 1980-2020',
                         'Life-years lived by untreated tuberculosis survivors 1980-2020'
                           ))

t1[,.N,by=quantity]

if(! 'Global' %in% wrk$g_whoregion) wrk <- rbind(wrk,data.table(g_whoregion='Global',name='Global'))
t1 <- merge(t1,namekey,by='rn')
t1 <- merge(t1,wrk,by='g_whoregion')

t1$name <- factor(t1$name,levels=wrk$name,ordered=TRUE)
t1$nm <- factor(t1$nm,levels=namekey$nm,ordered=TRUE)


t1 <- dcast(t1,nm ~ name,value.var = 'value')
T1


tmp <- T1[,lapply(.SD,fmtbig),.SDcols=2:ncol(T1)]
T1F <- cbind(quantity=T1[,nm],tmp)

fwrite(T1F,file=here('figs/table1.csv'))

## newversion with unc
T1.sd <- dcast(t1,nm ~ name,value.var = 'value.sd')
T1.sd

nmz <- names(T1)[2:ncol(T1)]

T1.sd <- as.matrix(T1.sd[,..nmz])
T1.mid <- as.matrix(T1[,..nmz])

T1.lo <- T1.mid - T1.sd*1.96
T1.lo[T1.lo<0] <- 0
T1.hi <- T1.mid + T1.sd*1.96

SF <- 4
T1.lo <- fmtbig(T1.lo,SF)
T1.hi <- fmtbig(T1.hi,SF)

BB <- matrix(paste0("(",T1.lo," - ",T1.hi,")"),nrow=nrow(T1.lo),ncol=ncol(T1.lo))
BB2 <- fmtbig(T1.mid,SF)
BB2 <- matrix(paste0(BB2,"  ",BB),nrow=nrow(T1.lo),ncol=ncol(T1.lo))

rownames(BB2) <- T1$nm
colnames(BB2) <- names(T1)[-1]
BB2
BB2[,7]

write.csv(BB2,file=here('figs/table1u.csv'))


## ======== Table 2 ===========
## see Results
load(here("../figdat/c1.Rdata")); load(here("../figdat/c1b.Rdata"))
load(here("../figdat/c2.Rdata")); load(here("../figdat/c2b.Rdata"))
load(here("../figdat/c4.Rdata")); load(here("../figdat/c4b.Rdata"))


## combine
c1[,totf:=fmtbig(total,SF)]
c1[,totf.lo:=fmtbig(total-total.sd*1.96,SF)]; c1[,totf.hi:=fmtbig(total+total.sd*1.96,SF)]
c1[,totf:=paste0(totf," (",totf.lo," - ",totf.hi,")")]
c2[,pcr:=round(pc)]

## c4[,age:=paste0(round(wm),"  (",round(wsd),")")]
c4[,age:=round(pck,1)]

c1b[,totf:=fmtbig(total,SF)]
c1b[,totf.lo:=fmtbig(total-total.sd*1.96,SF)]; c1b[,totf.hi:=fmtbig(total+total.sd*1.96,SF)]
c1b[,totf:=paste0(totf,"  (",totf.lo," - ",totf.hi,")")]
c2b[,pcr:=round(pc)]

## c4b[,age:=paste0(round(wm)," (",round(wsd),")")]
c4b[,age:=round(pck,1)]

tab2 <- merge(c1[,.(g_whoregion,total5=totf)],c2[,.(g_whoregion,pcr5=pcr)])
tab2 <- merge(tab2,c4[,.(g_whoregion,age5=age)])
tab2 <- merge(tab2,c1b[,.(g_whoregion,total2=totf)])
tab2 <- merge(tab2,c2b[,.(g_whoregion,pcr2=pcr)])
tab2 <- merge(tab2,c4b[,.(g_whoregion,age2=age)])
setkey(tab2,g_whoregion)
tab2 <- tab2[c('AFR','AMR','EMR','EUR','SEA','WPR','Global')] #reorder

write.csv(tab2,file=here::here('figs/table2.csv'))
