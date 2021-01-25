## untreated survivors
rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../figdat/t1r9.Rdata')
  if(file.exists(fn))stop('Not running as figdat/t1r9.Rdata exists!')
}

## load data
load(here('../tmpdata/estl.Rdata'))
load(here('../tmpdata/estg.Rdata'))
load(here('../tmpdata/LA.Rdata'))
load(here('../tmpdata/TBH.Rdata'))      #only for HIV %
load(here('../tmpdata/NNR.Rdata'))      #only for HIV %


## merge against life tables etc
setdiff(unique(estl$iso3),unique(LA$iso3))
estl <- merge(estl,LA,by=c('iso3','year','age','sex'),all.x = TRUE,all.y = FALSE)
estl[,range(year)]
estl

## calculations
estl[,alive.0:=S.0*gapls.0]
estl[,alive.h:=S.h*gapls.h]
estl[,alive:=alive.0+alive.h]
estl[,LYS.0:=LY.0*gapls.0]
estl[,LYS.h:=LY.h*gapls.h]
estl[,LYS:=LYS.h + LYS.0]

## SAs
estl[,alive.2:=(S.2m.0 * gapls.0 + S.2m.h * gapls.h)]
estl[,alive.4:=(S.4m.0 * gapls.0 + S.4m.h * gapls.h)]
estl[,alive.d:=(S.tm0 * gapls.0 + S.tmh * gapls.h)]


## uncertainty
estl[,alive.0.sd:=xfun(S.0,gapls.0,S.0.sd,gapls.0.sd)]
estl[,alive.h.sd:=xfun(S.h,gapls.h,S.h.sd,gapls.h.sd)]
estl[,alive.sd:=sqrt(alive.0.sd^2+alive.h.sd^2)]
estl[,LYS.0.sd:=xfun(LY.0,gapls.0,LY.0.sd,gapls.0.sd)]
estl[,LYS.h.sd:=xfun(LY.h,gapls.h,LY.h.sd,gapls.h.sd)]
estl[,LYS.sd:=sqrt(LYS.h.sd^2 + LYS.0.sd^2)]

## countries without life-table data
estl[!is.finite(alive.sd),sum(gapl)]
1e2*estl[!is.finite(alive.sd),sum(gapl)]/estl[is.finite(alive.sd),sum(gapl)]
estl <- estl[is.finite(alive.sd)]       #drop as tiny


summary(estl[,.(alive.sd,LYS.sd)])
## estl[!is.finite(alive.sd)]              #EGY 1980-1982 in 99 yo
## estl[!is.finite(alive.sd),c('S.sd','LY.sd','S.0.sd','LY.0.sd','alive.sd','LYS.sd','alive.0.sd','LYS.0.sd'):=0] #safety
## summary(N3[,.(alive.sd,LYS.sd)])

## sanity checks
estl[,sum(gapl,na.rm=TRUE)]/1e6
estl[,sum(gapls,na.rm=TRUE)]/1e6
estl[,sum(LYS,na.rm=TRUE)]/1e9
estl[,sum(alive,na.rm=TRUE)]/1e6
estl[,sum(alive.2,na.rm=TRUE)]/1e6
estl[,sum(alive.4,na.rm=TRUE)]/1e6
estl[,sum(alive.d,na.rm=TRUE)]/1e6

## --- plots
## age now! NOTE
estl[,agenow:=age + 2020-year]

estl

## figdat

## === table data
## estg
names(est)
est <- merge(est,NNR[,.(iso3,rat.sd)],by='iso3',all.x=TRUE)

## CDR & rat 
tmpc <- est[,.(value=sum(e_inc_num*(1-rat)),
               value.sd=sumxy(e_inc_num,(1-rat),
               (ocdr.sd/ocdr)*e_inc_num,rat.sd)), #variable-wise corr at country level
           by=.(iso3,g_whoregion)] #fractional unc for val same as rat

t1r1c <- tmpc[,.(iso3,value,value.sd)]
t1r1c[,quantity:='totnew']
save(t1r1c,file=here('../figdat/t1r1c.Rdata'))

t1r1c[,summary(1e2*value.sd/value)]

## regional level
t1r1 <- tmpc[,.(value=sum(value),
              value.sd=Ssum(value.sd)),
             by=g_whoregion]

tmp <- data.table(g_whoregion='Global',
                  value=t1r1[,sum(value)],
                  value.sd=t1r1[,Ssum(value.sd)])
t1r1 <- rbind(t1r1,tmp)
t1r1[,quantity:='totnew']

t1r1[,.(1e2*value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]

save(t1r1,file=here('../figdat/t1r1.Rdata')) #new total

## HIV
est <- merge(est,TBH,by=c('iso3','year'))

est[,value:=hs*e_inc_num*(1-rat)]       #country year

NHu <- est[,.(value=sum(value),
              value.sd=xfun3(hs,e_inc_num,(1-rat),0,(ocdr.sd/ocdr)*e_inc_num,rat.sd)
              ),by=.(iso3,year)]

NHu <- NHu[,.(value=sum(value),value.sd=sum(value.sd)),by=iso3] #cor across iso3
NHu <- NHu[,.(value=sum(value),value.sd=Ssum(value.sd))] #uncor across
NHu[,1e2*value.sd/value]


## NOTE gapl disaggregated over country years in uncorrelated fashion
## unc mainly CDR; neglect rat
tmp <- estl[,.(value=sum(gapl),value.sd=Ssum(gapls.sd)),by=.(iso3,g_whoregion,year)]
tmpc <- tmp[,.(value=sum(value),value.sd=sum(value.sd)),by=.(iso3,g_whoregion)] #NOTE cor iso3

t1r3c <- tmpc[,.(iso3,value,value.sd)]
t1r3c[,quantity:='totnotx']
save(t1r3c,file=here('../figdat/t1r3c.Rdata'))

t1r3c[,summary(1e2*value.sd/value)]

## regional level
t1r3 <- tmpc[,.(value=sum(value),
              value.sd=Ssum(value.sd)),
             by=g_whoregion]

tmp <- data.table(g_whoregion='Global',
                  value=t1r3[,sum(value)],
                  value.sd=t1r3[,Ssum(value.sd)])
t1r3 <- rbind(t1r3,tmp)
t1r3[,quantity:='totnotx']

t1r3[,.(1e2*value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]

save(t1r3,file=here('../figdat/t1r3.Rdata')) #new untreated

## NOTE age/sex disaggregation assumes no correlation <- gapl.sd by age sex
## alive = S * gapls = S * gapl * (1-CFR)
tmp <- estl[,.(value=sum(alive),
               value.2=sum(alive.2),value.4=sum(alive.4),value.d=sum(alive.d),
               value.sd=Ssum(alive.sd)),by=.(iso3,g_whoregion,year)]
tmpc <- tmp[,.(value=sum(value),
               value.2=sum(value.2),value.4=sum(value.4),value.d=sum(value.d),
               value.sd=sum(value.sd)),by=.(iso3,g_whoregion)] #assumed perfect correlation

t1r6c <- tmpc[,.(iso3,value,value.sd)]
t1r6c[,quantity:='totnotx2020']
save(t1r6c,file=here('../figdat/t1r6c.Rdata'))

t1r6c[,summary(1e2*value.sd/value)]

## regional level
t1r6 <- tmpc[,.(value=sum(value),
                value.2=sum(value.2),value.4=sum(value.4),value.d=sum(value.d),
                value.sd=Ssum(value.sd)),
             by=g_whoregion]

tmp <- data.table(g_whoregion='Global',
                  value=t1r6[,sum(value)],
                  value.2=sum(t1r6$value.2),value.4=sum(t1r6$value.4),
                  value.d=sum(t1r6$value.d),
                  value.sd=t1r6[,Ssum(value.sd)])
t1r6 <- rbind(t1r6,tmp)
t1r6[,quantity:='totnew']

t1r6[,.(1e2*value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]


save(t1r6,file=here('../figdat/t1r6.Rdata')) #untreated survivors

## LYS
## untreated LYS
tmp <- estl[,.(value=sum(LYS),value.sd=Ssum(LYS.sd)),by=.(iso3,g_whoregion,year)] #
tmpc <- tmp[,.(value=sum(value),value.sd=sum(value.sd)),by=.(iso3,g_whoregion)] #
t1r9c <- tmpc[,.(iso3,value,value.sd)]
t1r9c[,quantity:='LYnewnotx2020']

save(t1r9c,file=here('../figdat/t1r9c.Rdata'))

t1r9c[,summary(1e2*value.sd/(value+1e-10))]

t1r9 <- tmpc[,.(value=sum(value),value.sd=Ssum(value.sd)),by=g_whoregion] #

tmp <- data.table(g_whoregion='Global',value=t1r9[,sum(value)],
                  value.sd=t1r9[,Ssum(value.sd)])
t1r9 <- rbind(t1r9,tmp)
t1r9[,quantity:='LYnewnotx2020']

t1r9[,.(value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]

save(t1r9,file=here('../figdat/t1r9.Rdata')) #life-years untreated



## === figure data

## fig 3
NU <- merge(estl[,.(age=agenow,year,sex,iso3,alive,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
NU <- NU[!is.na(acat)]
NU <- NU[,.(alive=sum(alive)),by=.(acat,sex,g_whoregion)]
NU[,type:='untreated']

save(NU,file=here('../figdat/NU.Rdata'))

## figure 2
untx <- estl[,.(LYS.t=sum(LYS)),by=.(g_whoregion,sex,acat)]
untx[,type:="untreated"]

save(untx,file=here('../figdat/untx.Rdata'))

## fig 2b
estl <- merge(estl,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
estl[,YPT:=2020-year]
untxx <- estl[,.(YPT=weighted.mean(YPT,w=alive)),by=.(g_whoregion,sex,acats)] #stats for article
untxx[,type:="untreated"]

save(untxx,file=here('../figdat/untxx.Rdata'))

## fig 2c
estl[,acats:=NULL]
estl <- merge(estl,lamap[,.(age,acats)],by='age') #acats now category for age of TB
## estl[,YPT:=2020-year]
untxx2 <- estl[,.(YPT=weighted.mean(YPT,w=alive)),by=.(g_whoregion,sex,acats)] #stats for article
untxx2[,type:="untreated"]

save(untxx2,file=here('../figdat/untxx2.Rdata'))

## for age and yl summaries
utay <- estl[,.(wt=alive,YPT,agenow)]
save(utay,file=here('../figdat/utay.Rdata'))


## === extra stats
utxknum <- estl[acat %in% c('0-4','5-14'),sum(alive)]
utxkden <- estl[,sum(alive)]
utxknumy <- estl[acat %in% c('0-4','5-14'),sum(LYS)]
utxkdeny <- estl[,sum(LYS)]

save(utxknum,file=here('../tmpdata/utxknum.Rdata'))
save(utxkden,file=here('../tmpdata/utxkden.Rdata'))
save(utxknumy,file=here('../tmpdata/utxknumy.Rdata'))
save(utxkdeny,file=here('../tmpdata/utxkdeny.Rdata'))

## for paediatric LYs
tmp <- estl[acat %in% c('0-4','5-14'),.(value=sum(LYS),value.sd=Ssum(LYS.sd)),
            by=.(iso3,g_whoregion,year)]
tmp <- tmp[,.(value=sum(value),value.sd=sum(value.sd)),by=.(iso3,g_whoregion)]
LYpu <- tmp[,.(value=sum(value),value.sd=Ssum(value.sd))]
save(LYpu,file=here('../tmpdata/LYpu.Rdata'))


## for male survivors
tmp <- estl[sex=='Male',.(value=sum(alive),value.sd=Ssum(alive.sd)),
            by=.(iso3,g_whoregion,year)]
tmp <- tmp[,.(value=sum(value),value.sd=sum(value.sd)),by=.(iso3,g_whoregion)]
SMu <- tmp[,.(value=sum(value),value.sd=Ssum(value.sd))]
save(SMu,file=here('../tmpdata/SMu.Rdata'))

## for HIV survivors
tmp <- estl[sex=='Male',.(value=sum(alive.h),value.sd=Ssum(alive.h.sd)),
            by=.(iso3,g_whoregion,year)]
tmp <- tmp[,.(value=sum(value),value.sd=sum(value.sd)),by=.(iso3,g_whoregion)]
SHu <- tmp[,.(value=sum(value),value.sd=Ssum(value.sd))]
save(SHu,file=here('../tmpdata/SHu.Rdata'))

## for HIV new cases
save(NHu,file=here('../tmpdata/NHu.Rdata'))

