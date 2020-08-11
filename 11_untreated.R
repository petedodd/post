## untreated survivors
rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../figdat/IX.Rdata')
  if(file.exists(fn))stop('Not running as figdat/IX.Rdata exists!')
}

## load data
load(here('../tmpdata/estl.Rdata'))
load(here('../tmpdata/estg.Rdata'))
load(here('../tmpdata/LA.Rdata'))


## merge against life tables etc
estl <- merge(estl,LA,by=c('iso3','year','age','sex'),all.x = TRUE,all.y = FALSE)
estl[,range(year)]
estl

## calcualtions
estl[,alive.0:=S.0*gapls.0]
estl[,alive.h:=S.h*gapls.h]
estl[,alive:=alive.0+alive.h]
estl[,LYS.0:=LY.0*gapls.0]
estl[,LYS.h:=LY.h*gapls.h]
estl[,LYS:=LYS.h + LYS.0]

## uncertainty
estl[,alive.0.sd:=xfun(S.0,gapls.0,S.0.sd,gapls.0.sd)]
estl[,alive.h.sd:=xfun(S.h,gapls.h,S.h.sd,gapls.h.sd)]
estl[,alive.sd:=sqrt(alive.0.sd^2+alive.h.sd^2)]
estl[,LYS.0.sd:=xfun(LY.0,gapls.0,LY.0.sd,gapls.0.sd)]
estl[,LYS.h.sd:=xfun(LY.h,gapls.h,LY.h.sd,gapls.h.sd)]
estl[,LYS.sd:=sqrt(LYS.h.sd^2 + LYS.0.sd^2)]

summary(estl[,.(alive.sd,LYS.sd)])
estl[!is.finite(alive.sd)]              #EGY 1980-1982 in 99 yo
estl[!is.finite(alive.sd),c('S.sd','LY.sd','S.0.sd','LY.0.sd','alive.sd','LYS.sd','alive.0.sd','LYS.0.sd'):=0] #safety
## summary(N3[,.(alive.sd,LYS.sd)])

## ## totals
## N3[!iso3 %in% hivcountries,c('alive.t',
##                              'alive.t.sd',
##                              'LYS.t','LYS.t.sd'):=
##                              list(alive.0,
##                                   alive.0.sd,
##                                   LYS.0,
##                                   LYS.0.sd)]

## N3[iso3 %in% hivcountries,c('alive.t',
##                              'alive.t.sd',
##                              'LYS.t','LYS.t.sd'):=
##                              list(alive.0 + alive.h,
##                                   sqrt(alive.0.sd^2 + alive.h.sd^2),
##                                   LYS.0 + LYS.h,
##                                   sqrt(LYS.0.sd^2 + LYS.h.sd^2))]

## N3[!is.finite(LYS.t.sd),c('S.sd','LYS.sd','alive.sd',
##                           'LYS.t.sd','alive.t.sd',
##                           'S.0.sd','LYS.0.sd','alive.0.sd',
##                           'S.h.sd','LYS.h.sd','alive.h.sd'):=0]

## N3[!is.finite(alive.t.sd)]

estl[,sum(gapl,na.rm=TRUE)]/1e6                     #197m
estl[,sum(gapls,na.rm=TRUE)]/1e6                     #98m
estl[,sum(alive,na.rm=TRUE)]/1e6                     #52.4m
estl[,sum(LYS,na.rm=TRUE)]/1e9                     #1.2 bn


## Too small (uncorrelated)
estl[,Ssum(gapl.sd,na.rm=TRUE)]/1e6                     #0.76m
estl[,Ssum(gapls,na.rm=TRUE)]/1e6                     #0.66m
estl[,Ssum(alive,na.rm=TRUE)]/1e6                     #0.46m
estl[,Ssum(LYS,na.rm=TRUE)]/1e9                     #0.01bn

## Too large (perfectly correlated)
estl[,sum(gapl.sd,na.rm=TRUE)]/1e6                     #113m
estl[,sum(gapls,na.rm=TRUE)]/1e6                     #117m
estl[,sum(alive,na.rm=TRUE)]/1e6                     #66m
estl[,sum(LYS,na.rm=TRUE)]/1e9                     #1.7bn

## jjk 
## TODO think how to have correlated survival but uncorrelated CDR
## TODO CDR in sub-groups 

40 *                                    #years
  5 *                                   #million gap per year
  0.5 *                                 #CFR
  0.6                                   #mean survival?
## 60 miln

## --- plots
## age now! NOTE
estl[,agenow:=age + 2020-year]

estl

## figdat
## estg
t1r1 <- est[,.(value=sum(e_inc_num),
               value.sd=Ssum((ocdr.sd/ocdr)*e_inc_num)),
            by=g_whoregion] #NOTE this is new
tmp <- data.table(g_whoregion='Global',
                  value=t1r1[,sum(value)],
                  value.sd=t1r1[,Ssum(value.sd)])
t1r1 <- rbind(t1r1,tmp)
t1r1[,quantity:='totnew']

t1r1[,.(value.sd/value,value-value.sd,value+value.sd)]



save(t1r1,file=here('../figdat/t1r1.Rdata')) #new total

## check
estl[g_whoregion=='EMR',.(value=sum(gap),value.sd=Ssum(gap.sd),
                         Ssum(gap.sd)/sum(gap)),by=iso3] #

t1r3 <- estl[,.(value=sum(gap),value.sd=Ssum(gap.sd)),by=g_whoregion] #
tmp <- data.table(g_whoregion='Global',
                  value=t1r3[,sum(value)],
                  value.sd=t1r3[,Ssum(value.sd)])
t1r3 <- rbind(t1r3,tmp)
t1r3[,quantity:='totnotx']

t1r3[,.(value.sd/value,value-value.sd,value+value.sd)]

save(t1r3,file=here('../figdat/t1r3.Rdata')) #new untreated

if(cr){
  t1r6 <- estl[,.(value=sum(alive),value.sd=sum(alive.sd)),by=g_whoregion] #
  tmp <- data.table(g_whoregion='Global',value=t1r6[,sum(value)],
                    value.sd=t1r6[,sum(value.sd)])
} else{
    t1r6 <- estl[,.(value=sum(alive),value.sd=Ssum(alive.sd)),by=g_whoregion] #
    tmp <- data.table(g_whoregion='Global',value=t1r6[,Ssum(value)],
                    value.sd=t1r6[,sum(value.sd)])
}
t1r6 <- rbind(t1r6,tmp)
t1r6[,quantity:='totnotx2020']

t1r6[,.(value.sd/value,value-value.sd,value+value.sd)]

save(t1r6,file=here('../figdat/t1r6.Rdata')) #untreated survivors

if(cr){
  t1r9 <- estl[,.(value=sum(LYS),value.sd=sum(LYS.sd)),by=g_whoregion] #
  tmp <- data.table(g_whoregion='Global',value=t1r9[,sum(value)],
                    value.sd=t1r9[,sum(value.sd)])
} else {
  t1r9 <- estl[,.(value=sum(LYS),value.sd=Ssum(LYS.sd)),by=g_whoregion] #
  tmp <- data.table(g_whoregion='Global',value=t1r9[,sum(value)],
                    value.sd=t1r9[,Ssum(value.sd)])
}
t1r9 <- rbind(t1r9,tmp)
t1r9[,quantity:='LYnewnotx2020']

t1r9[,.(value.sd/value,value-value.sd,value+value.sd)]

save(t1r9,file=here('../figdat/t1r9.Rdata')) #life-years untreated


IX <- merge(estl[,.(age=agenow,sex,iso3,alive,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
IX <- IX[!is.na(acat)]
IX

IX$acat <- factor(IX$acat,levels=racts,ordered=TRUE)

tmp <- IX[,.(alive=sum(alive)),by=.(acat,sex,g_whoregion)]

save(IX,file=here('../figdat/IX.Rdata'))


