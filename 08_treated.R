## looking at treated survival

rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../figdat/t1r8.Rdata')
  if(file.exists(fn))stop('Not running as figdat/t1r8.Rdata exists!')
}

## load data
load(here('../tmpdata/N3.Rdata'))
load(here('../tmpdata/isokey.Rdata'))
load(here('../tmpdata/LA.Rdata'))


## merge with new LA here
N3 <- merge(N3,LA,by=c('iso3','year','age','sex'),all.x = TRUE,all.y = FALSE)
N3[,range(year)]

## previous version
N3[,alive:=S*value]
N3[,LYS:=LY*value]
N3[,sum(value)]/1e6                     #162

## HIV stratified
N3[,alive.0:=S.0*value*(1-H)]
N3[,LYS.0:=LY.0*value*(1-H)]
N3[,alive.h:=S.h*value*(H)]
N3[,LYS.h:=LY.h*value*(H)]

## uncertainty TODO unc check correlations UNC jak
N3[,alive.sd:=S.sd*value]               #not used
N3[,LYS.sd:=LY.sd*value]                #not used
N3[,alive.0.sd:=value * xfun(S.0,(1-H),S.0.sd,H.sd)]
N3[,LYS.0.sd:=value * xfun(LY.0,(1-H),LY.0.sd,H.sd)]
N3[,alive.h.sd:=value * xfun(S.h,H,S.h.sd,H.sd)]
N3[,LYS.h.sd:=value * xfun(LY.h,H,LY.h.sd,H.sd)]


summary(N3[,.(alive.sd,LYS.sd)])
N3[!is.finite(alive.sd)]                #odd EGY early 1980s in 99 yo
N3[!is.finite(alive.sd),c('hsa','haa','S.sd','LY.sd','S.0.sd','LY.0.sd','alive.sd','LYS.sd','alive.0.sd','LYS.0.sd'):=0] #safety
summary(N3[,.(alive.sd,LYS.sd)])

## N3[,alive.0.sd:=value*sqrt(S.0.sd^2*(1-H)^2 + S.0.sd^2*(H.sd^2) + S.0^2*(H.sd^2))]
## N3[,LYS.0.sd:=value*sqrt(LY.0.sd^2*(1-H)^2 + LY.0.sd^2*(H.sd^2) + LY.0^2*(H.sd^2))]
## N3[,alive.h.sd:=value*sqrt(S.h^2*(H.sd)^2 + S.h.sd^2*(H.sd)^2 + S.h.sd^2*(H)^2)]
## N3[,LYS.h.sd:=value*sqrt(LY.h.sd^2*(H^2) + LY.h.sd^2*(H.sd^2) + LY.h^2*(H.sd^2))]

## totals
N3[!iso3 %in% hivcountries,c('alive.t',
                             'alive.t.sd',
                             'LYS.t','LYS.t.sd'):=
                             list(alive.0,
                                  alive.0.sd,
                                  LYS.0,
                                  LYS.0.sd)]

N3[iso3 %in% hivcountries,c('alive.t',
                             'alive.t.sd',
                             'LYS.t','LYS.t.sd'):=
                             list(alive.0 + alive.h,
                                  sqrt(alive.0.sd^2 + alive.h.sd^2),
                                  LYS.0 + LYS.h,
                                  sqrt(LYS.0.sd^2 + LYS.h.sd^2))]

N3[!is.finite(LYS.t.sd),c('S.sd','LYS.sd','alive.sd',
                          'LYS.t.sd','alive.t.sd',
                          'S.0.sd','LYS.0.sd','alive.0.sd',
                          'S.h.sd','LYS.h.sd','alive.h.sd'):=0]

N3[!is.finite(alive.t.sd)]


N3 <- merge(N3,lamap[,.(age,acat=acats)],by='age',all.x=TRUE)

## figdat
N3 <- merge(N3,isokey,by = 'iso3')

if(cr){
  t1r5 <- N3[,.(value=sum(alive.t),
                value.sd=sum(alive.t.sd)),
             by=g_whoregion] #NOTE assuming perfect correlation?
  tmp <- data.table(g_whoregion='Global',
                    value=t1r5[,sum(value)],
                    value.sd=t1r5[,sum(value.sd)])
} else {
  t1r5 <- N3[,.(value=sum(alive.t),
                value.sd=Ssum(alive.t.sd)),
             by=g_whoregion] #NOTE assuming perfect correlation?
  tmp <- data.table(g_whoregion='Global',
                    value=t1r5[,sum(value)],
                    value.sd=t1r5[,Ssum(value.sd)])
}
t1r5 <- rbind(t1r5,tmp)
t1r5[,quantity:='totnewtx2020']

save(t1r5,file=here('../figdat/t1r5.Rdata'))

if(cr){
  t1r8 <- N3[,.(value=sum(LYS.t),value.sd=sum(LYS.t.sd)),by=g_whoregion]
  tmp <- data.table(g_whoregion='Global',
                    value=t1r8[,sum(value)],
                    value.sd=t1r8[,sum(value.sd)])
} else{
  t1r8 <- N3[,.(value=sum(LYS.t),value.sd=sum(LYS.t.sd)),by=g_whoregion]
  tmp <- data.table(g_whoregion='Global',
                    value=t1r8[,sum(value)],
                    value.sd=t1r8[,sum(value.sd)])
}
t1r8 <- rbind(t1r8,tmp)
t1r8[,quantity:='LYnewtx2020']

save(t1r8,file=here('../figdat/t1r8.Rdata'))

