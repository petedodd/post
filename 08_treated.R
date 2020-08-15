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
N3[,alive.sd:=xfun(value,S,value.sd,S.sd)]               #not used
N3[,LYS.sd:=xfun(LY,value,LY.sd,value.sd)]                #not used
N3[,alive.0.sd:=xfun3(value,S.0,(1-H),value.sd,S.0.sd,H.sd)]
N3[,LYS.0.sd:=xfun3(value,LY.0,(1-H),value.sd,LY.0.sd,H.sd)]
N3[,alive.h.sd:=xfun3(value,S.h,H,value.sd,S.h.sd,H.sd)]
N3[,LYS.h.sd:=xfun3(value,LY.h,H,value.sd,LY.h.sd,H.sd)]
## N3[,alive.sd:=S.sd*value]               #not used
## N3[,LYS.sd:=LY.sd*value]                #not used
## N3[,alive.0.sd:=value * xfun(S.0,(1-H),S.0.sd,H.sd)]
## N3[,LYS.0.sd:=value * xfun(LY.0,(1-H),LY.0.sd,H.sd)]
## N3[,alive.h.sd:=value * xfun(S.h,H,S.h.sd,H.sd)]
## N3[,LYS.h.sd:=value * xfun(LY.h,H,LY.h.sd,H.sd)]



summary(N3[,.(alive.sd,LYS.sd)])
N3[!is.finite(alive.sd)]                #odd EGY early 1980s in 99 yo
N3[!is.finite(alive.sd),c('hsa','haa','S.sd','LY.sd','S.0.sd','LY.0.sd','alive.sd','LYS.sd','alive.0.sd','LYS.0.sd'):=0] #safety
summary(N3[,.(alive.sd,LYS.sd)])


names(N3)

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

## ============  figdat =================
N3 <- merge(N3,isokey,by = 'iso3')


## === for table 1
t1r5 <- N3[,.(value=sum(alive.t),
              value.sd=Ssum(alive.t.sd)),
           by=g_whoregion]
tmp <- data.table(g_whoregion='Global',
                  value=t1r5[,sum(value)],
                  value.sd=t1r5[,Ssum(value.sd)])
t1r5 <- rbind(t1r5,tmp)
t1r5[,quantity:='totnewtx2020']

t1r5[,.(value.sd/value,value-value.sd,value+value.sd)]

save(t1r5,file=here('../figdat/t1r5.Rdata')) #treated survivors

t1r8 <- N3[,.(value=sum(LYS.t),value.sd=sum(LYS.t.sd)),by=g_whoregion]
tmp <- data.table(g_whoregion='Global',
                  value=t1r8[,sum(value)],
                  value.sd=t1r8[,sum(value.sd)])

t1r8 <- rbind(t1r8,tmp)
t1r8[,quantity:='LYnewtx2020']

t1r8[,.(value.sd/value,value-value.sd,value+value.sd)]

save(t1r8,file=here('../figdat/t1r8.Rdata')) #lifeyears among treated

## === for table 2

if(! 'agenow' %in% names(N3)) N3[,agenow:=age + 2020-year]

## EMR paed graph
tmp1 <- N3[year>=2015 & g_whoregion=='EMR',.(total=sum(alive)),by=.(iso3)]
tmp2 <- N3[year>=2015 & g_whoregion=='EMR' & agenow<15,.(kids=sum(alive)),by=.(iso3)]
tmp1 <- merge(tmp1,tmp2,by='iso3')
tmp1[,pcp:=1e2*kids/total]


GP <- ggpubr::ggdotchart(tmp1, x = "iso3", y = "total",
                         sorting = "descending",
                         add = "segments",
                         add.params = list(color = "black", size = 1),
                         rotate = TRUE,
                         dot.size = 6,
                         label = round(tmp1$pcp),
                         font.label = list(color = "white", size = 9,vjust = 0.5)
                         ) + ggpubr::grids() +
              scale_y_continuous(label=absspace) +
              ylab('Treated TB survivors alive 2020') +
              xlab('Country ISO3 code')

ggsave(GP,file=here('../plots/EMRpaed.pdf'),w=6,h=6)


## within 5 years
c1 <- rbind(N3[year>=2015,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
            data.table(g_whoregion='Global',
                       N3[year>=2015,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )

c1h <- N3[year>=2015 & g_whoregion=='AFR',
          .(total=sum(alive.h),total.sd=Ssum(alive.h.sd))]
c1m <- N3[year>=2015 & sex=='Male',
          .(total=sum(alive.t),total.sd=Ssum(alive.t.sd))]


tmp <- N3[year>=2015,.(total=sum(alive)),by=.(g_whoregion,sex)]
tmp[,tot:=sum(total),by=g_whoregion]
tmp[,pc:=1e2*total/tot]
tmp <- tmp[sex=='Male']
c2 <- rbind(tmp[,.(g_whoregion,pc)],
            data.table(g_whoregion='Global',pc=1e2*N3[year>=2015 & sex=='Male',
                                               sum(alive)]/
                                              N3[year>=2015,sum(alive)])
            )

tmp <- N3[year>=2015 & agenow<15,.(total=sum(alive)),by=.(g_whoregion)]
tmp2 <- N3[year>=2015,.(total=sum(alive)),by=.(g_whoregion)]
tmp <- merge(tmp,tmp2,by='g_whoregion')
tmp[,pck:=1e2*total.x/total.y]
c4 <- rbind(tmp[,.(g_whoregion,pck)],
            data.table(g_whoregion='Global',pck=1e2*N3[year>=2015 & agenow<15,
                                               sum(alive)]/
                                              N3[year>=2015,sum(alive)])
            )

## less than 2 years
c1b <- rbind(N3[year>=2018,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
            data.table(g_whoregion='Global',
                       N3[year>=2018,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )

c1hb <- N3[year>=2018 & g_whoregion=='AFR',
           .(total=sum(alive.h),total.sd=Ssum(alive.h.sd))]
c1mb <- N3[year>=2018 & sex=='Male',
          .(total=sum(alive.t),total.sd=Ssum(alive.t.sd))]

tmp <- N3[year>=2018,.(total=sum(alive)),by=.(g_whoregion,sex)]
tmp[,tot:=sum(total),by=g_whoregion]
tmp[,pc:=1e2*total/tot]
tmp <- tmp[sex=='Male']
c2b <- rbind(tmp[,.(g_whoregion,pc)],
            data.table(g_whoregion='Global',pc=1e2*N3[year>=2018 & sex=='Male',
                                               sum(alive)]/
                                              N3[year>=2018,sum(alive)])
            )



tmp <- N3[year>=2018 & agenow<15,.(total=sum(alive)),by=.(g_whoregion)]
tmp2 <- N3[year>=2018,.(total=sum(alive)),by=.(g_whoregion)]
tmp <- merge(tmp,tmp2,by='g_whoregion')
tmp[,pck:=1e2*total.x/total.y]
c4b <- rbind(tmp[,.(g_whoregion,pck)],
            data.table(g_whoregion='Global',pck=1e2*N3[year>=2018 & agenow<15,
                                               sum(alive)]/
                                              N3[year>=2018,sum(alive)])
            )

## age data
cage <- N3[year>=2015,.(agenow,alive)]
cageb <- N3[year>=2018,.(agenow,alive)]

## save data
save(c1,file=here("../figdat/c1.Rdata")); save(c1b,file=here("../figdat/c1b.Rdata"))
save(c2,file=here("../figdat/c2.Rdata")); save(c2b,file=here("../figdat/c2b.Rdata"))
save(c4,file=here("../figdat/c4.Rdata")); save(c4b,file=here("../figdat/c4b.Rdata"))
## HIV
save(c1h,file=here("../figdat/c1h.Rdata")); save(c1hb,file=here("../figdat/c1hb.Rdata"))
## M/F
save(c1m,file=here("../figdat/c1m.Rdata")); save(c1mb,file=here("../figdat/c1mb.Rdata"))
## age data
save(cage,file=here("../figdat/cage.Rdata")); save(cageb,file=here("../figdat/cageb.Rdata"))

## === for figures

## fig 3
NZ <- merge(N3[,.(age=agenow,year,sex,iso3,alive.t,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
NZ <- NZ[!is.na(acat)]
NZ
NZ$acat <- factor(NZ$acat,levels=racts,ordered=TRUE)
NZ[year>=2015,type:='treated within 5 years']
NZ[year<2015,type:='treated over 5 years ago']
NZ <- NZ[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion,type)]

save(NZ,file=here('../figdat/NZ.Rdata'))


## fig 2b
## ie Years Post TB
N3[,YPT:=2020-year]
N3 <- merge(N3,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
## N3[,acats:=NULL]
N3RYLx <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #TODO stats for article
N3RYLx
N3RYLx[,type:="treated"]

save(N3RYLx,file=here('../figdat/N3RYLx.Rdata'))

## fig 2c
N3[,acats:=NULL]
N3 <- merge(N3,lamap[,.(age=age,acats)],by='age') #acats now category for age of TB
N3RYLx2 <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #TODO stats for article
N3RYLx2
N3RYLx2[,type:="treated"]

save(N3RYLx2,file=here('../figdat/N3RYLx2.Rdata'))

## for age and yl summaries
txay <- N3[,.(wt=alive.t,YPT,agenow)]
save(txay,file=here('../figdat/txay.Rdata'))


## regional summary plots
N3R <- N3[,.(alive.t=sum(alive.t),LYS.t=sum(LYS.t),
             alive.h=sum(alive.h),LYS.h=sum(LYS.h)),
          by=.(g_whoregion,year,sex,acat)]

N3R$acat <- factor(N3R$acat,levels=racts,ordered=TRUE)
yy <- sort(unique(N3R$year))
N3R$year <- factor(N3R$year,levels=yy,ordered=TRUE)

save(N3R,file=here('../figdat/N3R.Rdata'))

N3RYL <- N3R[,.(LYS.t=sum(LYS.t)),by=.(g_whoregion,sex,acat)]
N3RYL[,type:="treated"]

save(N3RYL,file=here('../figdat/N3RYL.Rdata'))

## extra stats
txknum <- N3[acat %in% c('0-4','5-14'),sum(alive.t)]
txkden <- N3[,sum(alive.t)]
txknumy <- N3[acat %in% c('0-4','5-14'),sum(LY)]
txkdeny <- N3[,sum(LY)]

save(txknum,file=here('../tmpdata/txknum.Rdata'))
save(txkden,file=here('../tmpdata/txkden.Rdata'))
save(txknumy,file=here('../tmpdata/txknumy.Rdata'))
save(txkdeny,file=here('../tmpdata/txkdeny.Rdata'))

## for paediatric LYs
LYpt <- N3[acat %in% c('0-4','5-14'),.(value=sum(LYS.t),value.sd=sum(LYS.t.sd))]
save(LYpt,file=here('../tmpdata/LYpt.Rdata'))

## for male survivors
SMt <- N3[sex=='Male',.(value=sum(alive.t),value.sd=Ssum(alive.t.sd))] #
save(SMt,file=here('../tmpdata/SMt.Rdata'))

## for HIV survivors
SHt <- N3[,.(value=sum(alive.h),value.sd=Ssum(alive.h.sd))] #
save(SHt,file=here('../tmpdata/SHt.Rdata'))
