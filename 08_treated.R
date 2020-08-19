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
N3 <- merge(N3,isokey,by = 'iso3')

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


## checking unc
N3[,summary(1e2*value.sd/value)]
N3[,summary(1e2*pd.sd/pd)]
N3[,summary(1e2*rat.sd/rat)]


## value.sd
tmp <- N3[,.(value=sum(value),
             value.sd=Ssum(value.sd)),
          by=iso3]

tmp[,summary(1e2*value.sd/value)]

## S.sd
tmp <- N3[,.(value=sum(S),
             value.sd=Ssum(S.sd)),
          by=iso3]
tmp[,summary(1e2*value.sd/value)]



## uncertainty
N3[,alive.sd:=xfun(value,S,value.sd,S.sd)]               #not used
N3[,LYS.sd:=xfun(LY,value,LY.sd,value.sd)]                #not used

## HIV stratified unc
N3[,alive.0.sd:=xfun3(value,S.0,(1-H),value.sd,S.0.sd,H.sd)]
N3[,LYS.0.sd:=xfun3(value,LY.0,(1-H),value.sd,LY.0.sd,H.sd)]
N3[,alive.h.sd:=xfun3(value,S.h,H,value.sd,S.h.sd,H.sd)]
N3[,LYS.h.sd:=xfun3(value,LY.h,H,value.sd,LY.h.sd,H.sd)]




summary(N3[,.(alive.sd,LYS.sd)])
N3[!is.finite(alive.sd)]                #odd EGY early 1980s in 99 yo
N3[!is.finite(alive.sd),c('hsa','haa','S.sd','LY.sd','S.0.sd','LY.0.sd','alive.sd','LYS.sd','alive.0.sd','LYS.0.sd'):=0] #safety
summary(N3[,.(alive.sd,LYS.sd)])


names(N3)

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



##aggregate over ages/sex assuming value uncor, S, H corr in each country
u1h <- N3[,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2h <- N3[,.(LYS.h=sum(LYS.h),
             LYS.h.sd=Sssumxyz(value,LY.h,H,value.sd,LY.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]

u1 <- N3[,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2 <- N3[,.(LYS.0=sum(LYS.0),
             LYS.0.sd=Sssumxyz(value,LY.0,1-H,value.sd,LY.0.sd,H.sd)),
          by=.(iso3,g_whoregion)]

u1 <- merge(u1,u1h,by=c('iso3','g_whoregion')) #alive
u2 <- merge(u2,u2h,by=c('iso3','g_whoregion')) #LYs

u1[,alive:=alive.h+alive.0]
u2[,LYS:=LYS.h+LYS.0]
u1[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)]
u2[,LYS.sd:=sqrt(LYS.h.sd^2+LYS.0.sd^2)]

summary(u2)
u1[,see(sum(alive))]
u2[,see(sum(LYS))]

## ============  figdat =================


## === for table 1

## --- tx alive

t1r5c <- u1[,.(iso3,value=alive,value.sd=alive.sd)]
t1r5c[,quantity:='totnewtx2020']
save(t1r5c,file=here('../figdat/t1r5c.Rdata'))

t1r5c[,summary(1e2*value.sd/value)]
t1r5c[value.sd>value,.(iso3,value.sd/value)]

## regional level
t1r5 <- u1[,.(value=sum(alive),
              value.sd=Ssum(alive.sd)),
             by=g_whoregion]

tmp <- data.table(g_whoregion='Global',
                  value=t1r5[,sum(value)],
                  value.sd=t1r5[,Ssum(value.sd)])
t1r5 <- rbind(t1r5,tmp)
t1r5[,quantity:='totnewtx2020']

t1r5[,.(1e2*value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]

save(t1r5,file=here('../figdat/t1r5.Rdata')) #treated survivors


## --- LYS

## ## country level

t1r8c <- u2[,.(iso3,value=LYS,value.sd=LYS.sd)]

t1r8c[,quantity:='LYnewtx2020']
save(t1r8c,file=here('../figdat/t1r8c.Rdata'))

t1r8c[,summary(1e2*value.sd/value)]

## regional level
t1r8 <- u2[,.(value=sum(LYS),
              value.sd=Ssum(LYS.sd)),
             by=g_whoregion]

tmp <- data.table(g_whoregion='Global',
                  value=t1r8[,sum(value)],
                  value.sd=t1r8[,Ssum(value.sd)])
t1r8 <- rbind(t1r8,tmp)
t1r8[,quantity:='LYnewtx2020']

t1r8[,.(1e2*value.sd/value,see(value),see(value-2*value.sd),see(value+2*value.sd))]

save(t1r8,file=here('../figdat/t1r8.Rdata')) #lifeyears among treated


## for HIV survivors
SHt <- u1[,.(value=sum(alive.h),value.sd=Ssum(alive.h.sd))] #
save(SHt,file=here('../tmpdata/SHt.Rdata'))


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

cat(1e2*tmp1[iso3=='PAK',total]/tmp1[,sum(total)],file=here('texto/PAKinEMR.txt'))


## within 5 years
u1h <- N3[year>=2015,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2h <- N3[year>=2018,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive

u1 <- N3[year>=2015,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2 <- N3[year>=2018,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive

u1 <- merge(u1,u1h,by=c('iso3','g_whoregion')) #alive 5 year
u2 <- merge(u2,u2h,by=c('iso3','g_whoregion')) #alive 2 year

u1[,alive:=alive.h+alive.0];
u2[,alive:=alive.h+alive.0]
u1[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)];
u2[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)]

## HIV in AFR
c1h <- u1[g_whoregion=='AFR',.(total=sum(alive.h),total.sd=Ssum(alive.h.sd))]
c1hb <- u2[g_whoregion=='AFR',.(total=sum(alive.h),total.sd=Ssum(alive.h.sd))]

## combine
c1 <- rbind(u1[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
            data.table(g_whoregion='Global',
                       u1[,.(total=sum(alive),total.sd=Ssum(alive.sd))] )  )
c1b <- rbind(u2[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
            data.table(g_whoregion='Global',
                       u2[,.(total=sum(alive),total.sd=Ssum(alive.sd))] )  )

c1[,.(see(total),see(total-2*total.sd),see(total+2*total.sd))]



## for pc male in table 2
u1mh <- N3[sex=='Male' & year>=2015,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #alive
u1m <- N3[sex=='Male' & year>=2015,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u1m <- merge(u1m,u1mh,by=c('iso3','g_whoregion')) #alive 5 year
u2mh <- N3[sex=='Male' & year>=2018,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #alive
u2m <- N3[sex=='Male' & year>=2018,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2m <- merge(u2m,u2mh,by=c('iso3','g_whoregion')) #alive 5 year


u1m[,alive:=alive.h+alive.0];
u2m[,alive:=alive.h+alive.0]
u1m[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)];
u2m[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)]


## for male survivors
u1mha <- N3[sex=='Male',.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #alive
u1ma <- N3[sex=='Male',.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
SMt <- merge(u1ma,u1mha,by=c('iso3','g_whoregion'))
SMt[,alive:=alive.h+alive.0]
SMt[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)];
SMt <- SMt[,.(value=sum(alive),value.sd=Ssum(alive.sd))]

save(SMt,file=here('../tmpdata/SMt.Rdata'))


c1m <- rbind(u1m[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
             data.table(g_whoregion='Global',
                        u1m[,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )
c1mb <- rbind(u2m[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
             data.table(g_whoregion='Global',
                        u2m[,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )

c2 <- merge(c1m,c1,by='g_whoregion')
c2[,pc:=1e2*total.x/total.y]
c2 <- c2[,.(g_whoregion,pc)]

c2b <- merge(c1mb,c1b,by='g_whoregion')
c2b[,pc:=1e2*total.x/total.y]
c2b <- c2b[,.(g_whoregion,pc)]

## for pc kids in table 2
u1kh <- N3[agenow<15 & year>=2015,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #alive
u1k <- N3[agenow<15 & year>=2015,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u1k <- merge(u1k,u1kh,by=c('iso3','g_whoregion')) #alive 5 year
u2kh <- N3[agenow<15 & year>=2018,.(alive.h=sum(alive.h),
            alive.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #alive
u2k <- N3[agenow<15 & year>=2018,.(alive.0=sum(alive.0),
            alive.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
         by=.(iso3,g_whoregion)]   #alive
u2k <- merge(u2k,u2kh,by=c('iso3','g_whoregion')) #alive 5 year


u1k[,alive:=alive.h+alive.0];
u2k[,alive:=alive.h+alive.0]
u1k[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)];
u2k[,alive.sd:=sqrt(alive.h.sd^2+alive.0.sd^2)]



c4 <- rbind(u1k[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
             data.table(g_whoregion='Global',
                        u1k[,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )
c4b <- rbind(u2k[,.(total=sum(alive),total.sd=Ssum(alive.sd)),by=g_whoregion],
             data.table(g_whoregion='Global',
                        u2k[,.(total=sum(alive),total.sd=Ssum(alive.sd))]) )

c4 <- merge(c4,c1,by='g_whoregion')
c4[,pck:=1e2*total.x/total.y]
c4 <- c4[,.(g_whoregion,pck)]

c4b <- merge(c4b,c1b,by='g_whoregion')
c4b[,pck:=1e2*total.x/total.y]
c4b <- c4b[,.(g_whoregion,pck)]


## for paediatric LYs
u1kha <-  N3[age<15,.(LYS.h=sum(LYS.h),
            LYS.h.sd=Sssumxyz(value,S.h,H,value.sd,S.h.sd,H.sd)),
           by=.(iso3,g_whoregion)]   #LYS
u1ka <- N3[age<15,.(LYS.0=sum(LYS.0),
            LYS.0.sd=Sssumxyz(value,S.0,1-H,value.sd,S.0.sd,H.sd)),
          by=.(iso3,g_whoregion)]   #LYS
LYpt <- merge(u1ka,u1kha,by=c('iso3','g_whoregion'))
LYpt[,LYS:=LYS.h+LYS.0]
LYpt[,LYS.sd:=sqrt(LYS.h.sd^2+LYS.0.sd^2)];
LYpt <- LYpt[,.(value=sum(LYS),value.sd=Ssum(LYS.sd))]

save(LYpt,file=here('../tmpdata/LYpt.Rdata'))



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
N3[,acats:=NULL]
N3 <- merge(N3,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
N3RYLx <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #stats for article
N3RYLx
N3RYLx[,type:="treated"]

save(N3RYLx,file=here('../figdat/N3RYLx.Rdata'))

## fig 2c
N3[,acats:=NULL]
N3 <- merge(N3,lamap[,.(age=age,acats)],by='age') #acats now category for age of TB

N3RYLx2 <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #stats for article
N3RYLx2
N3RYLx2[,type:="treated"]

save(N3RYLx2,file=here('../figdat/N3RYLx2.Rdata'))

## for age and yl summaries
txay <- N3[,.(wt=alive.t,YPT,agenow)]
save(txay,file=here('../figdat/txay.Rdata'))


## regional summary plots
## N3[,acat:=NULL]
## N3 <- merge(N3,lamap[,.(age=age,acat=acats)],by='age') #acats as age TB
N3[,acats:=NULL]
N3R <- N3[,.(alive.t=sum(alive.t),LYS.t=sum(LYS.t), #acat not found
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
