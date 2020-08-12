## this file is for processing and presenting results

rm(list=ls())

library(ggplot2)
library(scales)
library(data.table)
library(ggthemes)
library(here)

absspace <- function(x,...) {             #works
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

whozg <- c('AFR','AMR','EMR','EUR','SEA','WPR')
whozt <- c('Africa','The Americas','Eastern Mediterranean','Europe','South-East Asia',
           'Western Pacific')
wrk <- data.table(g_whoregion=whozg,name=whozt)

## ====== demography turn N into 1 year age group interpolator for right years
acts <- c('04','514','1524','2534','3544','4554','5564','65')
racts <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')
rracts <- racts[c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8,8,8,8,8,8)]
rracts2 <- racts[c(rep(1,5),rep(2:7,each=10),rep(8,length(66:100)))]

## load data
load(here('../tmpdata/N3.Rdata'))
est <- fread(here('../indata/TB_burden_countries_2020-02-24.csv'))
est <- est[,.(iso3,year,e_inc_num)]
load(here('../indata/N_simple.Rdata'))           #5 year age groups

## ------- additional stats to move
## proportion covered
cnsdone <- N3[,(unique(iso3))]
cat(cnsdone,file=here('texto/cnsdone.txt'))


estr <- est[year==2018]
(pcnttbdone <- 1e2*estr[iso3 %in% cnsdone,sum(e_inc_num)]/estr[,sum(e_inc_num)])
cat(pcnttbdone,file=here('texto/pcnttbdone.txt'))


N[Year==2020,sum(PopTotal)]/1e6

(pcntpopdone <- 1e2*N[Year==2020 & iso3 %in% cnsdone,sum(PopTotal)]/N[Year==2020,sum(PopTotal)])
cat(pcntpopdone,file=here('texto/pcntpopdone.txt'))


## TB props in kids
load(file=here('../tmpdata/txknum.Rdata'))
load(file=here('../tmpdata/txkden.Rdata'))
load(file=here('../tmpdata/txknumy.Rdata'))
load(file=here('../tmpdata/txkdeny.Rdata'))

1e2*txknum/txkden #1e2*N3[acat %in% c('0-4','5-14'),sum(alive.t)]/N3[,sum(alive.t)]
1e2*estl[acat %in% c('0-4','5-14'),sum(alive)]/estl[,sum(alive)]

1e2*(estl[acat %in% c('0-4','5-14'),sum(alive)] + N3[acat %in% c('0-4','5-14'),sum(alive.t)]) / (N3[,sum(alive.t)]+estl[,sum(alive)])
## 10%

1e2*(estl[acat %in% c('0-4','5-14'),sum(LY)] + N3[acat %in% c('0-4','5-14'),sum(LY)]) / (N3[,sum(LY)]+estl[,sum(LY)])
## 26%

length(cnsdone)

## 1
## ## ## --- Results ---
## ## ## reporting
## N3S <- N3[,.(alive=sum(alive),LYS=sum(LYS)),by=.(iso3,year,sex,acat)]
## N3S[,sum(alive)]/1e6                  #98m
## N3S[,sum(LYS)]/1e9                    #1.97 bn

## N3S <- N3[,.(alive=sum(alive.t),LYS=sum(LYS.t)),by=.(iso3,year,sex,acat)]
## N3S[,sum(alive)]/1e6                  #98m
## N3S[,sum(LYS)]/1e9                    #2.0 bn

## N3S.h <- N3[,.(alive=sum(alive.h),LYS=sum(LYS.h)),by=.(iso3,year,sex,acat)]
## N3S.h[,sum(alive)]/1e6                    #3.3 mln
## N3S.h[,sum(LYS)]/1e6                    #54 mln

## ## regional aggregation
## N3S$acat <- factor(N3S$acat,levels=racts,ordered=TRUE)
## yy <- N3S[,unique(year)]
## yy <- sort(yy)
## N3S$year <- factor(N3S$year,levels=yy,ordered=TRUE)

## names(N3)
## N3[,unique(age)]
## racts

## N3[,acat:=rracts2[1+age]]

## N3$acat <- factor(N3$acat,levels=racts,ordered=TRUE)

## ## regional summary
## N3R <- N3[,.(alive.t=sum(alive.t),LYS.t=sum(LYS.t),
##              alive.h=sum(alive.h),LYS.h=sum(LYS.h)),
##           by=.(g_whoregion,year,sex,acat)]

## N3R$acat <- factor(N3R$acat,levels=racts,ordered=TRUE)
## N3R$year <- factor(N3R$year,levels=yy,ordered=TRUE)

## N3S <- merge(N3S,unique(N3[,.(iso3,g_whoregion)]),by='iso3')

## tmp <- N3S[,.(alive=sum(alive)),by=.(acat,sex,year,g_whoregion)]
## tmp[year==2010]

## N3R[g_whoregion=='AMR' & acat=='0-4' & sex=='Male',sort(year)]


## ## LY version
## N3RY <- N3R[,.(LYS.t=sum(LYS.t)),by=.(g_whoregion,sex,acat)]

load(here('estl.Rdata'))


## N3R[,last5:=ifelse(year>=2015,'yes','no')]
## N3RYL <- N3R[,.(LYS.t=sum(LYS.t)),by=.(last5,g_whoregion,sex,acat)]
load(file=here('../figdat/N3RYL.Rdata'))

## N3RYL <- N3R[,.(LYS.t=sum(LYS.t)),by=.(g_whoregion,sex,acat)]
untx <- estl[,.(LYS.t=sum(LYS)),by=.(g_whoregion,sex,acat)]
untx[,type:="untreated"]
## N3RYL[,type:="treated"]


names(untx); names(N3RYL)
N3RYL <- rbind(N3RYL,untx)


## N3RYL$type <- factor(N3RYL$type,
##                       levels=rev(c("treated within 5 years",
##                                    "treated over 5 years ago","untreated")))
N3RYL$type <- factor(N3RYL$type,levels=rev(c("treated","untreated")))


N3RYL <- merge(N3RYL,wrk,by='g_whoregion')

## show_col(colorblind_pal()(4))
clz <- colorblind_pal()(4)


ggplot(N3R[year>=2015],aes(acat,alive.t,fill=year)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('../plots/PostRegionYearLast5.pdf'),w=9,h=5)

## HIV
ggplot(N3R[g_whoregion=='AFR' & year>=2000],aes(acat,alive.h,fill=year)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('../plots/PostRegionYearHIV.pdf'),w=6,h=6)


## HIV
fwrite(N3R[,.(last5=format(round(sum(alive.h)),big.mark = ',',signif=3)),
           by=g_whoregion],
       file=here('figs/hivreg.csv'))

fwrite(N3R[,.(last5=format(round(sum(alive.h)),big.mark = ',',signif=3))],
       file=here('figs/hivglob.csv'))



ggplot(N3R,aes(year,alive.t,fill=acat)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Year') + ylab('Post-tuberculosis treatment, alive 2020')+
  ## theme_minimal()+
  scale_x_discrete(breaks=seq(from=1980,to=2020,by=5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),
        legend.position=c(0.05,0.8),
        legend.text=element_text(size=rel(0.75)),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(fill='Age') 

ggsave(here('../plots/PostRegionYear2.pdf'),w=10,h=7)


load(file=here('../tmpdata/lamap.Rdata'))

## ## age now! NOTE
## N3[,agenow:=age + 2020-year]

## NX <- merge(N3[,.(age=agenow,sex,iso3,alive.t,g_whoregion)],
##             lamap[,.(age,acat=acats)],
##             by='age',all.x=TRUE,all.y=FALSE)

## NX <- NX[!is.na(acat)]
## NX

## NX$acat <- factor(NX$acat,levels=racts,ordered=TRUE)

## tmp <- NX[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion)]



## LYs per person among those living
## ie Years Post TB
## N3[,YPT:=2020-year]
## N3 <- merge(N3,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
## ## N3[,acats:=NULL]
## N3RYLx <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #TODO stats for article
## N3RYL

names(estl)
estl <- merge(estl,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
estl[,YPT:=2020-year]

untxx <- estl[,.(YPT=weighted.mean(YPT,w=alive)),by=.(g_whoregion,sex,acats)] #TODO stats for article

load(file=here('../figdat/N3RYL.Rdata'))

untxx[,type:="untreated"]
## N3RYL[,type:="treated"]
names(untxx); names(N3RYL)
N3RYL <- rbind(N3RYL,untxx)

N3RYL$acats <- factor(N3RYL$acats,levels=racts,ordered=TRUE)

## TODO have overwritten figure 2!
## --- figure 2 TODO ----

## fig2 TODO
GP <- ggplot(N3RYL,aes(x=acats,y=LYS.t,fill=type)) +
  geom_bar(position='dodge',stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age now (years)') + ylab('Mean years lived post-tuberculosis among those alive 2020')+
  scale_fill_manual(values=clz[c(4,1)])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')
GP

ggsave(GP,file=here('figs/Figure2.pdf'),w=9,h=7)
ggsave(GP,file=here('figs/Figure2.eps'),w=9,h=7)
ggsave(GP,file=here('figs/Figure2.png'),w=9,h=7)


names(estl)
estl <- merge(estl,lamap[,.(agenow=age,acats)],by='agenow') #acats now category for age now
estl[,YPT:=2020-year]

untxx <- estl[,.(YPT=weighted.mean(YPT,w=alive)),by=.(g_whoregion,sex,acats)] #TODO stats for article

load(file=here('../figdat/N3RYLx.Rdata'))

untxx[,type:="untreated"]
## N3RYLx[,type:="treated"]
names(untxx); names(N3RYLx)
N3RYLx <- rbind(N3RYLx,untxx)

N3RYLx$acats <- factor(N3RYLx$acats,levels=racts,ordered=TRUE)




## fig2b
GP <- ggplot(N3RYLx,aes(x=acats,y=YPT,fill=type)) +
  geom_bar(position='dodge',stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age now (years)') + ylab('Mean years lived post-tuberculosis among those alive 2020')+
  scale_fill_manual(values=clz[c(4,1)])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')
GP

ggsave(GP,file=here('../plots/Figure2b.pdf'),w=9,h=7)


## version 2c
## LYs per person among those living
## ie Years Post TB
## N3[,acats:=NULL]
## N3 <- merge(N3,lamap[,.(age=age,acats)],by='age') #acats now category for age of TB
## N3RYLx2 <- N3[,.(YPT=weighted.mean(YPT,w=alive.t)),by=.(g_whoregion,sex,acats)] #TODO stats for article
## N3RYLx2

estl[,acats:=NULL]
estl2 <- merge(estl,lamap[,.(age,acats)],by='age') #acats now category for age of TB
estl2[,YPT:=2020-year]

untxx <- estl2[,.(YPT=weighted.mean(YPT,w=alive)),by=.(g_whoregion,sex,acats)] #TODO stats for article

load(file=here('../figdat/N3RYLx2.Rdata'))

untxx2[,type:="untreated"]
## N3RYLx2[,type:="treated"]
names(untxx2); names(N3RYLx2)
N3RYLx2 <- rbind(N3RYLx2,untxx2)

N3RYLx2$acats <- factor(N3RYLx2$acats,levels=racts,ordered=TRUE)

## fig2c
GP <- ggplot(N3RYLx2,aes(x=acats,y=YPT,fill=type)) +
  geom_bar(position='dodge',stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of tubeculosis (years)') + ylab('Mean years lived post-tuberculosis among those alive 2020')+
  scale_fill_manual(values=clz[c(4,1)])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')
GP

ggsave(GP,file=here('../plots/Figure2c.pdf'),w=9,h=7)




## ============== fig3 ====================================
## these are those treated
## NZ <- merge(N3[,.(age=agenow,year,sex,iso3,alive.t,g_whoregion)],
##             lamap[,.(age,acat=acats)],
##             by='age',all.x=TRUE,all.y=FALSE)
## NZ <- NZ[!is.na(acat)]
## NZ

## NZ$acat <- factor(NZ$acat,levels=racts,ordered=TRUE)
## NZ[year>=2015,type:='treated within 5 years']
## NZ[year<2015,type:='treated over 5 years ago']
## tmp <- NZ[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion,type)] #TODO check Figure 3


NU <- merge(estl[,.(age=agenow,year,sex,iso3,alive,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
NU <- NU[!is.na(acat)]
NU
tmp2 <- NU[,.(alive=sum(alive)),by=.(acat,sex,g_whoregion)]
tmp2[,type:='untreated']

## today <- NZ[,.(alive=sum(alive.t)),by=.(acat,sex,type,g_whoregion)]
## today <- rbind(tmp,tmp2)

load(file=here('../figdat/NZ.Rdata'))
load(file=here('../figdat/NU.Rdata'))

today <- rbind(NZ,NU)

today$type <- factor(today$type,
                     levels=rev(c("treated within 5 years",
                              "treated over 5 years ago","untreated")))

today <- merge(today,wrk,by='g_whoregion')

## fig3
ggplot(today,aes(x=acat,y=alive,fill=type)) +
  geom_bar(stat='identity') +
  facet_grid(sex~name) +
  scale_y_continuous(label=absspace) +
  xlab('Age in 2020 (years)') + ylab('Tuberculosis survivors alive 2020')+
  ggthemes::scale_fill_colorblind()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')

ggsave(here('figs/Figure3.pdf'),w=9,h=5)
ggsave(here('figs/Figure3.png'),w=9,h=5)
ggsave(here('figs/Figure3.eps'),w=9,h=5)
