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

## length(rracts2)
## data.frame(0:99,rracts2)

## load data
load(here('N3.Rdata'))
## TB <- fread('/Users/pjd/Documents/WHO_TBreports/data2018/TB_notifications_2019-02-14.csv')


## N3 <- merge(N3,unique(TB[,.(iso3,g_whoregion)]),by='iso3',all.x = TRUE)

## proportion covered
cnsdone <- N3[,(unique(iso3))]
cat(cnsdone,file=here('out/cnsdone.txt'))

est <- fread('/Users/pjd/Documents/WHO_TBreports/data2019/TB_burden_countries_2020-02-24.csv')
est <- est[,.(iso3,year,e_inc_num)]

estr <- est[year==2017]
(pcnttbdone <- 1e2*estr[iso3 %in% cnsdone,sum(e_inc_num)]/estr[,sum(e_inc_num)])
cat(pcnttbdone,file=here('texto/pcnttbdone.txt'))

load(here('../WPP/N_simple.Rdata'))           #5 year age groups
N[Year==2020,sum(PopTotal)]/1e6

(pcntpopdone <- 1e2*N[Year==2020 & iso3 %in% cnsdone,sum(PopTotal)]/N[Year==2020,sum(PopTotal)])
cat(pcntpopdone,file=here('texto/pcntpopdone.txt'))


length(cnsdone)

1
## ## --- Results ---
## ## reporting
N3S <- N3[,.(alive=sum(alive),LYS=sum(LYS)),by=.(iso3,year,sex,acat)]
N3S[,sum(alive)]/1e6                  #98m
N3S[,sum(LYS)]/1e9                    #1.97 bn

N3S <- N3[,.(alive=sum(alive.t),LYS=sum(LYS.t)),by=.(iso3,year,sex,acat)]
N3S[,sum(alive)]/1e6                  #98m
N3S[,sum(LYS)]/1e9                    #2.0 bn

N3S.h <- N3[,.(alive=sum(alive.h),LYS=sum(LYS.h)),by=.(iso3,year,sex,acat)]
N3S.h[,sum(alive)]/1e6                    #3.3 mln
N3S.h[,sum(LYS)]/1e6                    #54 mln

## regional aggregation
N3S$acat <- factor(N3S$acat,levels=racts,ordered=TRUE)
yy <- N3S[,unique(year)]
yy <- sort(yy)
N3S$year <- factor(N3S$year,levels=yy,ordered=TRUE)
## N3S <- merge(N3S,unique(TBA[,.(iso3,g_whoregion)]),by='iso3',all.x = TRUE)

names(N3)
N3[,unique(age)]
racts

N3[,acat:=rracts2[1+age]]

N3$acat <- factor(N3$acat,levels=racts,ordered=TRUE)

## regional summary
N3R <- N3[,.(alive.t=sum(alive.t),LYS.t=sum(LYS.t),
             alive.h=sum(alive.h),LYS.h=sum(LYS.h)),
          by=.(g_whoregion,year,sex,acat)]

N3R$acat <- factor(N3R$acat,levels=racts,ordered=TRUE)
N3R$year <- factor(N3R$year,levels=yy,ordered=TRUE)

N3S <- merge(N3S,unique(N3[,.(iso3,g_whoregion)]),by='iso3')

tmp <- N3S[,.(alive=sum(alive)),by=.(acat,sex,year,g_whoregion)]
tmp[year==2010]

N3R[g_whoregion=='AMR' & acat=='0-4' & sex=='Male',sort(year)]

ggplot(N3R,aes(acat,alive.t,fill=year)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment, alive 2020')+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYear.pdf'))


## LY version
N3RY <- N3R[,.(LYS.t=sum(LYS.t)),by=.(g_whoregion,sex,acat)]

ggplot(N3RY,aes(acat,LYS.t)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment LYs to 2020')+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYearLY.pdf')) #TODO LYV

load(here('estl.Rdata'))


N3R[,last5:=ifelse(year>=2015,'yes','no')]
N3RYL <- N3R[,.(LYS.t=sum(LYS.t)),by=.(last5,g_whoregion,sex,acat)]

untx <- estl[,.(LYS.t=sum(LYS)),by=.(g_whoregion,sex,acat)]

untx[,type:="untreated"]

N3RYL[,type:="treated over 5 years ago"]
N3RYL[last5=='yes',type:="treated within 5 years"]

N3RYL[,last5:=NULL]

names(untx); names(N3RYL)
N3RYL <- rbind(N3RYL,untx)


N3RYL$type <- factor(N3RYL$type,
                      levels=rev(c("treated within 5 years",
                                   "treated over 5 years ago","untreated")))

save(N3RYL,file=here('N3RYL.Rdata'))

N3RYL <- merge(N3RYL,wrk,by='g_whoregion')


## USE?
## fig3
ggplot(N3RYL,aes(x=acat,y=LYS.t,fill=type)) +
  geom_bar(stat='identity') +
  facet_grid(sex~name) +
  scale_y_continuous(label=absspace) +
  xlab('Age of tubeculosis (years)') + ylab('Post-tuberculosis life-years 1980-2020')+
  ggthemes::scale_fill_colorblind()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')

## ggsave('plots/3bar1LY.pdf',w=12,h=6)
ggsave(here('Neat/figs/Figure3.pdf'),w=9,h=5)
ggsave(here('Neat/figs/Figure3.png'),w=9,h=5)
ggsave(here('Neat/figs/Figure3.eps'),w=9,h=5)


ggplot(N3R[year>=2015],aes(acat,alive.t,fill=year)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYearLast5.pdf'))

## HIV
ggplot(N3R[g_whoregion=='AFR' & year>=2000],aes(acat,alive.h,fill=year)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYearHIV.pdf'))


## HIV
fwrite(N3R[,.(last5=format(round(sum(alive.h)),big.mark = ',',signif=3)),
           by=g_whoregion],
       file=here('plots/hivreg.csv'))

fwrite(N3R[,.(last5=format(round(sum(alive.h)),big.mark = ',',signif=3))],
       file=here('plots/hivglob.csv'))



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

ggsave(here('plots/PostRegionYear2.pdf'),w=10,h=7)
ggsave(here('plots/Figure2.pdf'))#,w=10,h=7)


## looking at SEA
ggplot(N3R[g_whoregion=='SEA'],aes(year,alive.t,fill=acat)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('year') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=4))

sea <- N3[g_whoregion=='SEA',
          .(value=sum(value),alive.t=sum(alive.t)),
          by=.(iso3,year)]

## IND
ggplot(sea[iso3=='IND'],aes(year,alive.t)) +
  geom_line() +
  facet_wrap(~iso3) +
  scale_y_continuous(label=absspace) +
  xlab('year') + ylab('Post-TB treatment, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('plots/IND1.pdf')

## ## IND
## ggplot(sea[iso3=='IND'],aes(year,value)) +
##   geom_line() +
##   geom_line(data=indc,aes(year,c_newinc),lty=2) + 
##   facet_wrap(~iso3) +
##   scale_y_continuous(label=absspace) +
##   ## annotate('rect',fill='red',col=NA,alpha=.2,xmin=1989,xmax=1993,ymin=0,ymax=Inf)+
##   xlab('year') + ylab('Notifications (less deaths on treatment)')+
##   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + expand_limits(y=0)

## ggsave('plots/IND2key.pdf')



## N3S[,range(year)]

load(file=here('out/lamap.Rdata'))

## age now! NOTE
N3[,agenow:=age + 2020-year]

NX <- merge(N3[,.(age=agenow,sex,iso3,alive.t,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)

NX <- NX[!is.na(acat)]
NX

NX$acat <- factor(NX$acat,levels=racts,ordered=TRUE)

tmp <- NX[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion)]


ggplot(tmp,aes(acat,alive)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age today') + ylab('Tuberculosis treatment survivors alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYear_today.pdf'))
ggsave(here('plots/Figure1.pdf'))



ggplot(tmp,aes(acat,alive)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age today') + ylab('Post-TB, un-treated, alive 2020')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/PostRegionYear_todayNoTx.pdf'))



## -------------
tmp

N3R[year>=2015,type:='treated within 5 years']
N3R[year<2015,type:='treated over 5 years ago']

night <- tmp
night[,type:='untreated']

## tmp2 <- NX[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion)]
## tmp[,type:='treated']
tmp2 <- N3R[,.(alive=sum(alive.t)),by=.(acat,sex,type,g_whoregion)]
night <- rbind(night,tmp2)

night$type <- factor(night$type,
                     levels=rev(c("treated within 5 years",
                              "treated over 5 years ago","untreated")))

save(night,file=here('night.Rdata'))

load(here('night.Rdata'))

night2 <- night[,.(alive=sum(alive)),by=.(g_whoregion,type)]
night2$type <- factor(night2$type,
                     levels=rev(c("treated within 5 years",
                              "treated over 5 years ago","untreated")))


ggplot(night2,aes(x=g_whoregion,y=alive,fill=(type))) +
  geom_bar(stat='identity',width = 1) +
  coord_polar(theta = "x", start = 11) +
  scale_y_continuous(label=absspace)+
  xlab('') + ylab('Number of tuberculosis survivors')+
  ggthemes::scale_fill_colorblind()+
  theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(here('plots/Polar.pdf'),w=6,h=8)



## USE
ggplot(night2,aes(x=g_whoregion,y=alive,fill=type)) +
  geom_bar(stat='identity') +
  scale_y_continuous(label=absspace)+
  xlab('') + ylab('Number of tuberculosis survivors')+
  ggthemes::scale_fill_colorblind()+
  theme_classic()+
  ggpubr::grids()+
  theme(legend.position = c(0.2,0.9),legend.title=element_blank())
  ## theme(legend.position = 'bottom',legend.title=element_blank())

ggsave(here('plots/PolarNotA.pdf'),w=6,h=6)


ggplot(night2,aes(x=g_whoregion,y=alive,fill=type)) +
  geom_bar(stat='identity',position='dodge') +
  scale_y_continuous(label=absspace)+
  xlab('') + ylab('Number of tuberculosis survivors')+
  ggthemes::scale_fill_colorblind()+
  theme_classic()+
  ggpubr::grids()+
  ## theme(legend.position = 'bottom',legend.title=element_blank())
  theme(legend.position = c(0.2,0.9),legend.title=element_blank())

ggsave(here('plots/PolarNotB.pdf'),w=6,h=6)


night2



night$type <- factor(night$type,levels=rev(levels(night$type)))

## USE
ggplot(night,aes(x=acat,y=alive,fill=type)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment numbers to 2020')+
  ggthemes::scale_fill_colorblind()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
 

ggsave(here('plots/PolarC.pdf'))
ggsave('plots/3bar1.pdf',w=12,h=6)
ggsave('plots/3bar1.png',w=12,h=6)




ggplot(night,aes(x=acat,y=alive,fill=type)) +
  geom_bar(stat='identity',position = 'dodge') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment numbers to 2020')+
  ggthemes::scale_fill_colorblind()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('plots/3bar2.pdf'),w=12,h=6)




night2[,sum(alive),by=type]             #97 vs 52 (26)

night[,tot:=sum(alive),by=.(acat,sex,g_whoregion)]
night[,frac:=alive/tot]
night


## NOTE - maybe this is correct - NB age at TB
ggplot(night,aes(x=acat,y=frac,fill=type)) +
  geom_bar(stat='identity') +
  facet_grid(sex~g_whoregion) +
  scale_y_continuous(label=absspace) +
  xlab('Age of TB') + ylab('Post-TB treatment numbers to 2020')+
  ggthemes::scale_fill_colorblind()+
  theme_classic()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(here('plots/3bar3.pdf'),w=12,h=6)


## fig2
## these are those treated
NZ <- merge(N3[,.(age=agenow,year,sex,iso3,alive.t,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
NZ <- NZ[!is.na(acat)]
NZ

NZ$acat <- factor(NZ$acat,levels=racts,ordered=TRUE)
NZ[year>=2015,type:='treated within 5 years']
NZ[year<2015,type:='treated over 5 years ago']
tmp <- NZ[,.(alive=sum(alive.t)),by=.(acat,sex,g_whoregion,type)] #TODO check Figure 3


NU <- merge(estl[,.(age=agenow,year,sex,iso3,alive,g_whoregion)],
            lamap[,.(age,acat=acats)],
            by='age',all.x=TRUE,all.y=FALSE)
NU <- NU[!is.na(acat)]
NU
tmp2 <- NU[,.(alive=sum(alive)),by=.(acat,sex,g_whoregion)]
tmp2[,type:='untreated']

today <- NZ[,.(alive=sum(alive.t)),by=.(acat,sex,type,g_whoregion)]
today <- rbind(tmp,tmp2)

today$type <- factor(today$type,
                     levels=rev(c("treated within 5 years",
                              "treated over 5 years ago","untreated")))

today <- merge(today,wrk,by='g_whoregion')

## fig2
ggplot(today,aes(x=acat,y=alive,fill=type)) +
  geom_bar(stat='identity') +
  facet_grid(sex~name) +
  scale_y_continuous(label=absspace) +
  xlab('Age in 2020 (years)') + ylab('Tuberculosis survivors alive 2020')+
  ggthemes::scale_fill_colorblind()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')

ggsave(here('Neat/figs/Figure2.pdf'),w=9,h=5)
ggsave(here('Neat/figs/Figure2.png'),w=9,h=5)
ggsave(here('Neat/figs/Figure2.eps'),w=9,h=5)


## product plots

## library(productplots)

## night2

## library(treemap)


## treemap(night2,index=c("g_whoregion", "type"),
##         vSize="alive",title = "",
##         position.legend = 'right',
##         fontsize.labels = c(0,10))


## treemap(night2,index=c("type","g_whoregion"),
##         vSize="alive",title = "",
##         position.legend = 'right',
##         fontsize.labels = c(0,10))

## treemap(night2,index=c("type","g_whoregion"),
##         vSize="alive",title = "",
##         position.legend = 'right',
##         align.labels=list(
##           c("left", "top"),
##           c("center", "center")
##         ),
##         fontsize.labels = c(20,10))


## night3 <- night2

## night3[,supertype:='untreated']
## night3[type!='untreated',supertype:='treated']
## night3[type=='untreated',supertype:='']

## pdf('plots/tm2.pdf')
## treemap(night3,index=c("supertype","type"),
##         vSize="alive",title = "",
##         position.legend = 'none',
##         align.labels=list(
##           c("left", "top"),
##           c("center", "center")
##         ),
##         fontsize.labels = c(20,10))
## dev.off()

## pdf('plots/tm3.pdf')

## treemap(night3,index=c("g_whoregion","supertype","type"),
##         vSize="alive",title = "",
##         position.legend = 'none',
##         align.labels=list(
##           c("left", "top"),
##           c("center", "center")
##         ),
##         fontsize.labels = c(20,10))

## dev.off()

## save(night3,file='night3.Rdata')



## mp <- prodplot(data=night2,alive ~ type+g_whoregion,divider=mosaic()) + aes(fill=type)
## mp

## mp <- mp + theme(axis.line=element_blank(),
##                  axis.text.x=element_text(angle=90),
##                  axis.text.y=element_text(),
##                  axis.ticks=element_blank(),
##                  axis.title.x=element_blank(),
##                  axis.title.y=element_blank(),
##                  legend.title=element_blank(),
##                  legend.position="none",
##                  panel.background=element_blank(),
##                  panel.border=element_blank(),
##                  panel.grid.major=element_blank(),
##                  panel.grid.minor=element_blank(),
##                  plot.background=element_blank()) +
##   ggtitle('Venn')
## mp

## ===== TODO =========
## plots to visualise the temporal and age missing data
## stability of notification patterns over time TODO plotNA?
## plotNA.distribution(x)
##  statsNA(x)


## TB props in kids
1e2*N3[acat %in% c('0-4','5-14'),sum(alive.t)]/N3[,sum(alive.t)]
1e2*estl[acat %in% c('0-4','5-14'),sum(alive)]/estl[,sum(alive)]

1e2*(estl[acat %in% c('0-4','5-14'),sum(alive)] + N3[acat %in% c('0-4','5-14'),sum(alive.t)]) / (N3[,sum(alive.t)]+estl[,sum(alive)])
## 10%

1e2*(estl[acat %in% c('0-4','5-14'),sum(LY)] + N3[acat %in% c('0-4','5-14'),sum(LY)]) / (N3[,sum(LY)]+estl[,sum(LY)])
## 26%
