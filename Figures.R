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


(pcntpopdone <- 1e2*N[Year==2020 & iso3 %in% cnsdone,sum(PopTotal)]/N[Year==2020,sum(PopTotal)])
cat(pcntpopdone,file=here('texto/pcntpopdone.txt'))


## TB props in kids
## treated
load(file=here('../tmpdata/txknum.Rdata'))
load(file=here('../tmpdata/txkden.Rdata'))
load(file=here('../tmpdata/txknumy.Rdata'))
load(file=here('../tmpdata/txkdeny.Rdata'))
## untreated
load(file=here('../tmpdata/utxknum.Rdata'))
load(file=here('../tmpdata/utxkden.Rdata'))
load(file=here('../tmpdata/utxknumy.Rdata'))
load(file=here('../tmpdata/utxkdeny.Rdata'))

## output?
1e2*txknum/txkden #1e2*N3[acat %in% c('0-4','5-14'),sum(alive.t)]/N3[,sum(alive.t)]
1e2*utxknum/utxkden#*estl[acat %in% c('0-4','5-14'),sum(alive)]/estl[,sum(alive)]

1e2*(txknum+utxknum)/(txkden+utxkden)   #10%
1e2*(txknumy+utxknumy)/(txkdeny+utxkdeny)   #26%

length(cnsdone)

## show_col(colorblind_pal()(4))
clz <- colorblind_pal()(4)
load(file=here('../figdat/N3R.Rdata'))


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



## --- figure 2 ----

load(file=here('../figdat/N3RYL.Rdata'))
load(file=here('../figdat/untx.Rdata'))

N3RYL <- rbind(N3RYL,untx)
N3RYL$acat <- factor(N3RYL$acat,levels=racts,ordered=TRUE)
N3RYL$type <- factor(N3RYL$type,levels=c('untreated','treated'),ordered=TRUE)
N3RYL <- merge(N3RYL,wrk,by='g_whoregion')

## fig2 
GP <- ggplot(N3RYL,aes(x=acat,y=LYS.t,fill=(type))) +
  geom_bar(stat='identity') +
  facet_grid(sex~name) +
  scale_y_continuous(label=absspace) +
  xlab('Age of tuberculosis (years)') + ylab('Post-tuberculosis life-years 1980-2020')+
  scale_fill_manual(values=clz[c(1,4)])+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank(),legend.position = 'top')
GP

ggsave(GP,file=here('figs/Figure2.pdf'),w=9.5,h=7)
ggsave(GP,file=here('figs/Figure2.eps'),w=9.5,h=7)
ggsave(GP,file=here('figs/Figure2.png'),w=9.5,h=7)

## fig2b
load(file=here('../figdat/untxx.Rdata'))
load(file=here('../figdat/N3RYLx.Rdata'))

N3RYLx <- rbind(N3RYLx,untxx)
N3RYLx$acats <- factor(N3RYLx$acats,levels=racts,ordered=TRUE)


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
load(file=here('../figdat/untxx2.Rdata'))
load(file=here('../figdat/N3RYLx2.Rdata'))

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

ggsave(here('figs/Figure3.pdf'),w=9.5,h=5)
ggsave(here('figs/Figure3.png'),w=9.5,h=5)
ggsave(here('figs/Figure3.eps'),w=9.5,h=5)
