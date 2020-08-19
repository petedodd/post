## extrapolation of HIV in TB
rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/TBH.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/TBH.Rdata exists!')
}

## load data
## WHO TB notifications
TB <- fread(here('../indata/TB_notifications_2019-02-14.csv'))
## WHO TB burden estimates
est <- fread(here('../indata/TB_burden_countries_2020-02-24.csv'))
## UNAIDS AIDSinfo ART coverage 4/8/2019
AC <- fread(here("../indata/Treatment_cascade_PLHIV_Allages.csv"),skip=1)

## replace this with estimated quantity
names(est)
est[,.(year,iso3,e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]

WH <- est[iso3 %in% hivcountries,.(year,iso3,
                                   e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]

## extrapolation see HIV
WH

xtrayr <- c(1980:1999,2019:2020)
WH2 <- rbind(WH,data.table(year=xtrayr,iso3=rep(unique(WH$iso3),
                                                   each = length(xtrayr)),
                           e_tbhiv_prct=NA,e_tbhiv_prct_lo=NA,e_tbhiv_prct_hi=NA))

WH2 <- WH2[order(iso3,year)]
WH2[year<1990,e_tbhiv_prct:=1e-10]

WH2[iso3=='BWA']

yrz <- unique(WH2$year)                 #last 2 years added

WH2[,e_tbhiv_prct:=na_kalman(e_tbhiv_prct),by=iso3]


GP <- ggplot(WH2,aes(year,e_tbhiv_prct)) +
  annotate('rect',alpha=.2,xmin=2000,xmax=2018,ymin=0,ymax=Inf)+
  geom_line(linetype=1) +
  xlab('Year')+
  ylab('HIV prevalence in TB (%)')+
  ylim(c(0,100))+
  rot45+
  facet_wrap(~iso3)#,scales='free_y')
GP

if(plt) ggsave(GP,file=here::here('../plots/HIVinTBinterp.pdf'),w=10,h=10)


TBH <- WH2[,.(iso3,year,hs=e_tbhiv_prct/100)] #using new version

ggplot(WH2,aes(year,e_tbhiv_prct,col=is.na(e_tbhiv_prct_lo))) + geom_point() + facet_wrap(~iso3)



TB[!is.na(hiv_art),.(iso3,year,newrel_hivpos,newrel_art,hivtest_pos,hiv_art)]

TBH <- merge(TBH,TB[!is.na(hiv_art),.(iso3,year,hsa=hiv_art/hivtest_pos)],
             by=c('iso3','year'),all.x = TRUE,all.y=FALSE)

TBH[year<=2000,hsa:=1e-10]
TBH[is.nan(hsa),hsa:=NA]
TBH[,hsa:=na_kalman(hsa),by=iso3]
TBH[hsa<0,hsa:=1e-10]
TBH[hsa>.90,hsa:=.90]
TBH[,hsa:=smooth(hsa),by=iso3]



## working with UNAIDS data
cz <- c(1,1+seq(from=1,by=3,len=10))
vz <- paste0('V',cz)
AC <- AC[,..vz]
names(AC)[1] <- 'country'
names(AC)[2:ncol(AC)] <- 2010:2019
AC

AC[,iso3:=countrycode::countrycode(as.character(country),
                                   origin='country.name',destination='iso3c')]
AC[is.na(iso3)]
AC <- AC[country!='Global']
AC[country=='Eswatini',iso3:='SWZ']
ACN <- cbind(iso3=AC[,iso3],AC[,lapply(.SD,as.numeric),.SDcols=2:(ncol(AC)-1)])
ACN <- ACN[iso3 %in% hivcountries]
ACN[iso3=='SWZ',`2019`:=95]
ACM <- melt(ACN,id='iso3')
ACM[,year:=as.integer(as.character(variable))]

TBH <- merge(TBH,ACM[,.(iso3,year,haa=value/1e2)],
             by=c('iso3','year'),all.x=TRUE,all.y=FALSE)

TBH[year<=2000,haa:=1e-10]
TBH[iso3=='PRT',haa:=hsa]               #all NA
TBH[,haa:=na_kalman(haa),by=iso3]

TBH[haa<0,haa:=1e-10]
TBH[haa>.90,haa:=.90]

GP <- ggplot(TBH,aes(year,1e2*hs)) +
  annotate('rect',alpha=.2,xmin=2000,xmax=2018,ymin=0,ymax=Inf)+
  geom_line(linetype=1) +
  geom_line(linetype=1,aes(year,1e2*hsa),col=2) +
  geom_line(linetype=1,aes(year,1e2*haa),col=4) +
  xlab('Year')+
  ylab('Percent')+
  ylim(c(0,100))+
  rot45+
  facet_wrap(~iso3)#,scales='free_y')
GP

if(plt)ggsave(GP,file=here::here('../plots/HIVinTBinterp2.pdf'),w=10,h=10)


save(TBH,file=here('../tmpdata/TBH.Rdata'))
