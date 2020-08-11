## TB notification extrapolation
rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/TBN.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/TBN.Rdata exists!')
}

## load data
## TODO understand!
TB <- fread(here('../indata/TB_notifications_2020-02-24.csv'))
TB <- fread('/Users/pjd/Documents/WHO_TBreports/data2018/TB_notifications_2019-02-14.csv')

## --- notifications -------
sexunks <- grep("sexunk",names(TB),value=TRUE) #
ageunks <- grep("ageunk",names(TB),value=TRUE) #
newrelz <- c(outer(c('newrel_m','newrel_f'),acts,paste0))
newz <- grep("new_",names(TB),value=TRUE)

newz_age <- grep("[sp|ep|sn]",newz,value=TRUE)
newz_age <- grep("\\_[f|m]",newz_age,value=TRUE)
(newz_age <- newz_age[!grepl("*u$",newz_age)])

nmz <- c('iso3','g_whoregion','year',newrelz,newz_age)

tmpw <- TB[,..nmz]

TBM <- melt(tmpw,id=c('iso3','g_whoregion','year'))
TBM[,Sex:='Male']
TBM[grepl('f',variable),Sex:='Female']

## age categories
TBM[,acat:='0-4']                       #this will need special treatment?
TBM[grepl('514',variable),acat:='5-14']
TBM[grepl('1524',variable),acat:='15-24']
TBM[grepl('2534',variable),acat:='25-34']
TBM[grepl('3544',variable),acat:='35-44']
TBM[grepl('4554',variable),acat:='45-54']
TBM[grepl('5564',variable),acat:='55-64']
TBM[grepl('65',variable),acat:='65+']   #this is different
TBM[grepl('15plus',variable),acat:='15+']        #different
TBM[grepl('014',variable),acat:='<15']           #different
TBM[,variable:=NULL]
TBM[,Time:=paste0(5*(year %/% 5),'-',5*(year %/% 5)+5)]
TBM[,range(year)]                                     #back to 1980, but lots missing
TBM <- TBM[!is.na(value)]
TBM[,range(year)]

TBM[,unique(acat)]
TBM

## === India

indc <- TB[iso3=='IND',.(year,c_newinc)]
indc <- indc[order(year)]

GPInd <- ggplot(indc,aes(x=year,y=c_newinc)) + geom_line() +
  scale_y_continuous(label=absspace) +
  geom_vline(xintercept = 1991,col=2) +
  geom_vline(xintercept = 2013,col=2) +
  expand_limits(y=0)
GPInd

## ggsave('plots/IND0c.pdf')


## correct this
## change in reporting meaning
indc2 <- copy(indc)
fac <- indc[year==1992,c_newinc]/indc[year==1991,c_newinc]
indc2[year<=1991,c_newinc:=round(c_newinc * fac)]
## NIKSHAY includes private sector
yz <- 1991:2013
c1 <- indc[year==2014,c_newinc]
c2 <- indc[year==2013,c_newinc]
facz <- seq(from = 1,to=(c1/c2),len=length(yz))
indc2[year %in% yz,c_newinc := round(facz * c_newinc)]

GP <- GPInd + geom_line(data=indc2,lty=2) +
  xlab('Year') + ylab('New and relapse notified TB cases') 

if(plt) ggsave(GP,file=here('../plots/IND0c2.pdf'),w=7,h=5)

## save in corrected
TB[iso3=='IND',.(year,c_newinc)]        #check
TB[iso3=='IND',c_newinc:=indc2$c_newinc]

## === end India

TBM[,sum(value)]/1e6                    #75
TBM[acat=='<15',sum(value)]/1e6         #1.5
TBM[acat=='15+',sum(value)]/1e6         #7.5

## compute a mean ratio of new:rel using ret_rel
TBN <- TB[,.(iso3,c_newinc,year,ret_rel,g_whoregion)]       #new and relapse
summary(TBN)

ggplot(TBN,aes(ret_rel/c_newinc)) + geom_histogram() + facet_wrap(~g_whoregion)


## hierarchy of means UNC
TBN[,rat:=sum(ret_rel,na.rm=TRUE)/sum(c_newinc,na.rm=TRUE),by=iso3] #country mean
TBN[,rat.sd2:=sqrt(rat*(1-rat)/sum(c_newinc,na.rm=TRUE)),by=iso3] #country mean sample unc
TBN[,rat.sd:=sd(ret_rel/c_newinc,na.rm=TRUE),by=iso3] #country variation
TBN[,.(rat,rat.sd/rat,rat.sd2/rat)]                   #first more important
TBN[,rat.sd2:=NULL]                                   #neglect second

TBN[is.na(rat),rat:=sum(ret_rel,na.rm=TRUE)/sum(c_newinc,na.rm=TRUE),by=g_whoregion] #regional mean if missing
TBN[is.na(rat),rat:=sum(ret_rel,na.rm=TRUE)/sum(c_newinc,na.rm=TRUE)] #global mean if missing
TBN[,rat.sd2:=sd(rat,na.rm=TRUE),by=g_whoregion] #regional mean if missing
TBN[,rat.sd1:=mean(rat.sd,na.rm=TRUE),by=g_whoregion] #regional mean if missing
unique(TBN[is.na(rat.sd),.(g_whoregion,rat.sd1,rat.sd2)]) #regional mean if missing
TBN[is.na(rat.sd),rat.sd:=sqrt(rat.sd1^2+rat.sd2^2)]
TBN[,c('rat.sd1',"rat.sd2"):=NULL]      #drop temp vars

TBN
summary(TBN)


## impute cnewinc & restrict to new
tmp <- TBN[year==2017]
tmp[,c('c_newinc','ret_rel'):=NA]
tmp[,year:=2018]
tmp1 <- copy(tmp)
tmp1[,year:=2019]
tmp <- rbind(tmp,tmp1)
TBN <- rbind(TBN,tmp)

TBN[,imp:=FALSE]
TBN[is.na(c_newinc),imp:=TRUE]
TBN[,c_newinc:=as.integer(na_kalman(c_newinc)),by=iso3]
TBN[,imputed:=imp]
if(plt){
  for(reg in TBN[,unique(g_whoregion)]){
    GP <- ggplot(TBN[g_whoregion==reg],aes(year,c_newinc,col=imputed))  +
      geom_point() + #geom_line() +
      facet_wrap(~iso3,scales='free_y') + theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(GP,filename=here(paste0('../plots/notes/IMPcheck_',reg,'.pdf')),w=10,h=10)
  }
}
(imps <- TBN[,table(imputed)])

cat(imps,file=here('texto/imps.txt'))

imps <- TBN[,sum(c_newinc),by=.(imputed)]
(imps[,V1:=1e2*V1/sum(V1)])

fwrite(imps,file=here('texto/imps1.csv'))


TBN[,summary(c_newinc)]
TBN[c_newinc<0,c_newinc:=0]

TBN[,val:=as.double(c_newinc-ret_rel)]  #new
TBN[is.na(val),val:=(1-rat)*c_newinc]   #restrict to new (val)

ggplot(TBN,aes(c_newinc,val,col=g_whoregion)) + geom_point()

TBN[rat.sd>rat,rat.sd:=rat]               #safety

NNR <- unique(TBN[,.(iso3,rat,rat.sd)])


save(NNR,file=here('../tmpdata/NNR.Rdata'))
summary(1e2*NNR$rat)                    #NOTE percent


## TODO consider uncertainty here
## UNC
tmp <- TBN[,.(value=sum(val),
              value.sd=Ssum((rat.sd/(rat+1e-9))*val)),
           by=g_whoregion] #fractional unc for val same as rat
tmp2 <- data.table(g_whoregion='Global',value=tmp[,sum(value)],
                   value.sd=tmp[,Ssum(value.sd)])
t1r2 <- rbind(tmp,tmp2)
t1r2[,quantity:='totnewtx']

## check
t1r2[,.(value.sd/value)]
t1r2[,.(see(value),see(value-value.sd*2),see(value+value.sd*2))]               #

save(t1r2,file=here('../figdat/t1r2.Rdata'))


TBN[,c_newinc0:=c_newinc]
TBN[,c_newinc:=NULL]
TBN[,c_newinc:=val]  #NOTE have here changed c_newinc to be only new!
TBN[,val.sd:=val * (rat.sd/(rat+1e-9))]  #UNC -- seems large, TODO also worry about disaggregation
TBN[,c('ret_rel','imp','imputed'):=NULL]

TBN[,sum(c_newinc,na.rm=FALSE)]/1e6      #171 m
TBN[,range(year)]

## TODO issue with the newer TB notification data
ggplot(TBN[,.(notes=sum(c_newinc,na.rm=TRUE)),by=year],aes(year,notes)) +
  geom_line()

ggplot(TBN[,.(nanotes=sum(is.na(c_newinc))),by=year],aes(year,nanotes)) +
  geom_line()


save(TBN,file=here('../tmpdata/TBN.Rdata'))



## drop duplicate age categories
TBM <- TBM[acat!='<15']
TBM <- TBM[acat!='15+']

save(TBM,file=here('../tmpdata/TBM.Rdata'))

## compare against TBN
TBN2 <- TBM[,.(anotes=sum(value)),by=.(iso3,year)]
TBN2[,range(year)]
TBB <- merge(TBN,TBN2,by=c('iso3','year'),all.x=TRUE)
TBB[,fracdisaggregated:=anotes/c_newinc]


tmp <- TBB[year<=2017,.(fraction=sum(anotes,na.rm = TRUE)/
                sum(c_newinc,na.rm = TRUE)),by=year]
tmp[,quantity:='Proportion of notified TB']
tmp1 <- TBB[year<=2017,.(fraction=mean(!is.na(fracdisaggregated))),by=year]
tmp1[,quantity:='Proportion of countries']
tmp <- rbind(tmp,tmp1)

GP <- ggplot(tmp,aes(year,fraction,lty=quantity)) +
  scale_y_continuous(label=scales::percent)+
  geom_line() + ylab('Proportion disaggregated by age') + xlab('Year')+
  theme_classic() + theme(legend.position = c(0.2,0.9)) +
  ggpubr::grids()

if(plt)ggsave(GP,file=here('../plots/AgeDisagByB.pdf'),w=7,h=5)
