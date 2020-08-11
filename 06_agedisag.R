## age disaggregation of TB notification

rm(list=ls())
library(here)
source(here('Neat/0_utilities.R'))

if(! overwrite ){
  fn <- here('tmpdata/TBA.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/TBA.Rdata exists!')
}

load(here('tmpdata/TBM.Rdata'))
load(here('tmpdata/TBN.Rdata'))
load(here('tmpdata/NSS.Rdata'))

## ---- disaggregated data
TBM$acat <- factor(TBM$acat,levels=TBM[,unique(acat)],ordered=TRUE)

## smooth 
TBM <- TBM[,.(value=sum(value)),by=.(iso3,g_whoregion,year,Sex,acat)] #NOTE agggregating


ggplot(TBM[iso3=='ZWE' & year == 2015],
       aes(acat,value,col=Sex,group=Sex)) + geom_line()

TBM[,ncats:=.N,by=.(iso3,year)]

## make averaged pattern by age and sex
YD <- TBM[ncats==16,.(iso3,g_whoregion,Sex,acat,value)]
str(YD)
YD <- YD[,.(prop=1.0*sum(value)),by=.(iso3,Sex,acat,g_whoregion)]
YD[,tot:=sum(prop),by=iso3]
YD[,prop:=prop/tot]
YD[iso3=='AFG']
YD[,test:=sum(prop),by=iso3]
YD <- YD[is.finite(test)]
YD[,test:=NULL]
(nodat <- setdiff(TBM[,unique(iso3)],YD[,unique(iso3)])) #8

cat(nodat,file=here('post/texto/notdat.txt'))


TBN[,.(mx=max(year),mn=min(year)),by=iso3]
YD                                      #split by TB age categories

## TBA
## TBA <- merge(TBN,YD,by='iso3',all.y=TRUE,allow.cartesian = TRUE)
TBA <- YD[rep(1:nrow(YD),each=length(1980:2019))]
TBA[,year:=rep(1980:2019,nrow(YD))]

TBA <- merge(TBA,TBN,by=c('iso3','year','g_whoregion'),all.x=TRUE)
TBA[,sum(prop),by=.(iso3,year)]
TBA[,notes:=c_newinc * prop]

TBA[,range(year)]

TBA <- TBA[!is.na(c_newinc)]

TBN[,sum(c_newinc,na.rm=TRUE)]/1e6      #170m
TBA[,sum(notes)]/1e6                    #nrly 
1e2*TBA[,sum(notes)]/TBN[,sum(c_newinc,na.rm=TRUE)]
TBA[,c('tot','c_newinc'):=NULL]
TBA[,prop:=NULL]
TBA$acat <- factor(TBA$acat,levels=racts,ordered=TRUE)

## TBA
## save(TBA,file=here('tmpdata/TBA.Rdata'))  #TB incidence disaggregated by TB age/sex


tmp <- TBA[,.(notes=sum(notes)),by=.(g_whoregion,Sex,acat)]

GP <- ggplot(tmp,aes(acat,notes,col=Sex,group=paste(g_whoregion,Sex))) +
  geom_line() +
  facet_wrap(~g_whoregion) +
  xlab('Age category') + ylab('Cumulative TB notifications') +
  scale_y_continuous(label=absspace)+
  rot45

if(plt)ggsave(GP,file=here('plots/NotesAgePattern.pdf'),h=7,w=10)


TBA <- merge(TBA,NSS,by=c('iso3','Sex','acat','year'),all.x=TRUE) #demography

TBA[,npc:=notes/value]
TBA$acat <- factor(TBA$acat,levels=racts,ordered=TRUE)


tmp <- TBA[,.(notes=mean(npc,na.rm=TRUE)),by=.(g_whoregion,Sex,acat)]


GP <- ggplot(tmp,aes(acat,notes,col=Sex,group=paste(g_whoregion,Sex))) +
  geom_line() +
  xlab('Age category') + ylab('Per capita TB notifications') +
  facet_wrap(~g_whoregion)

if(plt)ggsave(GP,filename=here('plots/npc.pdf'))

TBA

save(TBA,file=here('tmpdata/TBA.Rdata'))   #version with number in age group, popn, npc
