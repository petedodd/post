## interpolation of demography
rm(list=ls())
library(here)
source(here('Neat/0_utilities.R'))

if(! overwrite ){
  if(file.exists(here('tmpdata/N2.Rdata')))stop('Not running as tmpdata/N2.Rdata exists!')
}


N <- merge(N,amap,by='AgeGrp',all.x=TRUE)


## simple version with TB demography
NS <- N[,.(PopMale=sum(PopMale),PopFemale=sum(PopFemale)),by=.(iso3,Year,acats)]
NSS <- melt(NS[Year %in% 1980:2019],id.vars = c('iso3','Year','acats'))
NSS[,Sex:='Male']
NSS[variable=='PopFemale',Sex:='Female']
NSS[,acat:=acats]
NSS[,year:=Year]
NSS[,c('acats','variable','Year'):=NULL]

save(NSS,file=here('tmpdata/NSS.Rdata'))


N <- N[Year %in% 1980:2019]
N <- N[AgeGrp!='100+']
N <- N[order(iso3,Year,AgeGrp)]
nn <- nrow(N)
N <- N[rep(1:nrow(N),each=5)]
N[,age:=rep(0:99,nrow(N)/100)]
N[,amid:=mean(age),by=.(iso3,Year,AgeGrp)]


## --- example country
tmp <- N[iso3=='ABW' & Year==1980,.(amid,PopMale)]
utmp <- unique(tmp)
v <- approx(x=utmp$amid,y=utmp$PopMale,xout=0:99,rule=2)
v
tmp[,ap:=v$y]
tmp[,age:=v$x]

GP <- ggplot(utmp,aes(amid,PopMale)) +
  geom_line(data=tmp,aes(age,ap))+
  geom_point(data=tmp,aes(age,ap),size=1) +
  geom_point(shape=3,col=2,size=3) +
  xlab('Age') + ylab('Population')+
  theme_classic()

if(plt)ggsave(GP,file=here('plots/interp_5yr.pdf'),w=5,h=5)


## linear extrapolation of population across 5 year age groups
N2 <- N[,{
  vm <- approx(x=amid[ev5],y=PopMale[ev5],xout = age,rule=2)
  vf <- approx(x=amid[ev5],y=PopFemale[ev5],xout = age,rule=2)
  vt <- approx(x=amid[ev5],y=PopTotal[ev5],xout = age,rule=2)
  list(PMa=vm$y,PFa=vf$y,PTa=vf$y,
       PopMale=PopMale,PopFemale=PopFemale,PopTotal=PopTotal,
       acats=acats,age=age,amid=amid)
},by=.(iso3,Year)]
N2[,c('totm','totf','tott'):=list(sum(PMa),sum(PFa),sum(PTa)),by=.(iso3,Year,amid)]
N2[,c('pm','pf','pt'):=list(PMa/totm,PFa/totf,PTa/tott),by=.(iso3,Year,amid)]
N2[,sum(pt),by=.(iso3,Year,amid)]      #check - fraction from each TB age gp in 5 year ag
## multiply to get pop in 5 year age groups
N2[,c('PMa','PFa','PTa'):=list(PopMale*pm,PopFemale*pf,PopTotal*pt)]


GP <- ggplot(N2[iso3=="ZWE" & Year==1980,.(age,PMa)],aes(age,PMa)) +
  geom_point() + geom_line() + xlab('Age') + ylab('Population') + theme_classic()

if(plt)ggsave(GP,file=here('plots/interp_1yr.pdf'),w=5,h=5)

## check
N2[,sum(PTa)]/1e6                       #242
N2[,sum(PopTotal)/5]/1e6                #242

## extraoplating TB patterns to year age groups
## keep
N2 <- N2[,.(iso3,Year,age,acat=acats,PMa,PFa,PTa)]
N2 <- N2[,.(iso3,Year,age,acat,PMa,PFa,PTa)]

save(N2,file=here('tmpdata/N2.Rdata'))
