## combining the different notification elements


rm(list=ls())
library(here)

source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/N3.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/N3.Rdata exists!')
}

load(here('../tmpdata/TBA.Rdata'))
load(here('../tmpdata/N2.Rdata'))
load(here('../tmpdata/HH.Rdata'))
load(here('../tmpdata/TBOS.Rdata'))
load(here('../tmpdata/TBH.Rdata'))
load(here('../tmpdata/NNR.Rdata'))      #new sd

## merge against TBA and use npc
TBA[iso3=='ABW' & year==2010]

tbw <- dcast(TBA,iso3 + acat + year ~ Sex,value.var = 'npc')
tbw[,Year:=year]
tbw[,year:=NULL]

N2 <- merge(N2,tbw,by=c('iso3','Year','acat'),all.x = TRUE)

tmp <- N2[iso3=='ZAF']
ggplot(tmp,aes(Year,Female,col=acat)) + geom_line()

N2[,NF:=PFa*Female]                     #differ in Female/Male
N2[,NM:=PMa*Male]
N2 <- N2[!is.na(NF) & !is.na(NM)]

## check (this is restricted to new)
N2[,sum(NF)]/1e6 + N2[,sum(NM)]/1e6     #168m
TBA[,sum(notes)]/1e6                    #168m

N3 <- melt(N2[,.(iso3,year=Year,acat,age,Male=NM,Female=NF)],
           id.vars = c('iso3','year','acat','age'))

names(N3)[5] <- 'sex'
N3[,sum(value)]/1e6                     #168m
N3[,value0:=value]                      #including deaths


## merge TB outcomes here
N3 <- merge(N3,TBOS[,.(iso3,pd,pd.sd)],by='iso3',all.x=TRUE)
N3[,value:=value * (1-pd)]
N3[,value.sd:=value * pd.sd]            #UNC

N3 <- merge(N3,NNR,by='iso3',all.x=TRUE) #
N3[,value.sd:=xfun(value,1,value.sd,rat.sd)] #rat already included in value

## N3[,value.sd:=value * pd.sd] # bug! included 2x
N3[,sum(value)]/1e6                     #162m
N3[,summary(value/value0)]              #deaths out
N3[is.na(value/value0),.(iso3,value,value0)]

## merge in HIV data
N3 <- merge(N3,TBH,by=c('iso3','year'),all.x=TRUE)
N3[is.na(hs),hs:=0]
N3[,range(year)]

## merge against HIV data
lamaps <- merge(lamap,amap,by=c('AgeGrp','acats'),all.x=TRUE)


N3 <- merge(N3,lamaps,by='age',all.x=TRUE)
names(HH)[2] <- 'sex'
N3 <- merge(N3,HH,by=c('iso3','year','sex','age_name'),all.x=TRUE)
N3[is.na(val),c('val','upper','lower'):=0]


N3[,K:= hs * sum(value0)/sum(value0 * val)] #using TB notes before outcome

N3[,H:= K * val]
N3[,H.sd:=K * (upper-lower)/3.92]       #UNC

## safety
N3[H>1,H:=0.9]

## drop names
N3[,K:=NULL]
N3[,c('upper','lower'):=NULL]
N3[,c('val','hs'):=NULL]
N3[,c('acat','age_name'):=NULL]

N3


save(N3,file=here('../tmpdata/N3.Rdata'))
