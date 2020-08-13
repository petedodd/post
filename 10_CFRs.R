## Applying CFRs to those untreated

rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/estl.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/N3.Rdata exists!')
}

## load data
load(here('../tmpdata/estg.Rdata'))
load(here('../tmpdata/RRcdr.Rdata'))
load(here('../tmpdata/N3.Rdata'))
load(here('../tmpdata/TBH.Rdata'))


## NOTE this version TBN c_newinc is just new. Check OK with gap
if(! 'acat' %in% names(N3)) N3 <- merge(N3,lamap[,.(age,acat=acats)],by='age') #merge in coarse cats if not already there


N3[,totn:=sum(value),by=.(iso3,year,acat,sex)] #notifications across each coarse cat
N3[,sum(totn)]/1e6

## want to spread the gap proportional to the npc
N3[,totyr:=sum(value)+1e-9,by=.(iso3,year)] #notifications across country years
N3[,propnote:=value/totyr] #proportion of notifications in each age/sex

N3[,sum(propnote),by=.(iso3,year)]      #check
N3[,summary(propnote)]
N3[is.na(propnote)]                     #no cases

## long form estimates with age/sex: making estl
length(unique(est$iso3));length(unique(N3$iso3));
missed <- setdiff(est[,unique(iso3)],N3[,unique(iso3)])
## merge
estl <- merge(est,N3[,.(iso3,year,sex,age,acat,propnote)],
              by=c('iso3','year'),all.y=TRUE)

## global average
propnoteav <- estl[,.(pnav=mean(propnote)),by=.(sex,age)]
propnoteav[,pnav:=pnav/sum(pnav)]
propnoteav[,sum(pnav)]

length(unique(estl$iso3))
estl[is.na(gap)]
estl[!is.finite(gap),unique(iso3)]      #just SDN - no longer a country past 2000
estl <- estl[is.finite(gap)]             #drop these country years

## add back countries lost in merge using global average propnote
est[iso3 %in% missed,sum(gap)]
setdiff(names(estl),names(est))
xtra <- estl[iso3=='AFG' &  year==2000,.(sex,age,acat)]
xtra <- merge(xtra,propnoteav[,.(propnote=pnav,sex,age)],by=c('sex','age'))
cy <- nrow(est[iso3 %in% missed])       #country years missed
xtra2 <- est[iso3 %in% missed][rep(1:cy,each=nrow(xtra))]
xtra2 <- cbind(xtra2,xtra)
xtra2[iso3=='WLF' & year==2010]
xtra2[iso3=='WLF' & year==2010,sum(propnote)]

estl <- rbind(estl,xtra2)

## also add to countries with 0 cases to not muck up below
estl <- merge(estl,propnoteav,by=c('sex','age'),all.x=TRUE)
estl[,tprop:=sum(propnote),by=.(iso3,year)]
estl[tprop==0,propnote:=pnav]           #use average for these countries

## checks
est[year==2010,sum(gap,na.rm=FALSE)]
estl[year==2010,sum(propnote*gap,na.rm=FALSE)]

## merge in relative CDRs
estl <- merge(estl,RRcdr,by='iso3',all.x=TRUE)
summary(estl[,.(RRu5,RR514)])

## do by renormalizing propnote after scaling
estl[acat=='0-4',propnote:=propnote * RRu5]
estl[acat=='5-14',propnote:=propnote * RR514]
estl[,tprop:=sum(propnote),by=.(iso3,year)]
estl[,propnote:=propnote/tprop]
estl[,tprop:=sum(propnote),by=.(iso3,year)] #checks
estl[,summary(tprop)]

## checks
est[year==2010,sum(gap,na.rm=FALSE)]
estl[year==2010,sum(propnote*gap,na.rm=FALSE)]


## disaggregate
estl[,gapl:=gap*propnote]
## uncertainty
## F = sd/X
## F_i = sd_i/X_i = F / (sqrt(sum_jp_j^2)
estl[,tprop:=sqrt(sum(propnote^2)),by=.(iso3,year)]
estl[,summary(tprop)]
estl[,gapl.sd:=(gap.sd/(gap+1e-9)) * gapl / tprop]
estl[,summary(gapl.sd)]

## can ditch tprop and pnav now
estl[,c('tprop','pnav'):=NULL]

## checks
estl[,sum(gapl)]
est[,sum(gap)]
estl[,summary(gapl)]
estl[,Ssum(gapl.sd),by=year]
est[,Ssum(gap.sd),by=year]
estl[,Ssum(gapl.sd)]
est[,Ssum(gap.sd)]

## gapl sense testing
N3[year==2017,sum(value)]/1e6
estl[year==2017,sum(gapl)]/1e6
estl[,sum(gapl)/1e6,by=year]

ggplot(estl[,.(gap=sum(gapl,na.rm=TRUE)/1e6),by=year],aes(year,gap)) + geom_line() + scale_y_continuous(label=absspace)

ggplot(estl[,.(gap=sqrt(sum(gapl.sd^2,na.rm=TRUE))/sum(gapl,na.rm=TRUE)),by=year],aes(year,gap)) + geom_line() + scale_y_continuous(label=absspace)


## NOTE gap includes relapse
## TODO restrict to new


## what proportion of those on ART have been so for less than 1 year?
curve((x-1+exp(-x))/(x^2/2),from=0,to=20) #approximation with linear incidence
TBH[,p1:=pmax(0,year-2000)]
TBH[,p1:=(p1-1+exp(-p1))/(p1^2/2)]
TBH[year<=2000,p1:=1]

## 0.78 (0.65-0.94) HIV+/ART-   WHO appendix ref 40
hcfr <- 0.78
hcfr.sd <- abs(0.65-0.94)/3.92

## 0.62 (0.39-0.86) HIV+/ART<1y WHO appendix
hcfr1 <- 0.62
hcfr1.sd <- abs(0.39-0.86)/3.92

## 0.49 (0.31-0.70) HIV+/ART>1y WHO appendix
hcfrl <- 0.49
hcfrl.sd <- abs(0.31-0.70)/3.92

## weighted average CFR for untreated TB in PLHIV
TBH[,hfr:=haa * (p1*hcfr1 + (1-p1)*hcfrl) + (1-haa) * hcfr]
TBH[,hfr.sd:=sqrt(haa^2 * (p1^2*hcfr1.sd^2 + (1-p1)^2*hcfrl.sd^2) +
                  (1-haa)^2 * hcfr.sd^2)]

## 0.43 (0.28-0.53) WHO appendix
acfr <- 0.43
acfr.sd <- abs(0.28-0.53)/3.92

## 43·6%, 95% CI 36·8–50·6 Jenkins et al
kcfr <- 43.6/1e2
kcfr.sd <- abs(36.8-50.6)/392

## 14·9%, 11·5–19·1 Jenkins et al
ccfr <- 14.9/1e2
ccfr.sd <- abs(11.5-19.1)/392

## merge in HIV data
estl <- merge(estl,N3[,.(iso3,year,age,sex,H,H.sd)],
              by=c('iso3','year','age','sex'),
              all.x=TRUE,all.y = FALSE)

## merge in HIV CFR
estl <- merge(estl,TBH[,.(iso3,year,hfr,hfr.sd)],
              by=c('iso3','year'),all.x=TRUE,all.y=FALSE)
estl[is.na(hfr),c('hfr','hfr.sd'):=0]


## split gap by HIV
estl[,gapl.h:= H * gapl]
estl[,gapl.0:= (1-H) * gapl]
estl[,gapl.h.sd:= xfun(H,gapl,H.sd,gapl.sd)]
estl[,gapl.0.sd:= xfun((1-H),gapl,H.sd,gapl.sd)]
## estl[,gapl.h.sd:= sqrt(H^2 * gapl.sd^2+H.sd^2 * gapl.sd^2+H.sd^2 * gapl^2)]
## estl[,gapl.0.sd:= sqrt((1-H)^2 * gapl.sd^2+H.sd^2 * gapl.sd^2+H.sd^2 * gapl^2)]


## apply CFRs
## HIV -ve
estl[,gapls.0:=gapl.0 * (1-acfr)]
estl[acat=='0-4',gapls.0:=gapl.0 * (1-kcfr)] #kids
estl[acat=='5-14',gapls.0:=gapl.0 * (1-ccfr)] #kids
## HIV +ve
estl[,gapls.h:=gapl.h * (1-hfr)]
## total
estl[,gapls:=gapls.h + gapls.0]

## uncertainty
estl[,gapls.0.sd:=xfun(gapl.0,1-acfr,gapl.0.sd,acfr.sd)]
estl[acat=='0-4',gapls.0.sd:=xfun(gapl.0,1-kcfr,gapl.0.sd,kcfr.sd)]
estl[acat=='5-14',gapls.0.sd:=xfun(gapl.0,1-ccfr,gapl.0.sd,ccfr.sd)]
estl[,gapls.h.sd:=xfun(gapl.h,1-hfr,gapl.h.sd,hfr.sd)]

estl[,gapls.sd:=sqrt(gapls.h.sd^2 + gapls.0.sd^2)]


summary(estl[,.(gapls.sd)])

## estl[,gapls.0.sd:=sqrt(gapl.0^2 * acfr.sd^2 + gapl.0.sd^2 * acfr.sd^2 + gapl.0.sd^2 * (1-acfr)^2)]
## estl[acat=='0-4',gapls.0.sd:=sqrt(gapl.0^2 * kcfr.sd^2 + gapl.0.sd^2 * kcfr.sd^2 + gapl.0.sd^2 * (1-kcfr)^2)]
## estl[acat=='5-14',gapls.0.sd:=sqrt(gapl.0^2 * ccfr.sd^2 + gapl.0.sd^2 * ccfr.sd^2 + gapl.0.sd^2 * (1-ccfr)^2)]
## estl[,gapls.h.sd:=sqrt(gapl.h.sd^2 * (1-hfr)^2 + gapl.h.sd^2 * hfr.sd^2 + gapl.h^2 * hfr.sd^2)]

save(estl,file=here('../tmpdata/estl.Rdata'))
