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
## N3[,value0:=value]
## n3 <- N3[,.(iso3,year,age,sex,value=value0)] #this is correcting for death
if(! 'acat' %in% names(N3)) N3 <- merge(N3,lamap[,.(age,acat=acats)],by='age') #merge in coarse cats if not already there


N3[,totn:=sum(value),by=.(iso3,year,acat,sex)] #notifications across each coarse cat
N3[,sum(totn)]/1e6

## want to spread the gap proportional to the npc
N3[,totyr:=sum(value)+1e-9,by=.(iso3,year)] #notifications across country years
N3[,propnote:=value/totyr] #proportion of notifications in each age/sex

N3[,sum(propnote),by=.(iso3,year)]      #check
N3[,summary(propnote)]
N3[is.na(propnote)]                     #no cases

estl <- merge(est,N3[,.(iso3,year,sex,age,acat,propnote)],
              by=c('iso3','year'),all.y=TRUE)

estl <- merge(estl,RRcdr,by='iso3',all.x=TRUE)

## u5 fraction
estl[,p5:=NA_real_]
estl[acat=='0-4',p5:=propnote]
estl[,p5:=sum(p5,na.rm=TRUE),by=.(iso3,year)]
## 5-14 fraction
estl[,p15:=NA_real_]
estl[acat=='5-14',p15:=propnote]
estl[,p15:=sum(p15,na.rm=TRUE),by=.(iso3,year)]

## cdra = (p5*RR5 + p15*RR14 + pa) * ocdr
estl[,cdra := (p5*RRu5 + p15*RR514 + (1-p15-p5)) * ocdr]
estl[,cdru5 := cdra/RRu5]
estl[,cdr514 := cdra/RR514]
estl[,tcdr:=cdra]
estl[acat=='0-4',tcdr:=cdru5]
estl[acat=='5-14',tcdr:=cdr514]

## disaggregate gap TODO worry this loses IND correction (see above) jj
estl[,gapl2:=gap*propnote]               #1 years age/sex gap
estl[,gapl:=propnote * c_newinc * (1/ tcdr-1)]               #1 years age/sex gap
estl[gapl<0,gapl:=0]               #safety
## estl[,gapl.sd:=propnote * c_newinc * (ocdr.sd/( tcdr * ocdr))]   #Taylor expansion, perfect correlation

estl[!is.finite(gap),unique(iso3)]      #just SDN TODO
estl <- estl[iso3!='SDN']               #drop

## F = sd/X
## F_i = sd_i/X_i = F / (sqrt(sum_jp_j^2)
estl[,summary(ocdr.sd)]
estl[,summary(ocdr.sd/ocdr)]
estl[(ocdr.sd>ocdr),ocdr.sd:=ocdr]      #safety

estl[,pp:=gapl/sum(gapl),by=.(iso3,year)]
estl[!is.finite(pp),pp:=0]
estl[,den:=sqrt(sum(pp^2)),by=.(iso3,year)]

estl[,gapl.sd :=  gapl * (ocdr.sd/ocdr) / den]
estl[!is.finite(gapl.sd),gapl.sd:=0]


## TODO jj checks
tmp <- estl[,.(v1=sum(gapl.sd^2)/(sum(gapl)^2+1e-9),
               v2=(gap.sd[1]/(gap[1]+1e-9))^2),by=.(iso3,year)]
tmp
tmp[,summary(v1-v2)]
tmp[abs(v1-v2)>.01]

estl[,.(gapl=sum(gapl),gapl2=sum(gapl2)),by=.(iso3,year)] #check


estl[,bad:=FALSE]
estl[!is.finite(propnote),bad:=TRUE,by=.(iso3,year)]
ncts <- length(estl[,unique(age)])
estl[bad==TRUE,gapl:=gap/ncts]          #safety




## gap test
N3[year==2017,sum(value)]/1e6
estl[year==2017,sum(gapl)]/1e6
estl[year==2017,sum(gapl2,na.rm=TRUE)]/1e6 #approximate check

estl[,sum(gapl)/1e6,by=year]

ggplot(estl[,.(gap=sum(gapl,na.rm=TRUE)/1e6),by=year],aes(year,gap)) + geom_line() + scale_y_continuous(label=absspace)

ggplot(estl[,.(gap=sqrt(sum(gapl.sd^2,na.rm=TRUE))/sum(gapl,na.rm=TRUE)),by=year],aes(year,gap)) + geom_line() + scale_y_continuous(label=absspace)
## odd blow up TODO

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
