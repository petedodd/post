## making the age-specific CDRs

## TB notification imputations
## TB treatment outcomes
rm(list=ls())
library(here)
source(here('Neat/0_utilities.R'))

if(! overwrite ){
  fn <- here('tmpdata/RRcdr.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/RRcdr.Rdata exists!')
}

## load data
load(here('tmpdata/TBM.Rdata'))
est <- fread(here('indata/TB_burden_countries_2020-02-24.csv'))
## age/sex splits
AA <- fread(here("indata/TB_burden_age_sex_2020-02-24.csv"))

## calculating
AA <- AA[age_group %in% c('0-4','5-14'),.(iso3,acat=age_group,sex,best)]
AA2 <- TBM[acat %in% c('0-4','5-14') & year==2017,.(iso3,acat,Sex,value)]
AA[,Sex:='Male']
AA[sex=='f',Sex:='Female']
AA <- merge(AA,AA2,by=c('iso3','Sex','acat'),all.x=TRUE,all.y=FALSE)
AA <- AA[,.(value=sum(value,na.rm=TRUE),best=sum(best)),by=.(iso3,acat)]
AA[,cdr:=value/(best + 1e-9)]           #child CDRs
## AA <- AA[,.(cdr=sum(best*cdr)/sum(best+1e-9)),by=.(iso3,acat)] #weighted average
AA[cdr>1.1,cdr:=1]                      #safety
AA <- merge(AA,unique(TBM[,.(iso3,g_whoregion)]))
AA[,regm:=weighted.mean(x=cdr,w=best,na.rm=TRUE),by=g_whoregion]
AA[is.na(cdr)|cdr==0,cdr:=regm]
AA2 <- est[year==2017,.(iso3,ocdr=c_cdr/1e2,ocdr.sd=(c_cdr_hi-c_cdr_lo)/392)]
AA2[is.na(ocdr),c('ocdr','ocdr.sd'):=.(1,0)]
AA <- merge(AA,AA2,by='iso3',all.x=TRUE,all.y=TRUE)


## not much evidence for sex difference
## main age difference in AFR
## 
## ocdr = (n5+n15+na) / (i5+i15+ia)
##  = (n5+n15+na) / (n5/cdr5+n15/cdr15+na/cdra)
## onotes / ocdr = (n5/cdr5+n15/cdr15+na/cdra)
## 1 / ocdr = (p5/cdr5+p15/cdr15+pa/cdra)
## 1 / ocdr - (p5/cdr5+p15/cdr15)= pa/cdra
## pa / (1 / ocdr - (p5/cdr5+p15/cdr15))= cdra #<- use this
##
## 1 / ocdr = (p5*(cdra/cdr5) + p15*(cdra/cdr15) + pa) / cdra
## cdra = (p5*RR5 + p15*RR14 + pa) * ocdr

## calculating RRs for children
AA2 <- TBM[year==2017,.(iso3,acat,Sex,value)]
AA2[,tot:=sum(value),by=iso3]
AA2[acat=='0-4',u5:=value]
AA2[acat=='5-14',b515:=value]
AA2[,u5:=sum(u5,na.rm=TRUE),by=iso3]
AA2[,b515:=sum(b515,na.rm=TRUE),by=iso3]
AA2 <- unique(AA2[,.(iso3,tot,u5,b515)])
AA2[,c('pa','p5','p15'):=.((tot-u5-b515)/(tot+1e-9),
(u5)/(tot+1e-9),
(b515)/(tot+1e-9))]
AAM <- dcast(AA[,.(iso3,acat,cdr)],iso3 ~ acat,value.var = 'cdr')
AAM <- merge(AAM,unique(AA[,.(iso3,ocdr)])) #3 CDRs
AA2 <- merge(AA2,AAM,by='iso3')
AA2[,cdra:=pa / (1 / ocdr - (p5*(1/`0-4`)+p15*(1/`5-14`)))] #see above
AA2[,RRu5:=cdra/`0-4`]
AA2[,RR514:=cdra/`5-14`]

## make RR-CDR
RRcdr <- AA2[,.(iso3,RRu5,RR514)]
## safety -- mostly close to 1
RRcdr[RRu5<1,RRu5:=1]
RRcdr[RR514<1,RR514:=1]

## averages for missed or odd
(missed <- setdiff(AA$iso3,RRcdr$iso3))
cat(missed,file=here::here('post/texto/RRnone.txt'))
missed <- unique(AA[iso3 %in% missed,.(iso3,g_whoregion)])
RRcdr <- merge(RRcdr,unique(AA[,.(iso3,g_whoregion)]),
               by='iso3',all.x=TRUE,all.y=FALSE)
summary(RRcdr)
RRcdrR <- RRcdr[,.(RRu5a=median(RRu5),RR514a=median(RR514)),by=g_whoregion]
RRcdr <- rbind(RRcdr,missed,fill=TRUE)

RRcdr <- merge(RRcdr,RRcdrR,by='g_whoregion',all.x=TRUE)
RRcdr[is.na(RRu5),c('RRu5','RR514'):=.(RRu5a,RR514a)] #missing set to regional averages
## outliers
RRcdr[RRu5>quantile(RRu5,.975)]
RRcdr[RR514>quantile(RR514,.975)]

(missed <- RRcdr[RRu5>quantile(RRu5,.975),unique(iso3)])
cat(missed,file=here::here('post/texto/RRhiu5.txt'))
(missed <- RRcdr[RR514>quantile(RR514,.975),unique(iso3)])
cat(missed,file=here::here('post/texto/RRhi514.txt'))
RRcdr[RRu5>quantile(RRu5,.975),RRu5:=RRu5a]
RRcdr[RR514>quantile(RR514,.975),RR514:=RR514a]

RRcdr[,c('g_whoregion','RRu5a','RR514a'):=NULL]

summary(RRcdr)

save(RRcdr,file=here::here('tmpdata/RRcdr.Rdata'))
