## this will be the tables
rm(list=ls())
library(here)
library(data.table)

whozg <- c('AFR','AMR','EMR','EUR','SEA','WPR')
whozt <- c('Africa','The Americas','Eastern Mediterranean','Europe','South-East Asia',
           'Western Pacific')
wrk <- data.table(g_whoregion=whozg,name=whozt)

ftb <- Vectorize(function(x){
  ## formatter according to GTB rounding rules
  stopifnot (!is.character(x))
  if (!is.na(x)){
    smallpos <- x>0 & x<0.01
    dg <- ifelse(abs(x)>0.01 & abs(x)<100, 2, 3)
    x <- signif(x, 3)
    x <- format(x, digits=dg, nsmall=0, big.mark=" ", justify='right', drop0trailing = T)
    if (smallpos) x <- '<0.01'
  } else x <- '-'
  return(x)}, 'x')

## looking at big numbers
see <- function(x,ns=3)formatC(signif(x,ns),big.mark = ",",format='fg') #for reading big numbers

Ssum <- function(x,...) sqrt(sum(x^2,...))

## Function for rounding anything left with a decimal point and removing '<'
## (for number formatting for tables)
zapdotlt <- function(x){
  x <- gsub("<","",x)
  got <- grep("\\.",x)
  x[got] <- as.character(round(as.numeric(x[got])))
  x
}

## Function for formatting big numbers >=1e6, which are breaking with ftb
fmtbig <- function(x,sf=3){
  x <- signif(x, sf)
  x <- format(x, digits=3,
              nsmall=0, big.mark=" ", justify='right',
              drop0trailing = TRUE,scientific=FALSE)
  ## x <- stringi::stri_enc_toutf8(x)      #ensure UTF8
  x
}
fmtbig(123451000)
fmtbig(123451000,sf=4)

sda <- function(a.sd,b.sd) sqrt((a.sd)^2+(b.sd)^2)
sdq <- function(a,a.sd,b,b.sd) (a/b) * sda(a.sd/a,b.sd/b)

pcamong <- function(nmr,den,sf=4){
  out <- c(nmr[1]/den[1],sdq(nmr[1],nmr[2],den[1],den[2]))
  out <- 1e2*out
  c(fmtbig(out[1],sf=sf),
    fmtbig(out[1]-out[2]*1.96,sf=sf),
    fmtbig(out[1]+out[2]*1.96,sf=sf))
}


## ============ Table 1
load(here('../figdat/t1r1.Rdata'))         #Total new tuberculosis cases 1980-2019
load(here('../figdat/t1r2.Rdata'))         #New tuberculosis cases treated 1980-2019
load(here('../figdat/t1r3.Rdata'))         #New tuberculosis cases not treated 1980-2019
## load(here('figdat/t1r4.Rdata')) # = 5 + 6 #Total tuberculosis survivors alive in 2020
load(here('../figdat/t1r5.Rdata'))         #Treated tuberculosis survivors alive in 2020
load(here('../figdat/t1r6.Rdata'))       #Untreated tuberculosis survivors alive in 2020
## load(here('figdat/t1r7.Rdata')) # = 8 + 9 #Life-years post-tuberculosis
load(here('../figdat/t1r8.Rdata'))         #Life-years post-tuberculosis lived by treated
load(here('../figdat/t1r9.Rdata'))       #Life-years post-tuberculosis lived by untreated

t1r4 <- merge(t1r5[,.(g_whoregion,value,value.sd,quantity)],
              t1r6[,.(g_whoregion,value,value.sd,quantity)],by='g_whoregion')
t1r7 <- merge(t1r8,t1r9,by='g_whoregion')
t1r4[,value:=value.x+value.y]; t1r7[,value:=value.x+value.y]
t1r4[,value.sd:=sqrt(value.sd.x^2+value.sd.y^2)];
t1r7[,value.sd:=sqrt(value.sd.x^2+value.sd.y^2)]
t1r4[,quantity:='tot2020']; t1r7[,quantity:='LY2020'];
t1r4[,c('quantity.x','quantity.y'):=NULL];t1r7[,c('quantity.x','quantity.y'):=NULL];
t1r4[,c('value.x','value.y'):=NULL];t1r4[,c('value.sd.x','value.sd.y'):=NULL];
t1r7[,c('value.x','value.y'):=NULL];t1r7[,c('value.sd.x','value.sd.y'):=NULL];

## add row numbers
t1r1[,rn:=1]; t1r2[,rn:=2]; t1r3[,rn:=3]
t1r4[,rn:=4]; t1r5[,rn:=5]; t1r6[,rn:=6]
t1r7[,rn:=7]; t1r8[,rn:=8]; t1r9[,rn:=9]

## SAs: 4,5,6

t1 <- rbindlist(list(t1r1,t1r2,t1r3,t1r4,
                     t1r5[,.(g_whoregion,value,value.sd,quantity,rn)],
                     t1r6[,.(g_whoregion,value,value.sd,quantity,rn)],
                     t1r7,t1r8,t1r9),use.names = TRUE)

namekey <- data.table(rn=1:9,
                      nm=c('Total new tuberculosis cases 1980-2019',
                           'New tuberculosis cases treated 1980-2019',
                           'New tuberculosis cases not treated 1980-2019',
                           'Total tuberculosis survivors alive in 2020',
                           'Treated tuberculosis survivors alive in 2020',
                           'Untreated tuberculosis survivors alive in 2020',
                           'Life-years lived by tuberculosis survivors 1980-2020',
                         'Life-years lived by treated tuberculosis survivors 1980-2020',
                         'Life-years lived by untreated tuberculosis survivors 1980-2020'
                           ))

t1[,.N,by=quantity]

if(! 'Global' %in% wrk$g_whoregion) wrk <- rbind(wrk,data.table(g_whoregion='Global',name='Global'))
t1 <- merge(t1,namekey,by='rn')
t1 <- merge(t1,wrk,by='g_whoregion')

t1$name <- factor(t1$name,levels=wrk$name,ordered=TRUE)
t1$nm <- factor(t1$nm,levels=namekey$nm,ordered=TRUE)


T1 <- dcast(t1,nm ~ name,value.var = 'value')
T1


tmp <- T1[,lapply(.SD,fmtbig),.SDcols=2:ncol(T1)]
T1F <- cbind(quantity=T1[,nm],tmp)

fwrite(T1F,file=here('figs/table1.csv'))

## newversion with unc
T1.sd <- dcast(t1,nm ~ name,value.var = 'value.sd')
T1.sd

nmz <- names(T1)[2:ncol(T1)]

T1.sd <- as.matrix(T1.sd[,..nmz])
T1.mid <- as.matrix(T1[,..nmz])

T1.lo <- T1.mid - T1.sd*1.96
T1.lo[T1.lo<0] <- 0
T1.hi <- T1.mid + T1.sd*1.96

SF <- 3
T1.lo <- fmtbig(T1.lo,SF)
T1.hi <- fmtbig(T1.hi,SF)

BB <- matrix(paste0("(",T1.lo," - ",T1.hi,")"),nrow=nrow(T1.lo),ncol=ncol(T1.lo))
BB2 <- fmtbig(T1.mid,SF)
BB2 <- matrix(paste0(BB2,"  ",BB),nrow=nrow(T1.lo),ncol=ncol(T1.lo))

rownames(BB2) <- T1$nm
colnames(BB2) <- names(T1)[-1]
BB2
print(BB2[,7])

## percentage unc int
pcunc <- 392*T1.sd[,ncol(T1.sd)]/T1.mid[,ncol(T1.mid)]
print(data.frame(thing=T1$nm,pcunc=see(pcunc,ns=2)))


write.csv(BB2,file=here('figs/table1u.csv'))

## country level table
## load(here('../figdat/t1r1c.Rdata'))        #Total new tuberculosis cases 1980-2019
load(here('../figdat/t1r2c.Rdata'))       #New tuberculosis cases treated 1980-2019
load(here('../figdat/t1r3c.Rdata'))       #New tuberculosis cases not treated 1980-2019
## load(here('figdat/t1r4.Rdata')) # = 5 + 6 #Total tuberculosis survivors alive in 2020
load(here('../figdat/t1r5c.Rdata'))       #Treated tuberculosis survivors alive in 2020
load(here('../figdat/t1r6c.Rdata'))       #Untreated tuberculosis survivors alive in 2020
## load(here('figdat/t1r7.Rdata')) # = 8 + 9 #Life-years post-tuberculosis
load(here('../figdat/t1r8c.Rdata'))       #Life-years post-tuberculosis lived by treated
load(here('../figdat/t1r9c.Rdata'))       #Life-years post-tuberculosis lived by untreated

T1C <- rbindlist(list(t1r2c,t1r3c,t1r5c,t1r6c,t1r8c,t1r9c))
T1C[,unique(quantity)]                  #OK
fwrite(T1C,file=here('figs/countryest.csv'))

## ======== Table 2 ===========
## see Results
load(here("../figdat/c1.Rdata")); load(here("../figdat/c1b.Rdata"))
load(here("../figdat/c2.Rdata")); load(here("../figdat/c2b.Rdata"))
load(here("../figdat/c4.Rdata")); load(here("../figdat/c4b.Rdata"))


## more SF for table 2
SF <- 4

## combine
c1[,totf:=fmtbig(total,SF)]
c1[,totf.lo:=fmtbig(total-total.sd*1.96,SF)]; c1[,totf.hi:=fmtbig(total+total.sd*1.96,SF)]
c1[,totf:=paste0(totf," (",totf.lo," - ",totf.hi,")")]
c2[,pcr:=round(pc)]

## c4[,age:=paste0(round(wm),"  (",round(wsd),")")]
c4[,age:=round(pck,1)]

c1b[,totf:=fmtbig(total,SF)]
c1b[,totf.lo:=fmtbig(total-total.sd*1.96,SF)]; c1b[,totf.hi:=fmtbig(total+total.sd*1.96,SF)]
c1b[,totf:=paste0(totf,"  (",totf.lo," - ",totf.hi,")")]
c2b[,pcr:=round(pc)]

## c4b[,age:=paste0(round(wm)," (",round(wsd),")")]
c4b[,age:=round(pck,1)]

tab2 <- merge(c1[,.(g_whoregion,total5=totf)],c2[,.(g_whoregion,pcr5=pcr)],
              by='g_whoregion')
tab2 <- merge(tab2,c4[,.(g_whoregion,age5=age)], by='g_whoregion')
tab2 <- merge(tab2,c1b[,.(g_whoregion,total2=totf)], by='g_whoregion')
tab2 <- merge(tab2,c2b[,.(g_whoregion,pcr2=pcr)], by='g_whoregion')
tab2 <- merge(tab2,c4b[,.(g_whoregion,age2=age)], by='g_whoregion')
setkey(tab2,g_whoregion)
tab2 <- tab2[c('AFR','AMR','EMR','EUR','SEA','WPR','Global')] #reorder

write.csv(tab2,file=here::here('figs/table2.csv'))


## ======= Other stats =========

## ------- additional stats to move

## ------- stats for untx + tx -------
## HIV in new
load(file=here('../figdat/t1r1.Rdata'))
load(file=here('../tmpdata/NHu.Rdata'))

nmr <- c(NHu[,value],NHu[,value.sd])
den <- c(t1r1[g_whoregion=='Global',value],t1r1[g_whoregion=='Global',value.sd])

(out <- pcamong(nmr,den,sf=2))

cat(out,file=here('texto/s_HN.txt'),sep=' - ')

## paediatric LYS
## untreated
load(file=here('../tmpdata/LYpu.Rdata'))
load(file=here('../figdat/t1r9.Rdata'))
## treated
load(file=here('../tmpdata/LYpt.Rdata'))
load(file=here('../figdat/t1r8.Rdata'))

nmr <- c(LYpu[,value] + LYpt[,value],
         sda(LYpu[,value.sd],LYpt[,value.sd]))

den <- c(t1r8[g_whoregion=='Global',value] + t1r9[g_whoregion=='Global',value],
         sda(t1r8[g_whoregion=='Global',value.sd],t1r9[g_whoregion=='Global',value.sd]))


(out <- c(fmtbig(nmr[1],sf=4),
         fmtbig(nmr[1]-nmr[2]*1.96,sf=4),
         fmtbig(nmr[1]+nmr[2]*1.96,sf=4))
)

cat(out,file=here('texto/s_ply.txt'),sep=' - ')

(out <- pcamong(nmr,den))

cat(out,file=here('texto/s_plypc.txt'),sep=' - ')



## SEA of survivors
load(file=here('../figdat/t1r5.Rdata')) #treated survivors
load(file=here('../figdat/t1r6.Rdata')) #untreated survivors


nmr <- c(t1r5[g_whoregion=='SEA',value] + t1r6[g_whoregion=='SEA',value],
         sda(t1r5[g_whoregion=='SEA',value.sd],t1r6[g_whoregion=='SEA',value.sd]))

den <- c(t1r5[g_whoregion=='Global',value] + t1r6[g_whoregion=='Global',value],
         sda(t1r5[g_whoregion=='Global',value.sd],t1r6[g_whoregion=='Global',value.sd]))


(out <- pcamong(nmr,den))

cat(out,file=here('texto/s_SEApc.txt'),sep=' - ')

## male of survivors
load(file=here('../tmpdata/SMu.Rdata')) #untreated
load(file=here('../tmpdata/SMt.Rdata')) #treated


nmr <- c(SMu[,value] + SMt[,value],
         sda(SMu[,value.sd],SMt[,value.sd]))

(out <- pcamong(nmr,den))

cat(out,file=here('texto/s_SMpc.txt'),sep=' - ')


## HIV of survivors
load(file=here('../tmpdata/SHu.Rdata')) #untreated
load(file=here('../tmpdata/SHt.Rdata')) #treated


nmr <- c(SHu[,value]+SHt[,value],sda(SHu[,value.sd],SHu[,value.sd]))

(out <- pcamong(nmr,den,sf=2))

cat(out,file=here('texto/s_SHpc.txt'),sep=' - ')


load(file=here('../figdat/txay.Rdata'))
load(file=here('../figdat/utay.Rdata'))
both <- rbind(txay,utay)

## years lived among survivors
tmp <- both[,.(wt=sum(wt)),by=YPT]

## ggplot(tmp,aes(YPT,wt)) + geom_line()

out <- weighted.mean(tmp$YPT,tmp$wt)    #the mean
rsds <- tmp$YPT-out                     # residuals
ods <- sqrt(weighted.mean(rsds^2,tmp$wt)) #NOTE this doesn't give a meaninful UI for mean
out <- c(fmtbig(out,sf=4))

cat(out,file=here('texto/s_SYL0.txt'),sep=' - ')


long <- sample(1:nrow(tmp),1e4,replace=TRUE,prob=tmp$wt)
tmp <- tmp[long]

out <- c(quantile(tmp$YPT,0.5),
         quantile(tmp$YPT,0.25),
         quantile(tmp$YPT,0.75))

cat(out,file=here('texto/s_SYL.txt'),sep=' - ')


## age of survivors
tmp <- both[,.(wt=sum(wt)),by=agenow]
long <- sample(1:nrow(tmp),1e4,replace=TRUE,prob=tmp$wt)
tmp <- tmp[long]

out <- c(quantile(tmp$agenow,0.5),
         quantile(tmp$agenow,0.25),
         quantile(tmp$agenow,0.75))

cat(out,file=here('texto/s_SA.txt'),sep=' - ')

## ------- stats for 5+2  -------
load(here("../figdat/c1.Rdata")); load(here("../figdat/c1b.Rdata"))

## as % of all alive
nmr5 <- c(c1[g_whoregion=='Global',total],
          c1[g_whoregion=='Global',total.sd])

nmr2 <- c(c1b[g_whoregion=='Global',total],
          c1b[g_whoregion=='Global',total.sd])

(out <- pcamong(nmr5,den,sf=3))

cat(out,file=here('texto/s_5pc.txt'),sep=' - ')

(out <- pcamong(nmr2,den,sf=2))

cat(out,file=here('texto/s_2pc.txt'),sep=' - ')


## HIV among AFR last 5
load(file=here("../figdat/c1h.Rdata")); load(file=here("../figdat/c1hb.Rdata"))

nmr5 <- c(c1h[,total],
          c1h[,total.sd])

nmr2 <- c(c1hb[,total],
          c1hb[,total.sd])

den5 <- c(c1[g_whoregion=='AFR',total], #AFR
          c1[g_whoregion=='AFR',total.sd])

den2 <- c(c1b[g_whoregion=='AFR',total], #AFR
          c1b[g_whoregion=='AFR',total.sd])


(out2 <- pcamong(nmr2,den2,sf=2))
(out5 <- pcamong(nmr5,den5,sf=2))

cat(out5,file=here('texto/s_5Hpc.txt'),sep=' - ')
cat(out2,file=here('texto/s_2Hpc.txt'),sep=' - ')

## Male among
load(file=here("../figdat/c1m.Rdata")); load(file=here("../figdat/c1mb.Rdata"))

nmr5 <- c(c1m[g_whoregion=='Global',total],
          c1m[g_whoregion=='Global',total.sd])

nmr2 <- c(c1mb[g_whoregion=='Global',total],
          c1mb[g_whoregion=='Global',total.sd])

den5 <- c(c1[g_whoregion=='Global',total],
         c1[g_whoregion=='Global',total.sd])

den2 <- c(c1b[g_whoregion=='Global',total],
         c1b[g_whoregion=='Global',total.sd])


(out2 <- pcamong(nmr2,den2,sf=3))
(out5 <- pcamong(nmr5,den5,sf=3))       #OK

cat(out5,file=here('texto/s_5Mpc.txt'),sep=' - ')
cat(out2,file=here('texto/s_2Mpc.txt'),sep=' - ')

## age IQRs
load(file=here("../figdat/cage.Rdata")); load(file=here("../figdat/cageb.Rdata"))

## calculate weighted quantiles...
cage <- cage[,.(wt=sum(alive)),by=agenow]
long <- sample(1:nrow(cage),1e4,replace=TRUE,prob=cage$wt)
cagel <- cage[long]

cageb <- cageb[,.(wt=sum(alive)),by=agenow]
longb <- sample(1:nrow(cageb),1e4,replace=TRUE,prob=cageb$wt)
cagelb <- cageb[long]

out5 <- c(quantile(cagel$agenow,0.5),
         quantile(cagel$agenow,0.25),
         quantile(cagel$agenow,0.75))

out2 <- c(quantile(cagelb$agenow,0.5),
         quantile(cagelb$agenow,0.25),
         quantile(cagelb$agenow,0.75))

cat(out2,file=here('texto/s_2ageQ.txt'),sep=' - ')
cat(out5,file=here('texto/s_5ageQ.txt'),sep=' - ')


## SAs -----
## sensitivity ananlyses around HR
svvs1 <- t1r5[,.(g_whoregion,value,value.2,value.4,value.d)]
svvs2 <- t1r6[,.(g_whoregion,value,value.2,value.4,value.d)]
svvs1[,quantity:='Treated survivors']; svvs2[,quantity:='Untreated survivors']
svvs <- rbind(svvs1,svvs2)
names(svvs)[2:5] <- c('basecase','HR=2','HR=4','dynamic HR')
names(svvs)[1] <- 'Region'


svvo <- svvs[,lapply(.SD,fmtbig),.SDcols=2:5]
svvo <- cbind(Region=svvs$Region,svvo)

fwrite(svvo,file=here('figs/SAtable1.csv'))

## treated surv, untreated surv
## <5, <2
load(here("../figdat/c1sa.Rdata")); load(here("../figdat/c1bsa.Rdata"))

c1o <- c1[,.(g_whoregion,total)]
c1o <- merge(c1o,c1sa[,.(g_whoregion,total.2,total.4,total.d)],by='g_whoregion')
c1o[,quantity:='treated last 5 years']
c2o <- c1b[,.(g_whoregion,total)]
c2o <- merge(c2o,c1bsa[,.(g_whoregion,total.2,total.4,total.d)],by='g_whoregion')
c2o[,quantity:='treated last 2 years']
setkey(c1o,g_whoregion)
setkey(c2o,g_whoregion)

csab <- rbind(c1o[wrk$g_whoregion],c2o[wrk$g_whoregion])
csab <- csab[,.(g_whoregion,quantity,total,total.2,total.4,total.d)]
names(csab) <- c('Region','quantity','basecase','HR=2','HR=4','dynamic HR')

sa2 <- csab[,lapply(.SD,fmtbig),.SDcols=3:6]
sa2 <- cbind(csab[,.(Region,quantity)],sa2)

fwrite(sa2,file=here('figs/SAtable2.csv'))
