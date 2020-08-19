## making interpolated lifetables
library(countrycode)
library(glue)
library(akima)
library(data.table)
library(ggplot2)
library(here)

Here <- function(x) gsub("([^\\/\\.]+\\/\\.\\.\\/)","",here(x))


## The repo 'post' should be contained within another folder.
## The here package will identify a git repo as 'here'.
## Various other folders will be created
## at this top level (ie above 'here') to avoid storing data in the repo.
## 'indata' should be unzipped and placed at this top level also.
## 
## make directory structure
if(!file.exists(here('../tmpdata'))) dir.create(here('../tmpdata'))
if(!file.exists(here('../plots'))) dir.create(here('../plots')) #NOTE subfolders
if(!file.exists(here('texto'))) dir.create(here('texto'))
if(!file.exists(here('../figdat'))) dir.create(here('../figdat'))

## subfolders for plots
if(!file.exists(here('../plots/txo'))) dir.create(here('../plots/txo'))
if(!file.exists(here('../plots/hiv'))) dir.create(here('../plots/hiv'))
if(!file.exists(here('../plots/notes'))) dir.create(here('../plots/notes'))
if(!file.exists(here('../plots/inc'))) dir.create(here('../plots/inc'))


## lifetable data
load(here('../indata/LT.Rdata'))                  #data from UN WPP2019
load(here('../indata/H.Rdata'))           #IHME HIV data

LT[,iso3:=countrycode::countrycode(LocID,origin='un',destination='iso3c')]

nrow(LT)
nrow(LT[!is.na(iso3)])

LT <- LT[!is.na(iso3)]
LT[,.(iso3,Time,MidPeriod,Sex,AgeGrp,AgeGrpStart,AgeGrpSpan,lx)]
names(LT)


TT <- LT[,.(iso3,Time,MidPeriod,Sex,AgeGrp,qx,mx)]


agps <- c(0,1,seq(from=5,to=100,by=5))
ags <- c(0,
         paste(agps[2:(length(agps)-1)],
               agps[3:(length(agps))]-1,
               sep='-'),
         '100+')

TT$AgeGrp <- factor(TT$AgeGrp,levels=ags,ordered=TRUE)

tmp <- TT[iso3=='AFG' & MidPeriod==1953]

ggplot(tmp,aes(AgeGrp,mx,col=Sex,group=Sex)) + geom_line()


## testing consistency, recognozing this doesn't capture time trends
tmp[,qxe:=1-exp(-5*mx)]
tmp[AgeGrp=='0',qxe:=1-exp(-1*mx)]
tmp[AgeGrp=='1-4',qxe:=1-exp(-4*mx)]

GP <- ggplot(tmp,aes(AgeGrp,qx,col=Sex,group=Sex)) + geom_line() +
  geom_line(aes(AgeGrp,qxe,col=Sex,group=Sex),lty=2) +
  ggtitle('Dashed line use of mx as hazard, AFG 1953')
GP
ggsave(GP,filename = here('../plots/LTtest1.pdf'))

GP + scale_y_log10()
ggsave(GP,filename = here('../plots/LTtest1log.pdf'))

## inspect
TT[iso3=='AFG' & Sex=='Male' & MidPeriod==1953]

length(ags)
length(agps)

magps <- (c(agps[-1]) + c(rev(rev(agps)[-1])))/2
magps <- c(magps,100)

cbind(agps,magps)                       #midpoints

length(magps)

TT <- TT[order(iso3,Sex,MidPeriod,AgeGrp)]

## which countries are we talking about
allcns <- TT[,unique(iso3)]
length(allcns)
cat(allcns,file=here('texto/allcnsLT.txt'))

## smoothing experiment
mxf <- splinefun(x=agps,
                 y=TT[iso3=='AFG' & Sex=='Female' & MidPeriod==1953,mx])
qplot(agps,mxf(agps),geom='line')


tmp <- TT[iso3=='AFG' & Sex=='Male']

## choose an epoch:
## tmp[MidPeriod==1953]
## tmp[MidPeriod==1958]
tmp[MidPeriod==1963]

ymps <- tmp[,unique(MidPeriod)]
yfls <- min(ymps):max(ymps)

ZM <- matrix(ncol=length(ymps),nrow=length(magps),
             data=tmp[,mx])

## bilinear smoothing example
BL <- bilinear(x=magps,y=ymps,z=ZM,
               x0=1:10,y0=c(1961:1970))

## takes: matrix of age X years
## and a start year and end year
## and age at start year
## returns: someone's (diagonal) yearly hazard of death
getHzhist <- function(yr0,age0,yr1,M){
  ans <- bilinear(x=magps,y=ymps,z=M,
                  x0=age0 + 0:(yr1-yr0),y0=c(yr0:yr1))
  ans <- ans$z
  end <- -Inf
  if(!all(ans[-1]>0)){
    end <- which(ans[-1]>0)
    if(length(end)>0)
      end <- max(end)+1           #last-one-carried fwd for old
    else
      end <- -Inf
  }
  if(end > -Inf & end<length(ans))      #-Inf is if all >0
    ans[end:length(ans)] <- ans[end]
  ans
}

## test
getHzhist(1960,0,2020,ZM)
getHzhist(1960,90,2020,ZM)

ans <- getHzhist(1960,0,2020,ZM)
qplot(1:length(ans),exp(-cumsum(ans)),geom='line')
qplot(1:length(ans),ans,geom='line')
## look OK

## === separating out HIV effects

## find min year across age groups
oldies <- ags[18:length(ags)]


H$age_name <- factor(H$age_name,levels=unique(H$age_name),ordered=TRUE)
HT <- H[iso3=='GUY']                    #example

GP <- ggplot(HT,aes(year,1e2*val,col=sex_name,group=sex_name)) +
  geom_line(linetype=1) +
  geom_ribbon(aes(year,ymin=1e2*lower,ymax=1e2*upper,fill=sex_name),col=NA,alpha=.2) + 
  ylab('HIV prevalence (%)')+
  facet_wrap(~age_name,scales='free_y')
GP

H[val>5e-2,unique(iso3)]

H[,maxv:=max(val),by=iso3]

H[maxv>5e-2,unique((iso3))]
(hivcountries <- H[maxv>5e-2,unique(as.character(iso3))])
hivcountries <- sort(hivcountries)
cat(hivcountries,file=here('texto/hivcountries.txt'))

oldies <- ags[18:length(ags)]

## begin hiv country loop
plt <- TRUE                            #bother saving plots or not?
if(plt) file.remove(file.path(here('../plots/hiv'), list.files(here('../plots/hiv')))) #clear previous

HD <- HEL <- list()
## plt <- FALSE                            #bother saving plots or not?
for(cn in hivcountries){
  print(cn)
  ## pick country data
  tmp <- TT[iso3==cn & Sex != 'Total']

  ## find bump
  bump <- tmp[!AgeGrp %in% oldies & MidPeriod > 1975,{
    ii <- which.max(diff(mx)>0);
    ans <- ifelse(length(ii)>0 & ii < length(MidPeriod),MidPeriod[ii],max(MidPeriod));
    list(yrup=ans)
  },by=.(Sex,AgeGrp)]
  tmp <- merge(tmp,bump,by=c('Sex','AgeGrp'),all.x=TRUE)
  tmp[!is.na(yrup),yrdo:=2050]
  tmp[,ymax:=max(mx),by=.(Sex,AgeGrp)]

  ## approximate w/o bump
  tmpx <- tmp[!is.na(yrup),{
    who <- (MidPeriod<=yrup | MidPeriod>=yrdo);
    sf <- approxfun(x=MidPeriod[who],y=log(mx[who]))
    list(mx0=exp(sf(MidPeriod)),MidPeriod=MidPeriod)
  },
  by=.(AgeGrp,Sex)]

  ## join back in
  tmpb <- merge(tmp,tmpx,by=c('AgeGrp','Sex','MidPeriod'),all.x=TRUE)
  tmpb[!is.na(mx0),mx0:=pmin(mx0,mx)]   #doesn't make things worse

  HD[[cn]] <- tmpb                      #save LT without HIV effects

  ## plots for checking
  ## LT plot with bump smoothed
  plttitle <- glue(cn) + ', age-specific mortality by calendar time'
  plfn <- glue(here('../plots/hiv'))+'/check_'+cn+'.pdf'
  GP <- ggplot(tmpb,aes(MidPeriod,mx,col=Sex,group=Sex)) +
    geom_rect(aes(xmin=yrup,xmax=yrdo,ymin=0,ymax=ymax),alpha=0.1,fill='grey',col=NA)+
    geom_line(aes(y=mx0),linetype=2) +
    geom_line() +
    facet_wrap(~AgeGrp,scales='free_y') +
    ggtitle(plttitle)
  if(plt) ggsave(GP,file=plfn)

  ## HIV extend with cubic splines
  HT <- H[iso3==cn]
  HE <- as.data.table(expand.grid(age_name=HT[,unique(age_name)],
                                  year=c(1975:1989,2018:2020),
                                  sex_name=c('Male','Female'))) #expanded range
  for(sx in c('Male','Female'))
    for(an in HE[,unique(age_name)]){
      xz <- c(1975,HT[age_name==an & sex_name==sx,year])
      yz <- c(-10,HT[age_name==an & sex_name==sx,log(val)])
      sf <- splinefun(xz,yz)
      HE[age_name==an & sex_name==sx,val:=exp(sf(year))]
    }
  HE[,iso3:=cn]
  HEL[[cn]] <- HE

  ## HIVplot
  plttitle <- glue(cn)
  plfn <- glue(here('../plots/hiv'))+'/hiv_'+cn+'.pdf'
  GP <- ggplot(HT,aes(year,1e2*val,col=sex_name,group=sex_name)) +
    geom_line(linetype=1) +
    geom_line(data=HE[year<2000],linetype=2) +
    geom_line(data=HE[year>2000],linetype=2) +
    geom_ribbon(aes(year,ymin=1e2*lower,ymax=1e2*upper,fill=sex_name),
                col=NA,alpha=.2) +
    ylab('HIV prevalence (%)')+
    facet_wrap(~age_name,scales='free_y')+
    ggtitle(plttitle)

  if(plt) ggsave(GP,file=plfn)
}
## end HIV country loop

HD <- rbindlist(HD)
HD[is.na(mx0),mx0:=mx]
HEL <- rbindlist(HEL)                     #extra HIV data from extrap




## -- use the earliest/latest proportional error to propagate uncertainty
## calculate proportional SD
HH <- H[iso3 %in% hivcountries,.(iso3,sex_name,age_name,age_id,year,val,upper,lower)]
HH[,range(year)]
HH[year==1990,pus:=(upper-lower)/(2*val)]
HH[year==2017,puf:=(upper-lower)/(2*val)]
HH

## add to extra data
HEL <- merge(HEL,HH[year==1990,.(iso3,age_name,sex_name,pus)],
              by=c('iso3','age_name','sex_name'))
HEL <- merge(HEL,HH[year==2017,.(iso3,age_name,sex_name,puf)],
              by=c('iso3','age_name','sex_name'))

## calculate the bounds
HEL[year<1990,upper:=val*(1+pus)]; HEL[year<1990,lower:=val*(1-pus)]
HEL[year>2000,upper:=val*(1+puf)]; HEL[year>2000,lower:=val*(1-puf)]
HEL[lower<0,lower:=0]
HEL[,c('puf','pus'):=NULL]

## merge back into all HIV data
HH[,c('age_id','pus','puf'):=NULL]
HH <- rbind(HH,HEL)
HH <- HH[order(iso3,year,sex_name,age_name)]
HH

## check uncertainty procedure for 6 countries
unccheck <- tail(HH[,unique(iso3)])
for(cn in unccheck){
  HT <- HH[iso3==cn]
  plfn <- glue(here('../plots/hiv'))+'/uhiv_'+cn+'.pdf'
  GP <- ggplot(HT,aes(year,1e2*val,col=sex_name,group=sex_name)) +
    geom_line(linetype=1) +
    geom_ribbon(aes(year,ymin=1e2*lower,ymax=1e2*upper,fill=sex_name),
                col=NA,alpha=.2) +
    ylab('HIV prevalence (%)')+
    facet_wrap(~age_name,scales='free_y')+
    ggtitle(HT[,unique(iso3)])
    if(plt) ggsave(GP,file=plfn)
}

## add back to HD
save(HH,file=here('../tmpdata/HH.Rdata'))          #final HIV data by age/sex & year (all iso3)

HD[,c('qx','yrup','yrdo','ymax'):=NULL]
HD


## separate out HIV
## mx = (1-h) * mx0 + h * mxh
HD[,.(iso3,AgeGrp,Sex,MidPeriod,mx,mx0)]
HH[,.(iso3,age_name,Sex=sex_name,MidPeriod=year,h=val)]

## identify age groups across
(anmz <- as.character(HH[,unique(age_name)]))
(agpz <- as.character(HD[,unique(AgeGrp)]))

HD[AgeGrp==agpz[1],age_name:=anmz[1]]
for(i in 1:17)         #print(c(agpz[i+1],anmz[i]))
  HD[AgeGrp==agpz[i+1],age_name:=anmz[i]]
for(i in 18:21)        #print(c(agpz[i+1],anmz[17]))
  HD[AgeGrp==agpz[i+1],age_name:=anmz[17]]

HD <- merge(HD,
            HH[,.(iso3,age_name,Sex=sex_name,MidPeriod=year,h=val)],
            by=c('iso3','MidPeriod','Sex','age_name'),
            all.x=TRUE,all.y=FALSE)

HD[is.na(h),h:=0]
HD[,mxh := (mx - (1-h) * mx0) / h ]
HD[!is.finite(mxh),mxh := 0]
HD[h<1e-5,mxh := 0]                     #dodgy when very small HIV prevalence

HD[mxh>1e3]

HD[,qplot(mxh)] + scale_x_log10()
HD[,summary(mxh)]


HD$AgeGrp <- factor(HD$AgeGrp,levels=ags,ordered = TRUE)
HD <- HD[order(iso3,Sex,MidPeriod,AgeGrp)]

## checks
tst <- HD[iso3=='ZWE' & Sex=='Male' & AgeGrp=='25-29']
plot(tst$mx0,type='b')
plot(tst$mx,col='blue',type='l')
lines(tst$mx0,col=1)
plot(tst$mxh,type='b')

## === end HIV stuff


## === interpolation over years, ages
## NOTE include countries, sex below
PY <- expand.grid(age=0:100,year=1953:2020)
PY <- as.data.table(PY)
PY

## ## add iso3 to IHME data
## H[,iso3:=countrycode::countrycode(location_name,
##                                   origin='country.name',destination='iso3c')]

## single country example
HD[iso3=='ZWE']
HD[iso3=='AFG' & Sex=='Male']

tmp <- TT[iso3=='AFG' & Sex=='Male']
ymps <- tmp[,unique(MidPeriod)]         #1953:2098
yfls <- min(ymps):max(ymps)             #fine time points

ZM <- matrix(ncol=length(ymps),nrow=length(magps),
             data=TT[iso3=='ZAF' & Sex=='Female',mx])
 
ZM0 <- matrix(ncol=length(ymps),nrow=length(magps),
             data=HD[iso3=='ZAF' & Sex=='Female',mx0])

ZMh <- matrix(ncol=length(ymps),nrow=length(magps),
              data=HD[iso3=='ZAF' & Sex=='Female',mxh])

HD[1-mx0/mx > 5e-2,range(MidPeriod)]


H[iso3=='ZAF' & sex_name=='Female']
H[,unique(age_name)]        #different names
H[,length(unique(age_name))]            #just duplicate extras for now

## IHME data
htmp <- H[iso3=='ZAF' & sex_name=='Female']
htmp <- htmp[order(year,age_name)]
htmp[,mdpt:=(year %/% 5)*5 + 3]         #mid-point
htmps <- htmp[,.(val=mean(val)),by=.(mdpt,age_name)]
matrix(ncol=length(1990:2017),nrow=17,
       data=htmp[,age_name])

## tests
getHzhist(1960,0,2020,ZM)
getHzhist(1960,90,2020,ZM)

getHzhist(1960,20,2020,ZM)
getHzhist(1960,20,2020,ZM0)
getHzhist(1960,20,2020,ZMh)

## HR Romanowski
HR <- 2.91
HR.lo <- 2.21
HR.hi <- 3.84
HR.sd <- (HR.hi-HR.lo)/3.92

## get mx for relevant ages and years (diagonal) - single country example
PY[,{
  tmp <- getHzhist(year,age,2020,ZM)
  H <- sum(tmp)
  h <- cumsum(tmp)
  S.m <- exp(-HR*H - (HR.sd*H)^2/2)     #log normal formulae
  S.sd <- S.m * sqrt(exp((H*HR.sd)^2)-1)
  lyv.m <- exp(-HR*h - (HR.sd*h)^2/2)     #log normal formulae
  lyv.sd <- lyv.m * sqrt(exp((h*HR.sd)^2)-1)
  LY.m <- sum(lyv.m)
  LY.sd <- sqrt(sum(lyv.sd^2))
  list(S=S.m,LY=LY.m,S.sd=S.sd,LY.sd=LY.sd)
},by=.(age,year)]


## loop to compute survival tables across all countries
hivcs <- HD[,unique(iso3)]
LA <- list()
cnz <- TT[,unique(iso3)]
k <- 0
for(cn in cnz){
  print(cn)
  for(sx in c('Male','Female','Total')){
    k <- k+1
    ZM <- matrix(ncol=length(ymps),nrow=length(magps),
                 data=TT[iso3==cn & Sex==sx,mx])
    flg <- cn %in% hivcs
    if(flg){
      ZM0 <- matrix(ncol=length(ymps),nrow=length(magps),
                    data=HD[iso3=='ZAF' & Sex=='Female',mx0])
      ZMh <- matrix(ncol=length(ymps),nrow=length(magps),
                    data=HD[iso3=='ZAF' & Sex=='Female',mxh])
    }
    tdf <- PY[,{
      tmp <- getHzhist(year,age,2020,ZM)
      if(flg){
        tmp0 <- getHzhist(year,age,2020,ZM0)
        tmph <- getHzhist(year,age,2020,ZMh)
      } else {
        tmp0 <- tmp
        tmph <- tmp
        tmph <- 0*tmp
      }
      ## total
      H <- sum(tmp)
      h <- cumsum(tmp)
      S.m <- exp(-HR*H - (HR.sd*H)^2/2)     #log normal formulae
      S.sd <- S.m * sqrt(exp((H*HR.sd)^2)-1)
      lyv.m <- exp(-HR*h - (HR.sd*h)^2/2)     #log normal formulae
      lyv.sd <- lyv.m * sqrt(exp((h*HR.sd)^2)-1)
      LY.m <- sum(lyv.m)
      ## LY.sd <- sqrt(sum(lyv.sd^2))      #zero correlation
      LY.sd <- (sum(lyv.sd))            #perfect correlation
      ## HIV-ve
      H.0 <- sum(tmp0)
      h.0 <- cumsum(tmp0)
      S.0.m <- exp(-HR*H.0 - (HR.sd*H.0)^2/2)     #log normal formulae
      S.0.sd <- S.0.m * sqrt(exp((H.0*HR.sd)^2)-1)
      lyv.0.m <- exp(-HR*h.0 - (HR.sd*h.0)^2/2)     #log normal formulae
      lyv.0.sd <- lyv.0.m * sqrt(exp((h.0*HR.sd)^2)-1)
      LY.0.m <- sum(lyv.0.m)
      ## LY.0.sd <- sqrt(sum(lyv.0.sd^2))  #zero correlation
      LY.0.sd <- (sum(lyv.0.sd))  #perfect correlation
      ## HIV+ve
      H.h <- sum(tmph)
      h.h <- cumsum(tmph)
      S.h.m <- exp(-HR*H.h - (HR.sd*H.h)^2/2)     #log normal formulae
      S.h.sd <- S.h.m * sqrt(exp((H.h*HR.sd)^2)-1)
      lyv.h.m <- exp(-HR*h.h - (HR.sd*h.h)^2/2)     #log normal formulae
      lyv.h.sd <- lyv.h.m * sqrt(exp((h.h*HR.sd)^2)-1)
      LY.h.m <- sum(lyv.h.m)
      LY.h.sd <- (sum(lyv.h.sd))  #perfect correlation
      list(S=S.m,LY=LY.m,S.sd=S.sd,LY.sd=LY.sd,
           S.0=S.0.m,LY.0=LY.0.m,S.0.sd=S.0.sd,LY.0.sd=LY.0.sd,
           S.h=S.h.m,LY.h=LY.h.m,S.h.sd=S.h.sd,LY.h.sd=LY.h.sd)
    },by=.(age,year)]
    tdf[,c('iso3','sex'):=list(cn,sx)]
    LA[[k]] <- tdf
  }
}
LA <- rbindlist(LA)
LA <- LA[age<100]

LA[!iso3 %in% hivcs,c('S.h','LY.h','S.h.sd','LY.h.sd'):=0]


## safety
LA[,S.h:=min(S.h,S.0),by=.(iso3,year,sex,age)]

## save to file
save(LA,file=here('../tmpdata/LA.Rdata'))                #

## === checks
print(LA[iso3=='AFG' & year==2000 & sex=='Male'],n=Inf)


## NOTE be careful about 100 for LYs
tmp <- LA[iso3=='GBR' & year==2015 & sex=='Male']
ggplot(tmp,aes(age,S))  + geom_line()
ggplot(tmp,aes(age,LY))  + geom_line()


## face validity in ZAF
tmp <- LA[iso3=='ZAF' & year==2015 & sex=='Male']
ggplot(tmp,aes(age,S))  + geom_line()
ggplot(tmp,aes(age,S.0))  + geom_line()
ggplot(tmp,aes(age,S.h))  + geom_line()


ggplot(tmp) + geom_line(aes(age,S.0)) + geom_line(aes(age,S),col='blue') + geom_line(aes(age,S.h),col='red') + xlab('Age in 2015') + ylab('Probability of survival to 2020')
ggplot(tmp) + geom_line(aes(age,LY.0)) + geom_line(aes(age,LY),col='blue') + geom_line(aes(age,LY.h),col='red') + xlab('Age in 2015') + ylab('Expected LYs to 2020')
