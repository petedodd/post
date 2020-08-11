## imputing incidence & gaps


rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/estg.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/estg.Rdata exists!')
}

## load data
est <- fread(here("../indata/TB_burden_countries_2020-02-24.csv"))
load(here('../tmpdata/TBN.Rdata'))


## start work
est <- est[,.(iso3,year,e_inc_num,ocdr=c_cdr/1e2,ocdr.sd=(c_cdr_hi-c_cdr_lo)/392)]
xtray <- c(1980:1999)

est <- rbind(est,data.table(iso3=rep(est[,unique(iso3)],each=length(xtray)),
                            year=xtray,e_inc_num=NA,ocdr=NA,ocdr.sd=NA))
est <- est[order(year)]

## NOTE need to extrapolate back to 1980
TBN[year==2000,ratio:=as.numeric(c_newinc),by=iso3]
TBN[,ratio:=mean(ratio,na.rm=TRUE),by=iso3] #same for all years
TBN[,ratio:=c_newinc/ratio,by=iso3]

TBN[,.(iso3,year,c_newinc,ratio)]
TBN[is.finite(ratio),summary(ratio)]
TBN[is.finite(ratio),qplot(ratio)]
TBN[is.finite(ratio) & ratio > 10 & year<2000]
(badfac <- TBN[is.finite(ratio) & ratio > 5 & year<2000,unique(iso3)])

cat(badfac,file=here('texto/badfac.txt'))


TBN[iso3 %in% badfac,ratio:=1] #cap

## TODO
## uncertainty here in projections
## jj


est2000 <- est[year==2000]

TBNx <- merge(est2000[,.(iso3,e_inc_num)],
              TBN[,.(iso3,year,ratio)],by=c('iso3'),all.y=TRUE) #
TBNx[,e_inc_num2:=e_inc_num * ratio]                            #NOTE

est <- merge(est,TBNx[,.(iso3,year,e_inc_num2)],by=c('iso3','year'),all.x=TRUE)
## NOTE change type of inc_num
est[is.na(e_inc_num),e_inc_num:=e_inc_num2]
est[,e_inc_num2:=NULL]

## 2019 separately
est <- rbind(est[,.(year,iso3,e_inc_num,ocdr,ocdr.sd)],
             data.table(year=2019,iso3=unique(est$iso3),
                        e_inc_num=NA,ocdr=NA,ocdr.sd=NA))
est[,e_inc_num:=as.integer(na_interpolation(e_inc_num)),by=iso3]

est[,imp:=FALSE]
est[year %in% c(xtray,2019),imp:=TRUE]
est <- merge(est,isokey,by='iso3')

if(plt){
  for(reg in est[,unique(g_whoregion)]){
    GP <- ggplot(est[g_whoregion==reg],aes(year,e_inc_num,col=imp))  +
      geom_point() + #geom_line() +
      facet_wrap(~iso3,scales='free_y') + theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(GP,filename=here::here(paste0('../plots/inc/incIMPcheck_',reg,'.pdf')),
           w=9,h=7)
  }
}

## TODO there are some issues here
ggplot(est[iso3=='PAK'],aes(year,e_inc_num)) + geom_line()

est[year==2017,sum(e_inc_num)/1e6]
est[,sum(e_inc_num)/1e6,by=.(year)]
dcast(est[,sum(e_inc_num)/1e6,by=.(year,g_whoregion)],year~g_whoregion,value.var = 'V1')
ggplot(est[iso3=='IND'],aes(year,e_inc_num)) + geom_line()

ggplot(est[g_whoregion=='EMR'],aes(year,e_inc_num)) + geom_line() +
  facet_wrap(~iso3,scales='free')


## NOTE can't re-run without fresh TBN
tmpi <- est[,.(incnum=sum(e_inc_num)/1e6),by=.(year)] 

##
GP <- ggplot(tmpi,aes(year,incnum)) +
  geom_line() + ylab('Total estimated TB incidence') +
  xlab('Year') + 
  theme_classic() + ggpubr::grids()
GP

if(plt)ggsave(GP,filename=here(paste0('../plots/inc/incIMPcheck_Global.pdf')),w=7,h=5)

est <- merge(est,TBN[,.(iso3,year,c_newinc)],by=c('iso3','year'))
est[,gap:=e_inc_num - c_newinc]
est[gap< 0]
est[gap<0,summary(gap)]                 #mainly rounding
est[gap< -1e3]                          #except these
est[gap<0,gap:=0]

## CDRs
est[is.na(ocdr),ocdr:=c_newinc/e_inc_num] #TODO think dropped relapse jj
est[,tmp:=max(ocdr.sd,na.rm=TRUE),by=iso3] #use biggest sd
est[!is.finite(ocdr.sd),ocdr.sd:=tmp] #fill in
est[,tmp:=NULL]
est[,summary(ocdr)]
est[!is.finite(ocdr)]                   #0s
est[!is.finite(ocdr),ocdr:=1]                   #set to 1
est[ocdr==0]                   # rounding
est[ocdr==0,ocdr:=0.5]                   # rounding
est[,gap.sd:=gap * ocdr.sd / ocdr]
est[,summary(gap.sd/gap)]

tmp <- est[,.(gap=sum(gap)/1e6,gap2=sum(c_newinc*(1/ocdr-1))/1e6),by=year]
tmp                                     #TODO check new vs new + rel


GP <- ggplot(tmp,aes(year,gap)) + geom_line() +
  ylab('Undiagnosed TB incidence in millions')
GP

if(plt)ggsave(GP,filename=here::here('../plots/Gap.pdf'),w=7,h=5)

GP <- GP + geom_line(aes(year,gap2),col=2)    #TODO recheck
GP

if(plt)ggsave(GP,filename=here::here('../plots/GapCheck.pdf'),w=7,h=5) #think 1 lacks IND corr


tmpi <- merge(tmpi,est[,.(gap=sum(gap)/1e6),by=year],by='year')
tmpim <- melt(tmpi,id='year')
tmpim[variable=='incnum',variable:='total']
tmpim[variable=='gap',variable:='untreated']

GP2 <- ggplot(tmpim,aes(year,value,lty=variable)) + geom_line() +
  ylab('TB incidence (millions per year)') + expand_limits(y=0) +
  theme_classic() + ggpubr::grids()+
  theme(legend.position = c(0.1,0.9)) + xlab('Year')
GP2

if(plt)ggsave(GP2,file=here::here('../plots/Gap2.pdf'),w=7,h=5)


tmpi
GP3 <- ggplot(tmpi,aes(year,(1-gap/incnum))) + geom_line() +
  scale_y_continuous(label=percent) + expand_limits(y=0) +
  xlab('Year')+ ylab('Case Detection Ratio') +
  theme_classic() + ggpubr::grids()
GP3

if(plt)ggsave(GP3,file=here::here('../plots/Gap2cdr.pdf'),w=7,h=5)


save(est,file=here::here('../tmpdata/estg.Rdata'))
