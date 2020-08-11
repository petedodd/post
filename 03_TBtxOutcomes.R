## TB treatment outcomes
rm(list=ls())
library(here)
source(here('00_utilities.R'))

if(! overwrite ){
  fn <- here('../tmpdata/TBOS.Rdata')
  if(file.exists(fn))stop('Not running as tmpdata/TBOS.Rdata exists!')
}


## load data
## WHO TB outcome data
TBO <- fread(here('../indata/TB_outcomes_2020-02-24.csv'))

## treatment outcomes=================
names(TBO)

TBO[iso3=='ZWE',.(iso3,year,new_sp_coh,new_sp_died,new_sp_def,
                  new_snep_coh,new_snep_died,new_snep_def)]
if(plt){
  for(rg in c('AFR','AMR','EMR','EUR','SEA','WPR')){
    GP <- ggplot(TBO[g_whoregion==rg],
                 aes(year,new_sp_died/new_sp_coh,group=iso3)) +
      geom_line() +
      geom_line(aes(year,new_sp_def/new_sp_coh,group=iso3),lty=2) +
      geom_line(aes(year,new_snep_died/new_snep_coh,group=iso3),col=2) +
      geom_line(aes(year,new_snep_def/new_snep_coh,group=iso3),col=2,lty=2) +
      facet_wrap(~iso3)
    ggsave(here::here(paste0('../plots/txo/Txo_',rg,'.pdf')))
  }
}

TBOS <- TBO[,.(Nc=sum(new_snep_coh,na.rm=TRUE)+sum(new_sp_coh,na.rm=TRUE),
               dc=sum(new_snep_died,na.rm=TRUE)+sum(new_sp_died,na.rm=TRUE),
               ltfu=sum(new_snep_def,na.rm=TRUE)+sum(new_sp_def,na.rm=TRUE)),
            by=.(iso3,g_whoregion)]


## assume same frac of LTFU is death
TBOS[,pd:=(dc + ltfu * (dc/(Nc-ltfu))) / Nc]
TBOS[,pd.sd:= sqrt(pd*(1-pd)/ Nc)]      #UNC

TBOS[,summary(1e2*pd.sd/pd)]
TBOS[,summary(1e2*pd)]                  #  0.000   3.633   5.850   7.036   8.418  50.000 
TBOS[,qplot(pd)]
TBOS[pd>.3]
(tbonacns <- TBOS[is.na(pd),iso3])

cat(tbonacns,file=here::here('texto/tbonacns.txt'))

TBORS <- TBOS[is.finite(pd),.(pd=mean(pd),
                              pd.sd1=mean(pd.sd^2),pd.sd2=var(pd)),
              by=g_whoregion]
TBORS[,pd.sd:=sqrt(pd.sd1+pd.sd2)]       #total variance
TBORS[,c('pd.sd1','pd.sd2'):=NULL]

for(rg in TBOS[,unique(g_whoregion)]){
  TBOS[g_whoregion==rg & is.na(pd),pd.sd:=TBORS[g_whoregion==rg,pd.sd]]
  TBOS[g_whoregion==rg & is.na(pd),pd:=TBORS[g_whoregion==rg,pd]]
}
TBOS

TBOS[Nc==0,.(iso3,pd,g_whoregion)]

save(TBOS,file=here('../tmpdata/TBOS.Rdata'))
## end tx outcomes ======
