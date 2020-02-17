## this uses the stuff in Xdatanfunction.R
rm(list=ls())
library(here)
source(here('R/02children/Xpaeddatanfunctions.R'))

## running...
LHsample <- randomLHS(1e4,LHN)                #latin hypercube sample
LHparms <- makeParameters2(LHsample,metaparm) #turn the raw LHS into parameters

## qbeta warnings...
system.time({                           #around 100 sec for 1e4
    LAcl <- loopAns2(LHparms,metaparm,method=c('comm','lat'))
})

Dcl <- make2cat(LAcl)                      #aggregate output
testof <- as.data.table(Dcl)

## make df summing yng and old
testob <- testof[,list(value = sum(value)),by=list(country, replicate, iso3)]


## country data
testoc <- testob[,list(median=rme(value),LQ=rlq(value),UQ=ruq(value)),by=.(iso3,country)]


## ------- write out ---------
write.table(testoc,file=here('R/02children/plots/PJDout.csv'),
            sep=',',col.names=TRUE,row.names=FALSE)
