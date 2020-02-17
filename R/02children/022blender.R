##
## This is to be run after the model estimates have been produced
## Not a lot not now happens here
## 
## #Read in necessary data
##
## NB this is the old ensemble version which not only needs updating, but comparing with alternative approaches.
##
## Relevant libraries for loading and manipulating data:
##
rm(list=ls())
library(data.table)                     #data handling
library(ggplot2)                        #plotting
library(knitr)                          #for table formatting
library(here)

## Read in background data
load(here('Rdata/est.Rdata')) #estimates data
load(here('Rdata/pop.Rdata')) #population data
load(here('Rdata/tb.Rdata'))    #tb notifications
P <- fread(here('R/02children/plots/PJDout.csv')) #model estimates

est <- merge(est,pop[,.(iso3,year,e.pop.num)],by=c('iso3','year'))

## This year:
estyr <- 2018


## Generate SD and merge in:
P[,median.sd:=(UQ-LQ)/1.36]
## merge in
dim(P)

H <- merge(est[year==estyr,
               .(iso3,inc.num=inc*e.pop.num/1e5,
                 inc.num.sd=1e-5*inc*e.pop.num*inc.sd/inc,inc)],
           P[,.(iso3,median,median.sd)],by='iso3',all=TRUE) #NB -- check dims!

##
## # Calculations
##
## Deal with NAs
H[is.na(median) ,median:=0]



## ## Proportions in children
## 
H[,propu15b:=median/inc.num]
H[,propu15b.sd:=propu15b*sqrt(median.sd^2/median^2 + inc.num.sd/inc.num^2)]
kable(H[1:10,.(iso3,propu15b,propu15b.sd)])          #check

## Deal with NAs: replace with simple regression model against per capita incidence:
H[,sum(is.na(propu15b))]                 #how many NA?
mod1 <- lm(propu15b~inc,data=H)          #model
pp <- predict(mod1,newdata=H[is.na(propu15b),],interval = 'predict')
H[is.na(propu15b),propu15b.sd:=(pp[,3]-pp[,2])/3.92]
H[is.na(propu15b),propu15b:=pp[,'fit']]
kable(H[1:10,.(iso3,propu15b,propu15b.sd)])   #check


##
## ## Outputs
##
## Save out the relevant bits of H
##
K <- H[,.(iso3,inc.num,inc.num.sd,
          propu15b,propu15b.sd)]
save(K,file=here('R/02children/plots/K.Rdata'))
