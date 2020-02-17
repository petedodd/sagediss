## This file is for generating age/sex disaggregation of incidence.
##
## This is therefore to be run after the analysis of prevalence data and after the child TB data has been generated
##
##
##
## #Read in necessary data
##
## Relevant libraries for loading and manipulating data:
##

rm(list=ls())
library(knitr)                          # for table formatting
library(metafor)                        # for MA of sex splits in kids
library(here)                           # file locations
library(gridExtra)                      # for grid.arrange


## Flag for whether cairo device is functioning for outward plotting:
gypt <- TRUE
plotting <- TRUE                       #whether to plot ot not
estyr <- 2018
set.seed(1234)                          #PRNG seed
source(here('R/utilities.R'))

## Various pre-requisite data
load(here('Rdata/est.Rdata'))         # current burden estimates
load(here('Rdata/tb.Rdata'))          # includes notifications
load(here('Rdata/pop.Rdata'))          # population for denom
load(here('R/02children/plots/K.Rdata'))   # Child estimates
load(here('prevsplits/prevsplits.Rdata'))   # from prevalence split work
names(prevsplits)[2] <- 'age_group'
load(here('R/02children/mortinput/BB.Rdata'))     #for age splits

## a country hash
hbcsh <- merge(unique(est[,.(iso3,g.hbc)]),unique(tb[,.(iso3,name=country)]))
hbcsh <- hbcsh[g.hbc==TRUE]; hbcsh[,g.hbc:=NULL]
## shorter names for graphs
hbcsh[iso3=='COD',name:="DR Congo"]
hbcsh[iso3=='PRK',name:="DPR Korea"]
hbcsh[iso3=='TZA',name:="UR Tanzania"]



## Reforming BB for kid splits and adding in regional averages where missing
BB <- as.data.table(BB)[,.(iso3,aa,ab,g.whoregion)]
BB[,propu5:=aa/(aa+ab)];BB[,propu5.sd:=sqrt(aa*ab/(aa+ab+1))/(aa+ab)];
BBR <- BB[,.(propu5=mean(propu5),propu5.sd=mean(propu5.sd),
             aa=mean(aa),ab=mean(ab)),by=g.whoregion]
BBR <- merge(est[year==estyr][,.(iso3,g.whoregion)],BBR,by='g.whoregion',
             all.x=TRUE,all.y = TRUE,allow.cartesian = TRUE)
setkey(BBR,iso3)
(missing <- setdiff(BBR$iso3,BB$iso3))
BB <- rbind(BB,BBR[missing])            #add missing with regional averages
BB[,g.whoregion:=NULL]
setkey(BB,iso3)


## Having a look at using the model for kids prior. First needs some NAs filled
K <- merge(K,BBR[,.(iso3,g.whoregion)],by='iso3')
K[,rav1:=mean(propu15b,na.rm = TRUE),by=g.whoregion]
K[,rav2:=mean(propu15b.sd,na.rm = TRUE),by=g.whoregion]
K[is.na(propu15b),propu15b:=rav1]; K[is.na(propu15b.sd),propu15b.sd:=rav2]
K[propu15b==0,propu15b.sd:=rav2]; K[propu15b==0,propu15b:=rav1]; 
K[,propu15:=propu15b]; K[,propu15.sd:=propu15b.sd]; #using prog model


## First, for countries with `NA` in prevsplit (from method fails, especially in AMR)
## create a global median of proportions, normalised:
names(prevsplits)[names(prevsplits)=='age.group'] <- 'age_group'
(missed <- prevsplits[is.na(prop),unique(iso3)])
glop <- prevsplits[,.(propm=median(prop,na.rm=TRUE), #global median
                      propm.se=median(prop.se,na.rm=TRUE)),by=.(sex,age_group)]
tot <- glop[,sum(propm)]; glop[,propm:=propm/tot] #renormalize
prevsplits <- merge(prevsplits,glop,by=c('sex','age_group'),
                     all.x=TRUE,all.y=TRUE,allow.cartesian = TRUE)
prevsplits[iso3 %in% missed,
           c('prop','prop.se'):=list(propm,propm.se)] #fill out with global version
(missed <- prevsplits[is.na(prop.se),unique(iso3)])
prevsplits[iso3 %in% missed,prop.se:=propm.se] #fill out with global version
prevsplits[,c('propm','propm.se'):=NULL]
## prevsplits[iso3=='ZAF']                 #there
## prevsplits[iso3=='TKM']

prevsplits$age_group <- factor(prevsplits$age_group,levels=agz2,ordered=TRUE)
prevsplits$sex <- factor(prevsplits$sex,levels=c('M','F'),ordered=TRUE)
setkeyv(prevsplits,c('iso3','age_group','sex'))

##
## # Child sex ratios
## 
## Looking at sex ratios in the child notification data
## Make data:
tmpw <- tb[year==estyr,c('iso3','g.whoregion',
                        outer(c('newrel.m','newrel.f'),kdz,paste0)),with=FALSE]
tmpw[,N:=newrel.m04 + newrel.f04 + newrel.m514 + newrel.f514 ]
tmpw <- tmpw[!is.na(N)]
tmpw <- tmpw[N>0]

## Use metafor and conduct some meta-analyses for younger and older age groups
## 
## ======== 0-4 ===== proportion M
dat <- tmpw[,.(id=iso3,xi=newrel.m04,ni=newrel.m04 + newrel.f04,g.whoregion)]
dat <- dat[dat$ni>0,]
dat <- dat[dat$ni>dat$xi,]
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
citmp <- t(apply(dat[,c('xi','ni')],1,function(x) binom.test(x[1], x[2])$conf.int))
dat$ci.lb <- citmp[,1]
dat$ci.ub <- citmp[,2]
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
res.rma <- rma.glmm(measure="PLO", mods=~g.whoregion, xi=xi, ni=ni, data=dat) #regional
## ## plot
## with(dat, forest(yi, ci.lb=ci.lb, ci.ub=ci.ub, refline=predict(res, transf=transf.ilogit)$pred,slab = id))
## addpoly(res, row=-1, transf=transf.ilogit)
## title('Proportion male 0-4')
## abline(h=0)
(YS <- predict(res, transf=transf.ilogit, digits=3))
## by region
YSR <- predict(res.rma,transf=transf.ilogit, digits=3)
YSR <- data.table(mid=YSR$pred,hi=YSR$cr.ub,lo=YSR$cr.lb,g.whoregion=dat$g.whoregion)
(YSR <- YSR[,.(mid=mid[1],hi=hi[1],lo=lo[1]),by=g.whoregion])

## ======== 5-14 ===== proportion M
dat <- tmpw[,.(id=iso3,xi=newrel.m514,ni=newrel.m514 + newrel.f514,g.whoregion)]
dat <- dat[dat$ni>0,]
dat <- dat[dat$ni>dat$xi,]
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
citmp <- t(apply(dat[,c('xi','ni')],1,function(x) binom.test(x[1], x[2])$conf.int))
dat$ci.lb <- citmp[,1]
dat$ci.ub <- citmp[,2]
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
res.rma <- rma.glmm(measure="PLO", mods=~g.whoregion, xi=xi, ni=ni, data=dat) #regional

## ## plot
## with(dat, forest(yi, ci.lb=ci.lb, ci.ub=ci.ub, refline=predict(res, transf=transf.ilogit)$pred,slab = id))
## addpoly(res, row=-1, transf=transf.ilogit)
## title('Proportion male 5-14')
## abline(h=0)
(OS <- predict(res, transf=transf.ilogit, digits=3))
## by region
OSR <- predict(res.rma,transf=transf.ilogit, digits=3)
OSR <- data.table(mid=OSR$pred,hi=OSR$cr.ub,lo=OSR$cr.lb,g.whoregion=dat$g.whoregion)
(OSR <- OSR[,.(mid=mid[1],hi=hi[1],lo=lo[1]),by=g.whoregion])


##
## # Simplification of flow
##
## ## Pass1: countries using standard adjustments
##
## Make a list of countries to be done
isotodo <- est[,unique(iso3)]
length(isotodo)

## Use standard adjustment for standard adjustment or CR with good CDR
##
## CDRs of countries
cdrtab <- est[,.(iso3,year,inc.num=inc*pop)]
##newinc includes unk h
cdrtab <- merge(cdrtab,tb[,.(iso3,year,c.newunk,c.newinc)],by=c('iso3','year'))
cdrtab <- cdrtab[!is.na(c.newunk)]
cdrtab <- cdrtab[,.SD[year==max(year)],by=iso3]
dim(cdrtab)
cdrtab[,cdr:=c.newinc/inc.num*1e5]
cdrtab[,qplot(cdr)] + geom_vline(xintercept = 0.85,col=2)
cdrtab[!is.finite(cdr)]                      #
cdrtab[!is.finite(cdr),cdr:=1]
cdrtab[cdr>1]
setkey(cdrtab,iso3)

## Lists of countries todo
isosf <- est[source.inc=="Case notifications,\nStandard adjustment",
             as.character(unique(iso3))]
isocr <- est[source.inc=="Capture-recapture",
             as.character(unique(iso3))]
cdrtab[isocr]
cdrtab[isosf]


## Delete any existing plots in relevant folder!
## 
fz <- list.files(path=here('incsplits/plots'),pattern="*.pdf",
                 recursive = TRUE, full.names = TRUE)
if(length(fz)>0) do.call(file.remove,list(fz)) #clean up

## _Pass 1_ will be countries that are either:
## 
##  * SF adjustment & $CDR>0.85$
##  * CR &  have $CDR>0.85$
##  * fewer than 1000 notifications in total
## 
isopass1 <- c(cdrtab[isosf][cdr>0.85,iso3], #SF + CDR>0.85
              cdrtab[isocr][cdr>0.85,iso3], #CR + CDR>0.85 
              cdrtab[c.newinc<1e3,iso3]     #few notifications
              )
isopass1 <- unique(isopass1)
isopass1
cat(isopass1,file=here('incsplits/data/pass1.txt'))

## Work up incsplits for pass 1 countries
H <- W <- 7
ISL <- hi30 <- pltlist <- noteissue <- list()

## starting loop
for(cn in isopass1){
  tmp <- getNotes(cn)                   #most recent notifications
  if(tmp[,sum(newrel,na.rm=TRUE)]<50){
    tmp[is.na(newrel),newrel:=0] #correcting GRD, eg
    noteissue[[cn]] <- cn
  }
  cdr <- cdrtab[cn,cdr]
  tot.newrel <- tmp[,sum(newrel)]
  if(tot.newrel==0){
    tot.newrel <- 1
    cdr <- 1
  }
  tmp[,prop:=newrel/tot.newrel]
  tmp[,inc:=newrel/cdr]
  tmp[,c('iso3','prop.se','src','pass','problem'):=
         .(cn,NA,'Notifications x factor',1,FALSE)]
  tmp <- merge(tmp,AA[,.(age_group,age)],by='age_group')
  ISL[[cn]] <- tmp
  plt <- niplot(tmp)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- here('incsplits/plots/pass1') + glue('/{cn}1.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,h=H,w=W)
    } else
      ggsave(filename=fn,plot=plt,h=H,w=W)
  }
}

## Following countries with fewer than 50 notifications have some issues:
(noteissue <- unlist(noteissue))
cat(isopass1,file=here('incsplits/data/noteissue.txt'))

## Still left todo
isotodo <- setdiff(isotodo,isopass1)
length(isotodo)

incsplit <- rbindlist(ISL)
incsplit[,sum(inc)]*1e-6
incsplit[age_group %in% kds,sum(inc)]/incsplit[,sum(inc)]*1e2
incsplit[,length(unique(iso3))]

##
## ## Pass 2: sampling against prior
##
## 

priorsample <- function(cn,
                        nsim=1e6,
                        silent=FALSE,
                        tol=1e-2,
                        adultadjust=FALSE
                        ){
  tmp <- getNotes(cn,silent=silent)                 #Notifications
  if(tmp$year[1] < max(tb$year)){# if not most recent year then rescale
    EL <- est[cn][,.(iso3,year,inc,pop)]
    sf <- EL[year==max(tb$year),inc*pop]/EL[year==tmp$year[1],inc*pop]
    if(!is.finite(sf)) sf <- 1
    tmp$newrel <- tmp$newrel * sf
  }
  ## fetch adult child notifications, total incidence
  INC <- est[cn][year==estyr,inc*pop/1e5] #total incidence TODO check OK
  atmp <- tmp[age_group %ni% kds,]    #adult notifications
  tmp2 <- K[cn]
  ## --- make splits ----
  if(nsim>1){ ## replace with samples
    ## child splits
    lnp <- getLNparms(tmp2[,propu15],tmp2[,propu15.sd]^2)
    propu15 <- rlnorm(n=nsim,meanlog=lnp[1,1],sdlog = sqrt(lnp[2,1])) #adult/child
    aa <- BB[iso3==cn,aa]; ab <- BB[iso3==cn,ab]
    suppressWarnings({propu5 <- rbeta(nsim,aa,ab)})#child age split
    tt <- getLNparms(YS$pred,(YS$ci.lb-YS$ci.ub)^2/3.92^2)
    pu5m <- rlnorm(nsim,tt[1,1],sqrt(tt[2,1])) #young child sex split
    tt <- getLNparms(OS$pred,(OS$ci.lb-OS$ci.ub)^2/3.92^2)
    po5m <- rlnorm(nsim,tt[1,1],sqrt(tt[2,1])) #old child sex split
    ## adult splits
    if(!adultadjust){
      lnp <- getLNparms(prevsplits[cn][,prop],prevsplits[cn][,prop.se]^2)
      adultprops <- apply(t(lnp),1,
                          function(x)rlnorm(n=nsim,meanlog=x[1],sdlog = sqrt(x[2])))
    } else {
      adultprops <- matrix(tmp$newrel[-c(1:4)],nrow=nsim,ncol=12,byrow=TRUE)
    }
    adultprops <- adultprops / rowSums(adultprops) #normalize
    ## make matrix nsim x age/sex of proportions
    PZ <- cbind(
      propu15 * propu5 * pu5m,   #M<5
      propu15 * propu5 * (1-pu5m), #F<5
      propu15 * (1-propu5) * po5m, #M514
      propu15 * (1-propu5) * (1-po5m), #F514
      (1-propu15) * adultprops
                )
    PZ <- PZ/rowSums(PZ)                #safety
  } else {stop('nsim must be >1!');}
  ## --- compare to notifications ---
  NZ <- matrix(tmp$newrel,nrow=nsim,ncol=16,byrow=TRUE)
  NZ <- INC*PZ - NZ                     #hoping >0
  keep <- apply(NZ,1,function(x) all(x>0,na.rm = TRUE)) #NA for cases like MOZ
  accept <- TRUE
  if(all(!keep)){                       #all failed
    if(!silent)cat('None accepted, using tolerance!\n')
    accept <- FALSE
    NZ[is.na(NZ)] <- 0
    NZ[NZ>0] <- 0                       #only penalize undershoot
    penalty <- rowSums(NZ^2)
    threshold <- quantile(penalty,tol) #best 1%
    keep <- penalty < threshold
    if(!silent)cat('Using best ',sum(keep),' samples\n')
  } else{ if(!silent)cat('Accepted ',sum(keep),' samples\n');}
  if(sum(keep)>0){                      #safety
    if(sum(keep)>1){
      PZ <- PZ[keep,]
      PZ <- colMeans(PZ)
    } else PZ <- PZ[keep,]
    PZ <- PZ/sum(PZ)
    tmp[,prop:=PZ]
    PZ <- INC * PZ
    tmp[,inc:=PZ]
    tmp[,prop.se:=NA]
    if(accept) tmp[,src:='Rejection sampling from prior']
    if(!accept) tmp[,src:='Best samples from prior']
  } else {                              #error!
    tmp[,c('prop','inc','prop.se','src'):=NA]
  }
  tmp
}

## test
tmp <- priorsample('NGA')
niplot(tmp)


## Looping through countries
##
## 
for(cn in isotodo){
  cno <- priorsample(cn)
  cno$pass <- 2
  ISL[[cn]] <- cno
  plt <- niplot(cno)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- here('incsplits/plots/pass2') + glue('/{cn}2.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,w=W,h=H)
    } else
      ggsave(filename=fn,plot=plt,w=W,h=H)
  }
}



## # Pass 3: Ad hoc reversion to adjustments
##
##  Oddities among the pass2 countries deserving special treatment:
##
isopass3 <- c('BGD','PRK','MMR','IND','CHN')
cno <- priorsample('BGD',adultadjust = TRUE)
niplot(cno)
cno <- priorsample('PRK',adultadjust = TRUE)
niplot(cno)
cno <- priorsample('MMR',tol=1e-4,adultadjust = TRUE)
niplot(cno)
cno <- priorsample('IND',adultadjust = TRUE)
niplot(cno)
cno <- priorsample('CHN',adultadjust = TRUE)
niplot(cno)


## Sample following the adult pattern of notifications
for(cn in isopass3){
  cno <- priorsample(cn,adultadjust = TRUE)
  cno$pass <- 3
  ISL[[cn]] <- cno
  plt <- niplot(cno)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- here('incsplits/plots/pass3') + glue('/{cn}3.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,w=W,h=H)
    } else
      ggsave(filename=fn,plot=plt,w=W,h=H)
  }
}


## move them into the failed directory
for( cn in isopass3){
  ofn <- here('incsplits/plots/pass2') + glue('/{cn}2.pdf')
  nfn <- here('incsplits/plots/pass2/0pass2fails') + glue('/{cn}2.pdf')
  file.rename(from=ofn,to=nfn)
}

## # Finishing
##
## aggregations
incsplit <- rbindlist(ISL,fill=TRUE)

## checks
incsplit[,sum(inc,na.rm=TRUE)]*1e-6
incsplit[age_group %in% kds,sum(inc,na.rm=TRUE)]/incsplit[,sum(inc,na.rm=TRUE)]*1e2
incsplit[,length(unique(iso3))]

## proportions
incsplit[,testsum:=sum(prop),by=iso3]
(tf <- unique(incsplit[abs(testsum-1)>1e-2,iso3]))
incsplit[iso3 %in% tf,.(iso3,newrel,prop,src)]
incsplit[iso3 %in% tf,all(newrel==0)]   #only failing as zero
setkey(incsplit,iso3)

## incidence
tmp <- incsplit[,.(testinc=sum(inc)),by=iso3]
tmp2 <- merge(est[year==estyr,.(iso3,inc)],pop[year==estyr,.(iso3,e.pop.num)],by='iso3')
tmp2[,inc.num:=inc*e.pop.num/1e5]
tmp <- merge(tmp,tmp2[,.(iso3,inc.num)],by='iso3',all.x = TRUE,all.y = FALSE)
incsplit <- merge(incsplit,tmp2[,.(iso3,inc.num)],
                  by='iso3',all.x = TRUE,all.y = FALSE)
(tf <- unique(tmp[abs(testinc/inc.num-1)>1e-2,iso3]))
incsplit[iso3 %in% tf,inc:=prop*inc.num,by=iso3]
incsplit[,inc.num:=NULL]


## Save out data:
save(file=here('incsplits/data/incsplit.Rdata'),incsplit)

