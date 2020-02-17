library(data.table)                     # data handling
library(ggplot2)                        # plotting
library(scales)                         # plotting
library(ggpubr)                         # plot arrangents
library(glue)                           # strings

## A function for making absolute labels with space separators
absspace <- function(x,...) {             #works
   format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}

## not int
`%ni%` <- Negate(`%in%`)

## function for time stamps in right format:
getm <- function()glue(format(Sys.time(), "%Y_%m_%d_%H%M"))
getm()

## Function for date stamps in right format:
getdate <- function()glue(gsub('-','_',Sys.Date()))
getdate()


## Some useful age range vectors:
agz <- c('04','514','1524','2534','3544','4554','5564','65') #for extract
agz2 <- c('0_4','5_14','15_24','25_34','35_44','45_54','55_64','65plus') #for labels
agz3 <- gsub('_','-',agz2)
agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
kds <- c('0_4','5_14')
kdz <- c('04','514')
AA <- data.table(a1=agz,age_group=agz2,age=agz3,agegt=agz4) #for conversion


## A function get notification data (new & relapse) by age for a given country
##
## Needs tb.Rdata to be loaded and est.Rdata and pop.Rdata
## 
##
getNotes <- function(cn,m=TRUE,mostrecent=TRUE,silent=FALSE){
  if(!silent) cat('Country=',cn,'\n')
  notes <- matrix(nrow=2,ncol=8)
  yr <- myr <- max(tb$year)
  if(mostrecent){
    while(all(is.na(tb[iso3==cn&year==yr,paste0('newrel.m',agz),with=FALSE]))
          & yr>1980){
      yr <- yr-1
      if(!silent) cat('Check notification data! Using a previous year due to NAs!\n')
    }
  }
  here <- tb[cn]
  tmp <- unlist(here[year==yr,paste0('newrel.m',agz),with=FALSE])
  notes[1,] <- tmp
  tmp <- unlist(here[year==yr,paste0('newrel.f',agz),with=FALSE])
  notes[2,] <- tmp
  colnames(notes) <- agz2
  rownames(notes) <- c('M','F')
  if(m){                              #melt first
    notes <- data.table(newrel=c(notes),
                        sex=rep(c('M','F'),ncol(notes)),age_group=rep(agz2,each=2))
    notes$age_group <- factor(notes$age_group,levels=agz2,ordered=TRUE)
    notes$sex <- factor(notes$sex,levels=c('M','F'),ordered=TRUE)
  }
  notes[,c('iso3','year'):=.(cn,yr)]
  setkeyv(notes,c('age_group','sex'))
  notes
}

## ## test
## print(getNotes('ZAF'))
## print(getNotes('AUT'))
## print(getNotes('AUT',mostrecent=FALSE))
## print(getNotes('BES'))


niplot <- function(cndt,label.problem=FALSE){
  cn <- cndt[,iso3][1]
  cndt[,age:=gsub("_","-",age_group)];
  cndt$age <- factor(cndt$age,levels=agz3,ordered=TRUE)
  cndt$sex <- factor(cndt$sex,levels=c('M','F'),ordered=FALSE)
  ans <- ggplot(data=cndt,aes(x=age,y=newrel,fill=sex)) +
    coord_flip() + grids() +
    geom_bar(data=cndt[sex=='M'],stat='identity',aes(y=newrel))+
    geom_bar(data=cndt[sex=='F'],stat='identity',aes(y=-newrel))+
    ylab('New & relapse cases') + xlab('Age group')+
    scale_y_continuous(labels = absspace)
  if('inc' %in% names(cndt)){
    ans <- ans +
      geom_bar(data=cndt[sex=='M'],stat='identity',
               aes(x=age,y=inc),fill='transparent',col=1)+
      geom_bar(data=cndt[sex=='F'],stat='identity',
               aes(x=age,y=-inc),fill='transparent',col=1)
  }
  ans <- ans + scale_x_discrete(labels=agz4)
  nm <- tb[cn,country[1]]
  if(cn=='CIV') nm <- "Cote d'Ivoire"   #otherwise encoding pbm
  tlt <- nm
  if('problem' %in% names(cndt)){
    if(cndt[,any(problem)] & label.problem)
      tlt <- paste0(tlt,': -- PROBLEM!')
  }
  ans <- ans + ggtitle(tlt)
  ans
}

## niplot(getNotes('ZAF'))


rowSD <- function(M) sqrt(rowMeans((M-rowMeans(M))^2))
colSD <- function(M) sqrt(colMeans((M-colMeans(M))^2))

getLNparms <- function(m,v){
  lvv <- log(1+v/m^2)
  rbind(mu=log(m)-lvv/2,sig2=lvv)
}

## ## test
## getLNparms(1,.2)
## getLNparms(rep(1,3),rep(.2,3))
## tt <- getLNparms(1,1)
## mean(rlnorm(1e4,tt[1,1],sqrt(tt[2,1])))
## sd(rlnorm(1e4,tt[1,1],sqrt(tt[2,1])))


ftb <- Vectorize(function(x){
  stopifnot (!is.character(x))
  if (!is.na(x)){
    smallpos <- x>0 & x<0.01
    dg <- ifelse(abs(x)>0.01 & abs(x)<100, 2, 3)
    x <- signif(x, 3)
    x <- format(x, digits=dg, nsmall=0, big.mark=" ", justify='right', drop0trailing = T)
    if (smallpos) x <- '<0.01'
  } else x <- '-'
  return(x)}, 'x')



## Function for rounding anything left with a decimal point and removing '<'
## (for number formatting for tables)
zapdotlt <- function(x){
  x <- gsub("<","",x)
  got <- grep("\\.",x)
  x[got] <- as.character(round(as.numeric(x[got])))
  x
}

## Function for formatting big numbers >=1e6, which are breaking with ftb
fmtbig <- function(x){
  x <- signif(x, 3)
  x <- format(x, digits=3,
              nsmall=0, big.mark=" ", justify='right',
              drop0trailing = TRUE,scientific=FALSE)
  ## x <- stringi::stri_enc_toutf8(x)      #ensure UTF8
  x
}
## fmtbig(123451000)

