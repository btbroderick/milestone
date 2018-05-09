library("stats")
library("VGAM")


source('./utilsBayes1.r')
source('./utilsFreq.r')
source('./utilsWts.r')

rgompertz<-flexsurv::rgompertz
dgompertz<-flexsurv::dgompertz
pgompertz<-flexsurv::pgompertz
qgompertz<-flexsurv::qgompertz


N<-1500  #number to simulate for the full data set.
nE<-500 # landmark event number

lambda=.5
mnlog<-5
sdlog<-.25

#Create full data set
t<-cumsum( rexp(N,rate=lambda))
x<-rlnorm(N,meanlog=mnlog,sdlog=sdlog)
dat<-cbind(t,(t+x))


sort(dat[,2])[nE] #Landmark event time

#truncate to ceiling(P*nE) events
dat<-truncDatSim(dat,nE,P=.5)

#Priors
wP<-c(1,50,1,50)
lnP<-c(1,50,1,50)
gP<-c(1,50,1,50)
llP<-c(1,50,1,50)



cTime<-max(dat)

#Frequentist Predictions
freqRes<-getFreqInts(dat,nE,MM=200,totAcc=510)
# freqRes[[1]] are the 7 predictions in the same order 
#     as the paper tables
# freqRes[[2]] are the prediction intervals in the paper
# freqRes[[3]] is a symmetric prediction interval not included
#     in the paper

#Bayes predictions
BayesRes<-getBayesInt(dat,nE,wP,lnP,gP,llP,MM=800,totAcc=510)    
# BayesRes[[1]] are the 7 predictions in the same order 
#     as the paper tables
# BayesRes[[2]] are the prediction intervals in the paper
# BayesRes[[3]] is a symmetric prediction interval not included
#     in the paper
