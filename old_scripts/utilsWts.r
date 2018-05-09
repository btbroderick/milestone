library("quadprog")
library("Hmisc")


#######################################
#
# wtsFreq() - Estimates the prediction synthesis weights for the frequentist methods
#  
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   cTime - current time.
#   MM - Number of predictions to simulate based on past data
#   K - How many times to predict each simulated prediction for the mean computation
# Return - 
#   wts - prediction synthesis weights. First column is for the MSPE values and the 
#       second is for the vote method.
#
#######################################


wtsFreq<-function(dat,cTime,MM=100,K=10){
  pW<-pLN<-pG<-pLL<-y<-rep(0,MM)
  kpW<-kpLN<-kpG<-kpLL<-rep(0,K)
  wP<-lik(dat,cTime,"weibull")[c(2,1)]
  lnP<-lik(dat,cTime,"lnorm");
  gP<-lik(dat,cTime,"gompertz"); gP<-c(gP[2]/gP[1],gP[1])
  llP<-lik(dat,cTime,"llogis")[c(2,1)]
  lT<-length(dat[,1])/cTime
  ind<-which(dat[,2]>0)
  for(i in 1:MM){
    tms<-sample(dat[ind,2],2,FALSE)
    y[i]<-max(tms)
    cT<-min(tms)
    datTmp<-trnc(dat[,1:2],cT)
    nE<-sum(dat[ind,2]<=y[i])
    for(j in 1:K){
      jT<-cT+cumsum(rexp(2*nE,lT))
      kpW[j]<-weibTime(datTmp,params=wP, nE,cT,jT)
      kpLN[j]<-logNormTime(datTmp,params=lnP, nE,cT,jT)
      kpG[j]<-gompTime(datTmp,params=gP, nE,cT,jT)
      kpLL[j]<-logLogTime(datTmp,params=llP, nE,cT,jT)
    }
    pW[i]<-mean(kpW)
    pLN[i]<-mean(kpLN)
    pG[i]<-mean(kpG)
    pLL[i]<-mean(kpLL)
  }
  Am<-cbind(1,diag(4))
  b<-c(1,0,0,0,0)
  D<-cbind(pW,pLN,pG,pLL)
  wtsTmp<-table(apply(abs(sweep(D, 1, y, FUN = "-")),1,which.min))/MM  
  wts1<-rep(0,4)
  wts1[as.numeric(names(wtsTmp))]<-wtsTmp
  
  d<- t(D)%*%y
  D2<-t(D)%*%D  
  options(show.error.messages = FALSE)
  tst<-try( wts<-solve.QP(Dmat=D2, dvec=d, Amat=Am, bvec=b,meq=1)[[1]],T)
  options(show.error.messages = TRUE)
  if(!is.numeric(tst)){
    wts<-rep(0,4)
    wtsTmp<-apply(I(sweep(D, 1, y, FUN = "-")^2),2,mean)
    wts[which.min(wtsTmp)]=1
  }
  return(cbind(wts,wts1))
}


#######################################
#
# wtsBayes() - Estimates the prediction synthesis weights for the Bayesian methods
#  
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   cTime - current time.
#   MM - Number of predictions to simulate based on past data
#   K - How many times to predict each simulated prediction for the mean computation
#   W,LN,G,LL - MCMC chain of parameters for the weibull, LN, G, and LL methods 
#     respectively. Should be matrices of 2 columns. 
#     (First two columns generated from weibPost, logNormPost, gompPost or logLogPost.
#     perhaps with initial rows omited to allow for burn in.
# Return - 
#   wts - prediction synthesis weights. First column is for the MSPE values and the 
#       second is for the vote method.
#
#######################################


wtsBayes<-function(dat,cTime,MM=100,K=10,W,LN,G,LL){
  pW<-pLN<-pG<-pLL<-y<-rep(0,MM)
  kpW<-kpLN<-kpG<-kpLL<-rep(0,K)
  lnW<-length(W[,1])
  shpePost<-length(dat[,1])+1
  ind<-which(dat[,2]>0)
  for(i in 1:MM){
    tms<-sample(dat[ind,2],2,FALSE)
    y[i]<-max(tms)
    cT<-min(tms)
    datTmp<-trnc(dat[,1:2],cT)
    nE<-sum(dat[ind,2]<=y[i])
    for(j in 1:K){
      wh<-sample(1:lnW,1)
      wP<-W[wh,];  lnP<-LN[wh,];  gP<-G[wh,];  llP<-LL[wh,]; 
      lT<-rgamma(1,shape=shpePost,scale=1/cT)
      jT<-cT+cumsum(rexp(2*nE,lT))
      kpW[j]<-weibTime(datTmp,params=wP, nE,cT,jT)
      kpLN[j]<-logNormTime(datTmp,params=lnP, nE,cT,jT)
      kpG[j]<-gompTime(datTmp,params=gP, nE,cT,jT)
      kpLL[j]<-logLogTime(datTmp,params=llP, nE,cT,jT)
    }
    pW[i]<-mean(kpW)
    pLN[i]<-mean(kpLN)
    pG[i]<-mean(kpG)
    pLL[i]<-mean(kpLL)
  }
  Am<-cbind(1,diag(4))
  b<-c(1,0,0,0,0)
  D<-cbind(pW,pLN,pG,pLL)
  
  wtsTmp<-table(apply(abs(sweep(D, 1, y, FUN = "-")),1,which.min))/MM  
  wts1<-rep(0,4)
  wts1[as.numeric(names(wtsTmp))]<-wtsTmp
  
  d<- t(D)%*%y
  D2<-t(D)%*%D  
  options(show.error.messages = FALSE)
  tst<-try( wts<-solve.QP(Dmat=D2, dvec=d, Amat=Am, bvec=b,meq=1)[[1]],T)
  options(show.error.messages = TRUE)
  if(!is.numeric(tst)){
    wts<-rep(0,4)
    wtsTmp<-apply(I(sweep(D, 1, y, FUN = "-")^2),2,mean)
    wts[which.min(wtsTmp)]=1
  }
  return(cbind(wts,wts1))
} 


#######################################
#
# trnc() - Crops a dat matrix at thetm'th observed data point.
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   tM - Truncation point for prediction synthesis 
# Return - 
#   datRet: a new dat matrix with all arrivals after the cropping time removed, and any 
#      hitting times after the ceiling(P*nEvents)'th observation censored.
#######################################

trnc<-function(dat,tm){
  ind1<-which(dat[,1]<tm)
  datRet<-dat[ind1,]
  ind2<-which(datRet[,2]>tm)
  datRet[ind2,2]<-0
  
  return(datRet)
}

#######################################
#
# predInt() - Generates a confidence interval for the prediction synthesis method.
#
# Arguments -
#   wts - vector of 4 weights. All non-negative and sum to 1.
#   preds - MM x 4 matrix of predictions
#     Entries of wts and columns of pred must be in the same order.
#   alpha - size of the prediction interval.
# Return - 
#   confIntBounds - quantiles (alpha/2, 1-alpha/2) based.  
#
#######################################

predInt<-function(wts,preds,alpha=.05){
  its<-dim(preds)[1]
  W<-rep(wts,each=its)/its
  preds<-as.vector(preds)
  
  return( wtd.quantile(  preds,weights=W,probs=c(alpha/2,(1-alpha/2))   )   )
  
}

#######################################
#
# getQnts() - Generates quantiles for each method
#
# Arguments -
#   p - prediction values
#   wts - a 4 by 2 matrix with each column being a weight.
#   alpha - size of the prediction interval.
# Return - 
#   retInts - a 7 by 2 matrix of pred interval values for each individual and
#      prediction synthesis method.
#
#######################################

getQnts<-function(p,wts,alpha=.05){
  qs<-c(alpha/2,(1-alpha/2))
  its<-dim(p)[1]
  W1<-rep(wts[,1],each=its)
  W2<-rep(wts[,2],each=its)
  pp<-as.vector(p)
  retVal<-rbind( quantile(p[,1], probs = qs),
         quantile(p[,2], probs = qs),
         quantile(p[,3], probs = qs),
         quantile(p[,4], probs = qs),
         quantile(p, probs = qs),
         wtd.quantile(  pp,weights=W1,probs=qs),
         wtd.quantile(  pp,weights=W2,probs=qs) 
         )
  
  return(retVal)
  
}
#######################################
#
# symPredInt() - Generates a symmetric confidence interval for the prediction 
#     synthesis method.
#
# Arguments -
#   wts - vector of 4 weights. All non-negative and sum to 1.
#   preds - MM x 4 matrix of predictions
#     Entries of wts and columns of pred must be in the same order.
#   alpha - size of the prediction interval
# Return - 
#   confIntBounds - quantiles (alpha/2, 1-alpha/2) based.  
#######################################

symPredInt<-function(wts,preds,alpha=.05,MMM){
  mnPreds<-apply(preds,2,mean)
  mnVal<-mnPreds%*%wts
  devs<-abs(preds-mnVal)
  its<-dim(preds)[1]
  W<-rep(wts,each=its)/its
  preds<-as.vector(preds)
  marg<-wtd.quantile(  devs, weights=W, probs=(1-alpha)   ) 
  return(c(mnVal-marg, mnVal+marg )   )
  
}


#######################################
#
# getSymQnts() - Generates prediction intervals for each method using symmetric intervals
#
# Arguments -
#   p - prediction values
#   wts - a 4 by 2 matrix with each column being a weight.
#   alpha - size of the prediction interval.
# Return - 
#   retInts - a 7 by 2 matrix of pred interval values for each individual and
#      prediction synthesis method.
#
#######################################

getSymQnts<-function(p,wts, mns,alpha=.05){
  its<-dim(p)[1]
  W1<-rep(wts[,1],each=its)
  W2<-rep(wts[,2],each=its)
  pp<-as.vector(p)
  a<-1-alpha
  
  mrgns<-c( quantile(abs(p[,1]-mns[1]), probs = a),
              quantile(abs(p[,2]-mns[2]), probs = a),
              quantile(abs(p[,3]-mns[3]), probs = a),
              quantile(abs(p[,4]-mns[4]), probs = a),
              quantile(abs(p-mns[5]), probs = a),
              wtd.quantile(  abs(pp-mns[6]),weights=W1,probs=a),
              wtd.quantile(  abs(pp-mns[7]),weights=W2,probs=a) 
  )
  
  return(cbind(mns-mrgns, mns+mrgns))
  
}
