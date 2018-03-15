#######################################
#
# lik() -Generates the MLE of the survival data for the 
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   cTime - Current time
#   whchDist - which distribution to use for MLE model fitting
#         "weibull"	  Weibull	      
#         "llogis"	  Log-logistic	
#         "lnorm"	    Log-normal	 
#         "gompertz"	Gompertz	
# Return - 
#   param - MLE fit of parameters in the order and type returned by 'flexsurvreg' 
#
#######################################

lik<-function(dat,cTime,whchDist="weibull"){
  ind<-which(dat[,2]>0)
  y<-dat[ind,2]-dat[ind,1] 
  y<-c(y,cTime-dat[-ind,1])
  hit<-c(rep(1,length(ind)),c(rep(0,length(dat[-ind,1])) ))
  ind<-which(y==0)
  if(length(ind)>0){  y=y[-ind]; hit=hit[-ind]}
  
  return(flexsurvreg(Surv(y,hit)~1, dist=whchDist, method="Nelder-Mead")$res[,1])
}


#######################################
#
# getFreqVals() - Computes all frequentist predictions for the nE observation milestone
#     for dat.
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   nE - number of events needed to complete the study.
#   MM - number of predictions to simulate to estimate the mean prediction for the 
#       individual models.
# Return - 
#  preds:  a vector givin the predictions for each frequentist method in the following order
#         Weibull
#         log-normal
#         gompertz
#         log-logistic
#         Prediction Synthesis: - MSPE
#                                 Average
#                                 Vote
#######################################

getFreqVals<-function(dat,nE,MM=1000,shouldJT=1){
  cTime<-max(dat)
  wP<-lik(dat,cTime,"weibull")[c(2,1)]
  lnP<-lik(dat,cTime,"lnorm");
  gP<-lik(dat,cTime,"gompertz"); gP<-c(gP[2]/gP[1],gP[1])
  llP<-lik(dat,cTime,"llogis")[c(2,1)]
  pW<-pLn<-pG<-pLL<-0
  lT<-length(dat[,1])/cTime
  jT<- -1
  for(k in 1:MM){        
    if(shouldJT){
      jT<-cTime+cumsum(rexp(2*nE,lT))
    }
    pW<-pW+ weibTime(dat,params=wP, nE,cTime,jT)/MM
    pLn<-pLn+ logNormTime(dat,params=lnP, nE,cTime,jT)/MM
    pG<-pG+ gompTime(dat,params=gP, nE,cTime,jT)/MM
    pLL<-pLL+ logLogTime(dat,params=llP, nE,cTime,jT)/MM
  }
  tmp<-wtsFreq(dat,cTime,100,8)
  pWt<-tmp[,1]%*%c(pW,pLn,pG,pLL)
  pAve<-mean(c(pW,pLn,pG,pLL))
  pWt2<-tmp[,2]%*%c(pW,pLn,pG,pLL)
  preds<-c(pW,pLn,pG,pLL,pAve,pWt,pWt2)
  return(preds)
}



#######################################
#
# getFreqInts() - Computes all frequentist predictions for the nE observation milestone
#     for dat along with prediction intervals.
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   nE - number of events needed to complete the study.
#   MM - number of predictions to simulate to estimate the mean prediction for the 
#       individual models.
#   alpha - size of the prediction interval.
# Return - 
#  preds:  a vector givin the predictions for each frequentist method in the following order
#         Weibull
#         log-normal
#         gompertz
#         log-logistic
#         Prediction Synthesis: - MSPE
#                                 Average
#                                 Vote
#######################################

getFreqInts<-function(dat,nE,MM=1000,alpha=.05,shouldJT=1){
  cTime<-max(dat)
  wP<-lik(dat,cTime,"weibull")[c(2,1)]
  lnP<-lik(dat,cTime,"lnorm")
  gP<-lik(dat,cTime,"gompertz")
  gP<-abs(c(gP[2]/gP[1], gP[1]))
  llP<-lik(dat,cTime,"llogis")[c(2,1)]
  pW<-pLn<-pG<-pLL<-0
  lT<-length(dat[,1])/cTime
  p<-matrix(0,MM,4)
  jT<- -1
  for(k in 1:MM){        
    if(shouldJT){
      jT<-cTime+cumsum(rexp(2*nE,lT))
    }
    p[k,1]<-weibTime(dat,params=wP, nE,cTime,jT)
    p[k,2]<-logNormTime(dat,params=lnP, nE,cTime,jT)
    p[k,3]<-gompTime(dat,params=gP, nE,cTime,jT)
    p[k,4]<-logLogTime(dat,params=llP, nE,cTime,jT)
  }

  tmp<-wtsFreq(dat,cTime,100,8)
  mnPreds<-apply(p,2,mean)
  pWt<-tmp[,1]%*%mnPreds
  pAve<-mean(mnPreds)
  pWt2<-tmp[,2]%*%mnPreds
  preds<-c(mnPreds,pAve,pWt,pWt2)
  pInts1<-getQnts(p,tmp,alpha)
  pInts2<-getSymQnts(p,tmp,preds,alpha)
  
  
  return(list(preds,pInts1,pInts2))
}
