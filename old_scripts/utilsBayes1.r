
library("msm") #needed for truncated normal.


#######################################
#
# weibTime() - generates random time to achieve a given number of events 
#     based on a Weibull distribution for event occurance
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   params - (lambda, k)   the parameters for the Weibull distribution. 
#       'lambda' is the scale and 'k' is the shape for the Weibull distribution
#   nEvents - number of events needed to complete the study.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   joinTime - The times in which new subject arrivals occur. 
# Return - 
#   Time (in days) when nEvents is reached. 
#######################################

weibTime<-function(dat,params, nEvents,currTime,joinTime){
  if(length(params)!=2){
    return(0)
  }
  lambda<-params[1]; k<-params[2];
  indCens<-which(dat[,2]==0);
  chkVal<- currTime-dat[indCens,1]
  t<-length(indCens)
  
  P<-pweibull(chkVal, shape=k,scale=lambda) 
  dat[indCens,2]<-dat[indCens,1] + qweibull(I(P+(1-P)*runif(t)), shape=k,scale=lambda) 
  if(min(joinTime)>=0){
    dat<-rbind(dat, cbind(joinTime,(joinTime+rweibull(length(joinTime),shape=k,scale=lambda))))
  }
  return(sort(dat[,2])[nEvents])
}


#######################################
#
# weibPost() - generates an MCMC chain of parameters and predictions using the Metropolis Hastings method
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   priParams - [muLambda,varLambda, muk,vark] mean and variance for the prior parameters
#     of the Weibull. (for a gamma distribution)
#   M - the number of Monte Carlo iterations
#   nEvents - number of events needed to complete the study.
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - A matrix of size M x 3 
#   Columns 1 & 2 : lam, k - sequence of parameter values of the MCMC with the starting point 
#     being the MLE.
#   Column 3: rndTime - sequence of random times for milestone completion based on a single 
#     simulation for each of the parameter values in lam,k
#######################################

weibPost<-function(dat,currTime=max(dat),priParams,M=1000,nEvents=length(dat[,1]),shouldjT="T",totAcc=2*nEvents){
  Btr<-200
  shpePost<-length(dat[,1])+1
  ind_i<-which(dat[,2]>0)  #ind_i are uncensored values
  n_i=length(ind_i)
  needed=nEvents-n_i
  x_i<-dat[ind_i,2]-dat[ind_i,1]
  
  lpx_i=sum(log(x_i))
  x_j<-currTime-dat[-ind_i,1]
  x<-c(x_i,x_j)
  
  # Starting values
  wP<-lik(dat,currTime,"weibull")[c(2,1)]
  K<-K_c<-wP[2]
  lam<-l_c<-wP[1]
  
  B_L<-priParams[1]/priParams[2];  A_L<-priParams[1]^2/priParams[2];
  B_K<-priParams[3]/priParams[4];  A_K<-priParams[3]^2/priParams[4];
  f_theta_c<-dgamma(l_c, A_L, B_L)*dgamma(K_c, A_K, B_K)
  lT<-rgamma(1,shape=shpePost,scale=1/currTime)
  simN=totAcc-dim(dat)[1]
  simN=max(c(simN,0))
  if(simN==0){shouldJT=0}
  if(shouldjT){ jT<-currTime+cumsum(rexp(simN,rate=lT)) }
  else{jT=-1}
  
  rndTime<-weibTime(dat,params=c(l_c,K_c), nEvents,currTime,jT)
  
  
  for(i in 1:M){
    K_t<-rgamma(1,(K_c*Btr),Btr)
    l_t<-rgamma(1,(l_c*Btr),Btr)
    
    gt_c<-dgamma(K_c,(K_t*Btr),Btr)*dgamma(l_c,(l_t*Btr),Btr)
    gc_t<-dgamma(K_t,(K_c*Btr),Btr)*dgamma(l_t,(l_c*Btr),Btr)
    f_theta_t<-dgamma(l_t, A_L, B_L)*dgamma(K_t, A_K, B_K)
    
    op<-sum(I((x/l_c)^K_c ))-sum(I((x/l_t)^K_t ))
    logRat<-n_i*log(K_t*l_c^K_c/(K_c*l_t^K_t))+(K_t-K_c)*lpx_i
    fRat<-exp(logRat+op)
    
    p=fRat*f_theta_t*gc_t/(f_theta_c*gt_c)
    u<-runif(1,0,1)
    if(is.finite(p)){
      if(u<p){
        K_c=K_t; l_c=l_t;  f_theta_c<-f_theta_t
      }
    }
    K<-c(K,K_c);    lam=c(lam,l_c)
    lT<-rgamma(1,shape=shpePost,scale=1/currTime)
    mnLambda<-1/lT
    mostWeib<-qweibull(.99, shape=K_c,scale=l_c) 
    if(shouldjT){    jT<-currTime+cumsum(rexp(simN,rate=lT))}
    rndTime<-c(rndTime,weibTime(dat,params=c(l_c,K_c), nEvents,currTime,jT) )
    if(i/100==round(i/100)){
      c<-(mean(diff(lam[(i-99):i])==0)/.7)
      Btr<-c*Btr
      
    }
    
  }
  return(cbind(lam,K,rndTime))
}


#######################################
#
# loglogTime() - generates random time to achieve a given number of events 
#     based on a log-logistic distribution for event occurance
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is the time of an event.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   params - (alpha,beta)   the parameters for the log-logistic distribution. 
#       'alpha' is the scale and 'beta' is the shape parameter
#   nEvents - number of events needed to complete the study.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   joinTime - The times in which new subject arrivals occur. 
# Return - 
#   Time (in days) when nEvents is reached. 
#######################################

logLogTime<-function(dat,params, nEvents,currTime,jT){
  if(length(params)!=2){
    return(0)
  }
  alpha<-params[1]; beta<-params[2];
  indCens<-which(dat[,2]==0);
  chkVal<- currTime-dat[indCens,1]
  P<-pllogis(chkVal, scale=alpha, shape=beta)
  t<-length(indCens)
  dat[indCens,2]<-dat[indCens,1] + qllogis(I(P+(1-P)*runif(t)), scale=alpha, shape=beta)
  if(min(jT)>=0){
    dat<-rbind(dat, cbind(jT,(jT+ rllogis(length(jT), scale=alpha, shape=beta)  )))
  }
  return(sort(dat[,2])[nEvents])
}


#######################################
#
# logLogPost() - generates an MCMC chain of parameters and predictions using the Metropolis Hastings method
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   priParams - [mualpha,varalpha, mub,varb] mean and variance for the prior parameters
#     of the log-logistic distribution. (for a gamma prior)
#   M - the number of Monte Carlo iterations
#   nEvents - number of events needed to complete the study.
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - A matrix of size M x 3 
#   Columns 1 & 2 alpha, b - sequence of parameter values of the MCMC with the 
#       starting point being the MLE.
#   Column 3: rndTime - sequence of random times for milestone completion based on a single simulation 
#     for each of the parameter values in alpha,b
#######################################

logLogPost<-function(dat,currTime=max(dat),priParams,M=1000,nEvents=length(dat[,1]) ,shouldjT="T",totAcc=2*nEvents){
  ind_i<-which(dat[,2]>0)  #ind_i are uncensored values
  shpePost<-length(dat[,1])+1
  n_i=length(ind_i)
  x_i<-dat[ind_i,2]-dat[ind_i,1]
  sx_i=sum(x_i)
  x_j<-currTime-dat[-ind_i,1]
  Btr<-200  
  # Starting values
  llP<-lik(dat,currTime,"llogis")[c(2,1)]
  a<-a_c<-llP[1]; b<-b_c<-llP[2]
  
  B_a<-priParams[1]/priParams[2];  A_a<-priParams[1]^2/priParams[2];
  B_b<-priParams[3]/priParams[4];  A_b<-priParams[3]^2/priParams[4];
  
  f_theta_c<-dgamma(b_c, A_b, B_b)*dgamma(a_c, A_a, B_a)
  lT<-rgamma(1,shape=shpePost,scale=1/currTime)
  simN=totAcc-dim(dat)[1]
  simN=max(c(simN,0))
  if(simN==0){shouldJT=0}
  if(shouldjT){ jT<-currTime+cumsum(rexp(simN,rate=lT)) }
  else{jT=-1}
  rndTime<-logLogTime(dat,params=c(a_c,b_c), nEvents,currTime,jT)
  
  for(i in 1:M){
    a_t<-rgamma(1,(a_c*Btr),Btr)
    b_t<-rgamma(1,(b_c*Btr),Btr)
    
    gt_c<-dgamma(a_c,(a_t*Btr),Btr)*dgamma(b_c,(b_t*Btr),Btr)
    gc_t<-dgamma(a_t,(a_c*Btr),Btr)*dgamma(b_t,(b_c*Btr),Btr)
    f_theta_t<-dgamma(b_t, A_b, B_b)*dgamma(a_t, A_a, B_a)
    
    op<-n_i*(log(b_t*a_c/(b_c*a_t))+(1-b_t)*log(a_t)+(b_c-1)*log(a_c))
    op2<-sum(I((b_t-b_c)*log(x_i)+2*log(1+(x_i/a_c)^b_c )-2*log(1+(x_i/a_t)^b_t ))) 
    op3<-sum(I( log(1+(x_j/a_c)^b_c )-log(1+(x_j/a_t)^b_t ) ))
    fRat<-exp(op+op2+op3)
    
    p=fRat*f_theta_t*gc_t/(f_theta_c*gt_c)
    u<-runif(1,0,1)
    
    if(is.finite(p)){
      if(u<p){
        a_c=a_t; b_c=b_t;  f_theta_c<-f_theta_t
      }
    }
    a<-c(a,a_c)
    b<-c(b,b_c) 
    lT<-rgamma(1,shape=shpePost,scale=1/currTime)
    if(shouldjT){    jT<-currTime+cumsum(rexp(simN,rate=lT))}
    rndTime<-c(rndTime,logLogTime(dat,params=c(a_c,b_c), nEvents,currTime,jT))
    if(i/100==round(i/100)){
      c<-(mean(diff(b[(i-99):i])==0)/.7)
      Btr<-c*Btr
      
    }
  }
  return(cbind(a,b,rndTime))
  
}



#######################################
#
# gompTime() - generates random time to achieve a given number of events 
#     based on a Gompertz distribution for event occurance
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is the time of an event.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   params - (eta,b)   the parameters for the Gompertz distribution. 
#       'b' is the scale and 'eta' is the shape parameter
#   nEvents - number of events needed to complete the study.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   joinTime - The times in which new subject arrivals occur. 
# Return - 
#   Time (in days) when nEvents is reached. 
#######################################

gompTime<-function(dat,params, nEvents,currTime,jT){
  if(length(params)!=2){
    return(0)
  }
  eta<-params[1]; b<-params[2];
  indCens<-which(dat[,2]==0);
  chkVal<- currTime-dat[indCens,1]
  t<-length(indCens)
  P<-pgompertz(chkVal, shape=b, rate = (eta*b))
  dat[indCens,2]<-dat[indCens,1] + qgompertz(I(P+(1-P)*runif(t)), shape=b, rate = (eta*b))
  if(min(jT)>=0){
    dat<-rbind(dat, cbind(jT,(jT+ rgompertz(length(jT), shape=b, rate = (eta*b))  )))
  }
  return(sort(dat[,2])[nEvents])
}

#######################################
#
# gompPost() - generates an MCMC chain of parameters and predictions using the Metropolis Hastings method
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   priParams - [muEta,varEta, mub,varb] mean and variance for the prior parameters
#     of the Gompertz distribution. (for a gamma prior)
#   M - the number of Monte Carlo iterations
#   nEvents - number of events needed to complete the study.
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - A matrix of size M x 3 
#   Columns 1 & 2: eta,b - sequence of parameter values of the MCMC with the 
#       starting point being the MLE.
#   Column 3: rndTime - sequence of random times for milestone completion based on a single simulation 
#     for each of the parameter values in eta,b
#######################################
gompPost<-function(dat,currTime=max(dat),priParams,M=1000,nEvents=length(dat[,1]),shouldjT="T",totAcc=2*nEvents){
  ind_i<-which(dat[,2]>0)  #ind_i are uncensored values
  n_i=length(ind_i)
  shpePost<-length(dat[,1])+1
  x_i<-dat[ind_i,2]-dat[ind_i,1]
  Btr<-200
  sx_i=sum(x_i)
  x_j<-currTime-dat[-ind_i,1]
  x<-c(x_i,x_j)
  N<-length(x)
  
  # Starting values fomr MLE 
  
  gP<-lik(dat,currTime,"gompertz");
  gP<-abs(c(gP[2]/gP[1],gP[1]))
  
  eta<-eta_c<-gP[1]
  b<-b_c<-gP[2]
  
  B_eta<-priParams[1]/priParams[2];  A_eta<-priParams[1]^2/priParams[2];
  B_b<-priParams[3]/priParams[4];  A_b<-priParams[3]^2/priParams[4];
  
  
  f_theta_c<-dgamma(b_c, A_b, B_b)*dgamma(eta_c, A_eta, B_eta)
  lT<-rgamma(1,shape=shpePost,scale=1/currTime)
  
  simN=totAcc-dim(dat)[1]
  simN=max(c(simN,0))
  if(simN==0){shouldJT=0}
  
  if(shouldjT){ jT<-currTime+cumsum(rexp(simN,rate=lT)) }
  else{jT=-1}
  rndTime<-gompTime(dat,params=c(eta_c,b_c), nEvents,currTime,jT=jT)
  
  for(i in 1:M){
    eta_t<-rgamma(1,(eta_c*Btr),Btr)
    b_t<-rgamma(1,(b_c*Btr),Btr)
    
    gt_c<-dgamma(eta_c,(eta_t*Btr),Btr)*dgamma(b_c,(b_t*Btr),Btr)
    gc_t<-dgamma(eta_t,(eta_c*Btr),Btr)*dgamma(b_t,(b_c*Btr),Btr)
    f_theta_t<-dgamma(eta_t, A_eta, B_eta)*dgamma(b_t, A_b, B_b)
    
    op<-(b_t-b_c)*sx_i+eta_c*sum(I(exp(b_c*x)))-eta_t*sum(I(exp(b_t*x)))
    op2<-n_i*(log((b_t*eta_t/(b_c*eta_c))))
    
    fRat<-exp(N*(eta_t-eta_c) +op+op2)
    p=fRat*f_theta_t*gc_t/(f_theta_c*gt_c)
    u<-runif(1,0,1)
    
    
    if(is.finite(p)){
      if(u<p){
        eta_c=eta_t; b_c=b_t;  f_theta_c<-f_theta_t
      }
    }
    eta<-c(eta,eta_c)
    b<-c(b,b_c)
    lT<-rgamma(1,shape=shpePost,scale=1/currTime)
    if(shouldjT){    jT<-currTime+cumsum(rexp(simN,rate=lT))}
    rndTime<-c(rndTime,gompTime(dat,params=c(eta_c,b_c), nEvents,currTime,jT=jT))
    if(i/100==round(i/100)){
      c<-(mean(diff(eta[(i-99):i])==0)/.7)
      Btr<-c*Btr
      
    }
  }
  return(cbind(eta,b,rndTime))
  
}


#######################################
#
# logNormTime() - generates random time to achieve a given number of events 
#     based on a log-normal distribution for event occurance
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is the time of an event.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   params - (mu,sigma)   the parameters for the log-logistic distribution. 
#       mu and sigma are the mean and standard deviation for the logarithm of
#       the event time.
#   nEvents - number of events needed to complete the study.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   joinTime - The times in which new subject arrivals occur. 
# Return - 
#   Time (in days) when nEvents is reached. 
#######################################


logNormTime<-function(dat,params, nEvents,currTime,jT){
  if(length(params)!=2){
    return(0)
  }
  mu<-params[1]; sigma<-params[2];
  
  indCens<-which(dat[,2]==0);
  chkVal<- currTime-dat[indCens,1]
  P<-plnorm(chkVal, meanlog =mu, sdlog = sigma)
  t<-length(indCens)
  dat[indCens,2]<-dat[indCens,1] + qlnorm(I(P+(1-P)*runif(t)), meanlog =mu, sdlog = sigma)
  if(min(jT)>=0){
    dat<-rbind(dat,cbind(jT,(jT+rlnorm(length(jT), meanlog =mu, sdlog = sigma))))
  }
  return(sort(dat[,2])[nEvents])
}

#######################################
#
# logNormPost() - generates an MCMC chain of parameters and predictions using the Metropolis Hastings method
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   currTime - current Time. Must be greater than any value in dat, 
#       and less than any time in joinTime  
#   priParams - [muMu,varMu, musigma,varsigma] mean and variance for the prior parameters
#     of the log-normal  distribution. (normal prior for mu, gamma prior for sigma)
#   M - the number of Monte Carlo iterations
#   nEvents - number of events needed to complete the study.
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - A matrix of size M x 3 
#   Columns 1 & 2 mu,sigma - sequence of parameter values of the MCMC with the 
#       starting point being the MLE.
#   Column 3: rndTime - sequence of random times for milestone completion based on a single simulation 
#     for each of the parameter values in mu, sigma
#######################################

logNormPost<-function(dat,currTime=max(dat),priParams,M=1000,nEvents=length(dat[,1]),shouldjT="T",totAcc=2*nEvents){
  ind_i<-which(dat[,2]>0)  #ind_i are uncensored values
  shpePost<-length(dat[,1])+1
  n_i=length(ind_i)
  x_i<-dat[ind_i,2]-dat[ind_i,1]
  Btr=200;sTr=.25
  
  lpx_i=sum(log(x_i))
  lpx_i2=sum(I(log(x_i)^2))
  lx_j<-log(currTime-dat[-ind_i,1])
  
  # Starting values
  lnP<-lik(dat,currTime,"lnorm");
  sig<-sig_c<-lnP[2]
  mu<-mu_c<-lnP[1]
  
  muPrior<-priParams[1]
  sigPrior<-priParams[2]
  
  B_sig<-priParams[3]/priParams[4];  A_sig<-priParams[3]^2/priParams[4];
  DD_c<-prod( pnorm(I(mu_c-lx_j)/(sqrt(2)*sig_c)))
  f_theta_c<-dnorm(mu,muPrior,sigPrior)*dgamma(sig,A_sig,B_sig)
  lT<-rgamma(1,shape=shpePost,scale=1/currTime)
  
  simN=totAcc-dim(dat)[1]
  simN=max(c(simN,0))
  if(simN==0){shouldJT=0}
  
  if(shouldjT){ jT<-currTime+cumsum(rexp(simN,rate=lT)) }
  else{jT=-1}
  rndTime<-logNormTime(dat,params=c(mu_c,sig_c), nEvents,currTime,jT)
  for(i in 1:M){
    mu_t<-rnorm(1,mu_c,sTr)
    sig_t<-rgamma(1,(sig_c*Btr),Btr)
    g<-dgamma(sig_t,(sig_c*Btr),Btr)/dgamma(sig_t,(sig_c*Btr),Btr)
    f_theta_t<-dnorm(mu_t,muPrior,sigPrior)*dgamma(sig_t,A_sig,B_sig)
    op<-n_i*(mu_c^2/(2*sig_c^2)-mu_t^2/(2*sig_t^2))+(mu_t/(sig_t^2)-mu_c/(sig_c^2))*lpx_i +
      (1/(2*sig_c^2)-1/(2*sig_t^2))*lpx_i2
    DD_t<-prod( pnorm(I(mu_t-lx_j)/(sqrt(2)*sig_t)))
    
    fRat<-exp(op)*DD_t/DD_c*(sig_c/sig_t)^n_i
    
    p=fRat*f_theta_t*g/f_theta_c
    u<-runif(1,0,1)
    
    if(is.finite(p)){
      if(u<p){
        mu_c=mu_t; sig_c=sig_t; f_theta_c<-f_theta_t; DD_c=DD_t
      }
    }
    mu<-c(mu,mu_c)
    sig<-c(sig,sig_c)
    lT<-rgamma(1,shape=shpePost,scale=1/currTime)
    if(shouldjT){    jT<-currTime+cumsum(rexp(simN,rate=lT))}
    rndTime<-c(rndTime,logNormTime(dat,params=c(mu_c,sig_c), nEvents,currTime,jT))
    if(i/100==round(i/100)){
      c<-(mean(diff(sig[(i-99):i])==0)/.7)
      Btr<-c*Btr
      sTr<-sTr/c
      
    }
    
  }
  return(cbind(mu,sig,rndTime))
}

#######################################
#
# truncDatSim() - Crops a dat matrix at the ceiling(nEvents*P)'th observed data point.
#  (Not used in published program)
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   nEvents - number of events needed to complete the study.
#   P - The proportion of events for cropping
# Return - 
#   dat: a new dat matrix with all arrivals after the cropping time removed, and any 
#      hitting times after the ceiling(P*nEvents)'th observation censored.
#######################################

truncDatSim<-function(dat,nEvents,P){
  N=ceiling(nEvents*P)
  tot<-sum(dat[,2]>0)
  if(N>tot){
    return(dat[,1:2])
  }
  ind<-which(dat[,2]>0)
  currTime<-sort(dat[ind,2])[N]
  ind2<-which(dat[,1]>currTime)
  if(length(ind2)){
    dat<-dat[-ind2,]
  }
  ind3<-which(dat[,2]>currTime)
  dat[ind3,2]<-0
  return(dat  )
}



#######################################
#
# getBayes() - Computes all Bayesian predictions for the nE observation using dat.
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   nE - number of events needed to complete the study.
#   wP,lnP,gP,llP - Prior parameters for Weibull, log-normal, gompertz, and log-logistic 
#       distributions (respectfully).
#   MM - number of observations in the MCMC chain for each posterior distribution
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - 
#  preds:  a vector givin the predictions for each Bayesian method in the following order
#         Weibull
#         log-normal
#         gompertz
#         log-logistic
#         Prediction Synthesis: - MSPE
#                                 Average
#                                 Vote
#######################################

getBayes<-function(dat,nE,wP,lnP,gP,llP,MM=5000,totAcc=2*nE){
  MMstrt<-round(MM/3)
  cTime<-max(dat)
  preds<-rep(0,7)
  predsA<-weibPost(dat,currTime=cTime,priParams=wP,M=MM,nEvents=nE,1,totAcc )
  preds[1]<-mean(predsA[MMstrt:MM,3])
  predsB<-logNormPost(dat,currTime=cTime,priParams=lnP,M=MM,nEvents=nE,1,totAcc )
  preds[2]<-mean(predsB[MMstrt:MM,3])
  predsC<-gompPost(dat,currTime=cTime,priParams=gP,M=MM,nEvents=nE,1,totAcc)
  preds[3]<-mean(predsC[MMstrt:MM,3])
  predsD<-logLogPost(dat,currTime=cTime,priParams=llP,M=MM,nEvents=nE,1,totAcc )
  preds[4]<-mean(predsD[MMstrt:MM,3])
  
  tmp<-wtsBayes(dat,cTime,MM=100,K=6,predsA[,1:2],predsB[,1:2],predsC[,1:2],predsD[,1:2])
  preds[5]<-mean(preds[1:4])
  preds[6]<-tmp[,1]%*%preds[1:4]
  preds[7]<-tmp[,2]%*%preds[1:4]
  return(preds)
}


#######################################
#
# getBayesInt() - Computes all Bayesian predictions with prediction intervals
#   for the nE observation using dat.
#
# Arguments - 
#    dat - matrix(n, 2) where n are the number of individuals currently in the study.
#         First column is the time (in days) the individual joined the study.
#         Second column is time the hitting time was reached.  
#         If the hitting time has not currently been reached,
#         the value should be set to 0.
#   nE - number of events needed to complete the study.
#   wP,lnP,gP,llP - Prior parameters for Weibull, log-normal, gompertz, and log-logistic 
#       distributions (respectfully).
#   MM - number of observations in the MCMC chain for each posterior distribution
#   totAcc - the total number of subjects to accrue. Must be greater than nE, otherwise the 
#       program will not terminate.
# Return - 
#  preds:  a vector givin the predictions for each Bayesian method in the following order
#         Weibull
#         log-normal
#         gompertz
#         log-logistic
#         Prediction Synthesis: - MSPE
#                                 Average
#                                 Vote
#######################################

getBayesInt<-function(dat,nE,wP,lnP,gP,llP,MM=5000,totAcc=2*nE,alpha=.05,addjT=1){
  MMstrt<-round(MM/3)
  cTime<-max(dat)
  preds<-rep(0,7)
  
  predsA<-weibPost(dat,currTime=cTime,priParams=wP,M=MM,nEvents=nE,shouldjT=addjT,totAcc)
  preds[1]<-mean(predsA[MMstrt:MM,3])
  predsB<-logNormPost(dat,currTime=cTime,priParams=lnP,M=MM,nEvents=nE,shouldjT=addjT,totAcc)
  preds[2]<-mean(predsB[MMstrt:MM,3])
  predsC<-gompPost(dat,currTime=cTime,priParams=gP,M=MM,nEvents=nE,shouldjT=addjT,totAcc )
  preds[3]<-mean(predsC[MMstrt:MM,3])
  predsD<-logLogPost(dat,currTime=cTime,priParams=llP,M=MM,nEvents=nE ,shouldjT=addjT,totAcc)
  preds[4]<-mean(predsD[MMstrt:MM,3])
  
  tmp<-wtsBayes(dat,cTime,MM=100,K=6,predsA[,1:2],predsB[,1:2],predsC[,1:2],predsD[,1:2])
  preds[5]<-mean(preds[1:4])
  preds[6]<-tmp[,1]%*%preds[1:4]
  preds[7]<-tmp[,2]%*%preds[1:4]
  
  pInts1<-getQnts(cbind(predsA[MMstrt:MM,3],predsB[MMstrt:MM,3],
                        predsC[MMstrt:MM,3],predsD[MMstrt:MM,3]),tmp,alpha)
  
  pInts2<-getSymQnts(cbind(predsA[MMstrt:MM,3],predsB[MMstrt:MM,3],
                        predsC[MMstrt:MM,3],predsD[MMstrt:MM,3]),tmp,preds,alpha)

  
  return(list(preds,pInts1,pInts2))
}


