#**********************************************************
# A Bayesian Approach to Diagnosing Covariance Matrix Shifts
#**********************************************************
whitening <- function(xI = x0, xII = x){
  muI <- apply(xI,2,mean)
  SigmaI <- cov(xI)
  eigenSigmaI <- eigen(SigmaI)
  Lambda <- eigenSigmaI$values
  TV <- eigenSigmaI$vectors
  x0 = t(apply(xI,1,function(xI) TV%*%solve(diag(sqrt(Lambda)))%*%t(TV)%*%(xI-muI)))
  muII <- apply(xII,2,mean)
  x = t(apply(xII,1,function(xII) TV%*%solve(diag(sqrt(Lambda)))%*%t(TV)%*%(xII-muII)))
  barx <- round(apply(x,2,mean),3)
  hatSigma <- abs(round(cov(x),2))
  return(list('xI'=x0,'xII'=x))
}
BDC_diag <- function(x0,x,h,Nruns=5000,tau=2000,ifwhitening = TRUE){
  p = dim(x0)[2]
  N = dim(x0)[1]
  n = dim(x)[1]
  #*************************************************************
  if(ifwhitening){
    whiteningRes <- whitening(xI = x0, xII = x)
    x0 = whiteningRes$xI
    x = whiteningRes$xII
  }
  barx0 = apply(x0,2,mean);
  hatSigma0 = abs(cov(x0));
  barx = apply(x,2,mean);
  hatSigma = abs(cov(x));
  #*************************************************************
  # list(barx0,barx,hatSigma0,hatSigma)
  invhatSigma0 = solve(hatSigma0)
  S0 = (N-1)*hatSigma0;S0
  S = (n-1)*hatSigma;S                         
  #*************************************************************
  # Gibbs sampling       *****************************
  #*************************************************************
  # h = (n-1)/qchisq(0.05,n-1) - 1;h
  S2 <- ifelse(S<(n-1),n-1,S)
  Delt <- S2*h;
  # origin value         *****************************
  invSigma = invhatSigma0                        
  delta.M = matrix(rep(0,p*p),p,p);
  prob.Delta <- matrix(rep(0.5,p*p),p,p)
  # start iteration      *****************************
  counter.delta.M = list()
  for(i in 1:Nruns){
    # sample  mu^i         *****************************
    V = solve(n*invSigma+N*invhatSigma0)
    g = V%*%(n*invSigma%*%barx+N*invhatSigma0%*%barx0)
    mu = rmvnorm(1,g,V)
    # sample invSigma     *****************************
    v=N+n-1
    delt <- ifelse(delta.M==1,Delt,0)
    S_delta = S0 + delt
    Scale = solve(n*t(barx-mu)%*%(barx-mu)+S+S_delta)
    invSigma <- rwish(v,Scale)
    # sample Delta        *****************************
    for(l in 1:p){      # row
      for(j in l:p){    # column
        #******************************************
        pr = vector()
        for(k in 1:2){
          delta.M[l,j] = delta.M[j,l] = k-1
          delt2 <- ifelse(delta.M==1,Delt,0)
          S_delta2 = S0 + delt2
          prob.Delta2 = ifelse(delta.M==0,1-prob.Delta,prob.Delta)
          log.S.part = ((N-1)/2)*log(det(S_delta2))
          log.tr.part = -0.5*sum(diag(invSigma%*%S_delta2))
          log.pd.part = log(prod(prob.Delta2))
          r_k = log.S.part + log.tr.part + log.pd.part
          pr[k] = r_k
        }
        prob.r_0 = 1/(1+exp(pr[2]-pr[1]));
        delta.smaple = sample(c(0,1),1,replace = TRUE,prob = c(prob.r_0,1-prob.r_0))
        delta.M[l,j] = delta.smaple
        delta.M[j,l] = delta.smaple
      }
    }
    counter.delta.M[[i]] <- delta.M
  }
  #*************************************************************
  #   calculate probability of delta****************************
  #*************************************************************
  res = matrix(rep(0,p*p),p,p)
  for(m in (tau+1):Nruns){
    res = counter.delta.M[[m]] + res
  }
  res.prob <- res/length((tau+1):Nruns)
  # if(probability){
  #   res = res.prob
  # }else{
  #   res0 = res.prob[!upper.tri(res.prob)]
  #   res = as.vector(ifelse(res0>0.5,1,0))
  # }
  return(res.prob)
}
