#***********************************************
#  
#***********************************************
# p = 6;delta=1
Sigma_OC1 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(p*delta^2)-Delta
  }
  delta <- uniroot(func, c(0,5), Delta = Delta, p = p, tol = 0.01)$root
  D <- diag(delta,p)
  D <- round(D,4)
  return(D)
}
Sigma_OC2 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(2*(delta*0.5)^2+(p-2)*delta^2)-Delta
  }
  delta <- uniroot(func, c(0,5), Delta = Delta, p = p, tol = 0.01)$root
  c(delta/2,delta/2,rep(delta,p-2))
  D <- diag(c(delta/2,delta/2,rep(delta,p-2)))
  D <- round(D,4)
  return(D)
}
Sigma_OC3 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(12*delta^2)-Delta
  }
  delta <- uniroot(func,c(0,5),Delta = Delta, p = p,tol=0.01)$root
  D = matrix(0,p,p)
  
  D[1,2] = D[1,3] = delta
  D[2,1] = D[2,3] = D[2,4] = delta
  D[3,1] = D[3,2] = D[3,4] = delta
  D[4,2] = D[4,3] = delta
  
  D[1,4] = D[4,1] = delta
  # D[1,3] = delta
  # D[2,4] = delta
  # D[3,1] = delta
  # D[4,2] = delta
  # 
  # D[1,2] = D[2,3] = D[3,4] = delta
  # D[2,1] = D[3,2] = D[4,3] = delta
  
  D <- round(D,4)
  return(D)
}
# Sigma_OC3(p = 5,Delta = 1.5)
Sigma_OC4 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(10*delta^2)-Delta
  }
  delta <- uniroot(func,c(0,5),Delta = Delta, p = p,tol=0.01)$root
  D = matrix(0,p,p)
  D[1,2] = D[1,3] = delta
  D[2,1] = D[2,4] = delta
  D[3,1] = D[3,4] = delta
  D[4,2] = D[4,3] = delta
  D[2,3] = D[3,2] = delta
  
  # D[1,2] = delta
  # D[2,1] = D[2,3] = delta
  # D[3,2] = D[3,4] = delta
  # D[4,3] = delta
  D <- round(D,4)
  return(D)
}
Sigma_OC5 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(9*delta^2)-Delta
  }
  delta <- uniroot(func,c(0,5),Delta = Delta, p = p,tol=0.01)$root
  D = matrix(0,p,p)
  
  D[1,1] = D[1,2] = D[1,3] = delta
  D[2,1] = D[2,2] = D[2,3] = delta
  D[3,1] = D[3,2] = D[3,3] = delta
  D <- round(D,4)
  return(D)
}
Sigma_OC6 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(6*delta^2+3*(delta*0.5)^2)-Delta
  }
  delta <- uniroot(func,c(0,5),Delta = Delta, p = p,tol=0.01)$root
  D = matrix(0,p,p)
  D[1,1] = D[1,2] = delta
  D[2,1] = D[2,3] = delta
  D[3,2] = D[3,3] = delta
  D[1,3] = D[2,2] = D[3,1] = delta*0.5
  D <- round(D,4)
  return(D)
}
Sigma_OC7 <- function(p,Delta){
  func <- function(delta,Delta,p){
    sqrt(8*delta^2)-Delta
  }
  delta <- uniroot(func,c(0,5),Delta = Delta, p = p,tol=0.01)$root
  D = matrix(0,p,p)
  D[1,1] = D[1,2] = delta
  D[2,1] = D[2,2] = delta
  D[p-1,p-1] = D[p-1,p] =delta
  D[p,p-1] = D[p,p] = delta
  D <- round(D,4)
  return(D)
}
# p = 7; N = 200; n = 200; Delta = 1;  OC=7
GDATA <- function(p, N, n, Delta, OC, MeanCov=FALSE){
  SigmaOC <- switch(OC,
                    Sigma_OC1(p,Delta),
                    Sigma_OC2(p,Delta),
                    Sigma_OC3(p,Delta),
                    Sigma_OC4(p,Delta),
                    Sigma_OC5(p,Delta),
                    Sigma_OC6(p,Delta),
                    Sigma_OC7(p,Delta)
  );SigmaOC
  mu = rep(0,p);mu                                 
  Sigma0 = diag(rep(1,p));Sigma0
  x0 = rmvnorm(N,mu,Sigma0)
  
  Sigma = Sigma0 + SigmaOC;Sigma                   
  x = rmvnorm(n,mu,Sigma);
  if(MeanCov){
    res <- list('xI' =x0,'xII'=x, 'delta' = c(rep(0,p),SigmaOC[upper.tri(SigmaOC,diag = TRUE)]))
  }else{
    res <- list('xI' =x0,'xII'=x, 'delta' = SigmaOC[upper.tri(SigmaOC,diag = TRUE)])
  }
  return(res)
}
simBDC <- function(h, p, N, n, Delta, OC, index = c("CENE","ErrorRate")){
  dat <- GDATA(p , N, n, Delta, OC, MeanCov = FALSE)
  x0 <- dat$xI
  x <- dat$xII
  shift = ifelse(dat$delta!=0,1,0)
  hNum = length(h)
  # ress = matrix(0,hNum,p*(p+1)/2)
  res = rep(0,hNum*2)
  for(i in 1:hNum){
    BDC0 <- BDC_diag(x0,x,h=h[i],Nruns=5000,tau=2000,ifwhitening = FALSE)
    BDC0 = BDC0[upper.tri(BDC0,diag = TRUE)]
    BDC0 = as.vector(ifelse(BDC0>0.5,1,0))
    if(index=="ErrorRate"){
      res[i*2-1] = sum(BDC0[which(shift==0)])/length(which(shift==0))       # in-control identified out-control
      res[i*2] = 1-sum(BDC0[which(shift==1)])/length(which(shift==1))       # out-control identified in-control
    }
    if(index=="CENE"){
      res[i*2] = sum(abs(BDC0-shift))
      res[i*2-1] = ifelse(res[i*2]==0,1,0)
    }
  }
  return(res)
}
simLEB_SD <- function(p, N, n, Delta, OC){
  dat <- GDATA(p, N, n, Delta, OC, MeanCov = TRUE)
  x1 <- dat$xI
  x2 <- dat$xII
  shift = ifelse(dat$delta!=0,1,0)
  #_________________________
  deltaM0 <- deltaMC(x1,x2)
  delta <- deltaM0$delta
  deltaSig <- deltaM0$deltaSig
  #__________  LEB _______________
  LEB <- LEB_diag(N,n,delta,deltaSig,r = 1)
  LEBe <- sum(abs(LEB-shift))
  if(LEBe==0){LEBc = 1}else{LEBc = 0}
  #__________  step-Down _______________
  SD5 <- stepDown_diag(alpha0 = 0.05,delta,deltaSig)
  SD5e <- sum(abs(SD5-shift))
  if(SD5e==0){SD5c = 1}else{SD5c=0}
  
  SD1 <- stepDown_diag(alpha0 = 0.01,delta,deltaSig)
  SD1e <- sum(abs(SD1-shift))
  if(SD1e==0){SD1c = 1}else{SD1c=0}
  
  SD2 <- stepDown_diag(alpha0 = 0.002,delta,deltaSig)
  SD2e <- sum(abs(SD2-shift))
  if(SD2e==0){SD2c = 1}else{SD2c=0}
  #_________________________
  CENE <- c(LEBc,LEBe,SD5c,SD5e,SD1c,SD1e,SD2c,SD2e)
  return(CENE)
}
#***********************************************
runTime <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
CENE.fig <- function(data){
  table <- table[,-1];
  op <- par(oma = c(2.7,1.5,1.5,0.75), mar = c(0.1,1,0.5,0.1)+0.1, mgp = c(1,0,0))
  table[,c('C_S5','C_S1','C_S2','C_LEB','C_B6','C_B7')] = 
    table[,c('C_S5','C_S1','C_S2','C_LEB','C_B6','C_B7')]+1
  table[,c('ENE_S5','ENE_S1','ENE_S2','ENE_LEB','ENE_B6','ENE_B7')] = 
    table[,c('ENE_S5','ENE_S1','ENE_S2','ENE_LEB','ENE_B6','ENE_B7')]/25
  
  plot(0,0,type = 'l', axes=FALSE,xlab = '',ylab = '', xlim=c(0.8,35),ylim = c(0,2))
  xy <- par("usr");xy
  
  points(c(xy[1],5.4),c(xy[4],xy[4]),type='l')
  points(c(xy[1],5.4),c(0.95,0.95),type='l')
  points(c(xy[1],5.4),c(xy[3],xy[3]),type='l')
  
  lineaxis = c(0.6,5.4)
  for(i in 1:5){
    lineaxis = lineaxis + 5
    points(lineaxis,c(xy[4],xy[4]),type='l')
    points(lineaxis,c(0.95,0.95),type='l')
    points(lineaxis,c(xy[3],xy[3]),type='l')
  }
  points(c(xy[2],30.6),c(xy[4],xy[4]),type='l')
  points(c(xy[2],30.6),c(0.95,0.95),type='l')
  points(c(xy[2],30.6),c(xy[3],xy[3]),type='l')
  
  abline(v=c(5,10,15,20,25,30)+0.4,lwd=1,lty = 1)
  abline(v=c(5,10,15,20,25,30)+0.6,lwd=1,lty = 1)
  
  abline(v=xy[1],lwd=1,lty = 1)
  abline(v=xy[2],lwd=1,lty = 1)
  
  abline(h=seq(0,2,by=0.1),lwd=1,col = "lightgray", lty = "dotted")
  
  xyz  = xy[2] - xy[1]
  xyz0 = xy[4] - xy[3]
  text(1:35,rep(xy[3]-xyz0*0.03,35),rep(c("0.50","0.75","1.00","1.25","1.50"),7),xpd = NA,srt = 45)
  for(i in 1:35) points(c(0,0)+i,c(xy[3],xy[3]+0.02),type='l')
  
  # axis(side = 1,at = 1:35, labels = rep(NA,35),las= 0, tck= 0.01)
  
  tick0 = c("0","2.5","5","7.5","10","12.5","15","17.5","20",NA,
            "0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0");
  axis(side = 2,at = seq(0,2,by=0.1), labels = tick0,las= 0, tck= 0.01)
  
  text(xy[1]-xyz*0.04,1.5,"C",xpd = NA, cex = 1.2,adj = c(0.5,0.5),srt = 90)
  text(xy[1]-xyz*0.04,0.5,"ENE",xpd = NA, cex = 1.2,adj = c(0.5,0.5),srt = 90)
  
  rabRow = 1:5
  for(ii in 1:7){
    points(rabRow,table[rabRow,'C_S5'],type='b', pch = 1, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'C_S1'],type='b', pch = 2, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'C_S2'],type='b', pch = 22, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'C_LEB'],type='b', pch = 4, lty = 2, lwd = 2, col = "blue")
    points(rabRow,table[rabRow,'C_B6'],type='b', pch = 15, lty = 1, lwd = 2, col = "black")
    points(rabRow,table[rabRow,'C_B7'],type='b', pch = 16, lty = 1, lwd = 2, col = "black")
    rabRow = rabRow + 5
  }
  rabRow = 1:5
  for(jj in 1:7){
    points(rabRow,table[rabRow,'ENE_S5'],type='b', pch = 1, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'ENE_S1'],type='b', pch = 2, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'ENE_S2'],type='b', pch = 22, lty = 3, lwd = 2, col = "red")
    points(rabRow,table[rabRow,'ENE_LEB'],type='b', pch = 4, lty = 2, lwd = 2, col = "blue")
    points(rabRow,table[rabRow,'ENE_B6'],type='b', pch = 15, lty = 1, lwd = 2, col = "black")
    points(rabRow,table[rabRow,'ENE_B7'],type='b', pch = 16, lty = 1, lwd = 2, col = "black")
    rabRow = rabRow + 5
  }
  text(3,-0.28,expression(Sigma[OC[1]]),xpd = NA, cex = 1.2)
  text(8,-0.28,expression(Sigma[OC[2]]),xpd = NA, cex = 1.2)
  text(13,-0.28,expression(Sigma[OC[3]]),xpd = NA, cex = 1.2)
  text(18,-0.28,expression(Sigma[OC[4]]),xpd = NA, cex = 1.2)
  text(23,-0.28,expression(Sigma[OC[5]]),xpd = NA, cex = 1.2)
  text(28,-0.28,expression(Sigma[OC[6]]),xpd = NA, cex = 1.2)
  text(33,-0.28,expression(Sigma[OC[7]]),xpd = NA, cex = 1.2)
  text(36.5,xy[3]-xyz0*0.02,expression((Delta)),xpd = NA, cex = 1.2)
  
  polygon(c(0,10.4,10.4,0),c(xy[4],xy[4],xy[4]+xyz0*0.05,xy[4]+xyz0*0.05),
          col = "lightgray", border = "black",xpd = NA)
  polygon(c(10.6,25.4,25.4,10.6),c(xy[4],xy[4],xy[4]+xyz0*0.05,xy[4]+xyz0*0.05),
          col = "lightgray", border = "black",xpd = NA)
  polygon(c(25.6,36,36,25.6),c(xy[4],xy[4],xy[4]+xyz0*0.05,xy[4]+xyz0*0.05),
          col = "lightgray", border = "black",xpd = NA)
  text(5.5,xy[4]+xyz0*0.027,"Scenario I",xpd = NA, cex = 1.2)
  text(18,xy[4]+xyz0*0.027,"Scenario II",xpd = NA, cex = 1.2)
  text(30.5,xy[4]+xyz0*0.027,"Scenario III",xpd = NA, cex = 1.2)
  legend("topright", yjust=0.5, inset = .015,
         legend = c("h = 0.6(BDC)", "h = 0.7(BDC)", "LEB",
                    expression(paste( alpha[0]," = 0.05  (S-D)  ")),
                    # "(Step-down)",
                    expression(paste( alpha[0]," = 0.01  (S-D)  ")),
                    # "(Step-down)",
                    expression(paste( alpha[0]," = 0.002(S-D)  "))
                    # "(Step-down)"
         ),ncol = 1,
         col = c("black","black","blue",rep("red",3)), y.intersp = 1, #lty = c(3,3,3,2,1,1), lwd = 1.5,
         pch = c(15,16,4,1,2,22), cex = 0.75, bg = "white", xpd = NA, title = " Methods", title.adj = 0)  
  par(op)
}
