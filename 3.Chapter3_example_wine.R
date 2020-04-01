setwd("C:\\XuF\\1.Ph.D\\SPD\\SPD_code\\chapter3")
# section 3.4 a real-dataset  ----
#***********************************************
winedata <- read.csv("winequality-white.csv")
nameC <- names(winedata)
indices <- c("FA","VA","CA","RS","CH","FSD","TSD","DE","PH","SU","AL","quality")

names(winedata) <- indices
head(winedata)
dim(winedata)
LV6 <- winedata[winedata$quality==6,]
dim(LV6)
LV4 <- winedata[winedata$quality==4,]
dim(LV4)
#***********************************************
removeobs <- read.csv("removeobs.csv",header = T)
# removeobs <- read.table("clipboard",header = T)
robs1 <- removeobs[1:100,]
LV6 <- LV6[-robs1[,2],]
robs2 <- removeobs[101:148,]
LV6 <- LV6[-robs2[,2],]
robs3 <- removeobs[149:268,]
LV6 <- LV6[-robs3[,2],]
robs4 <- removeobs[269:328,]
LV6 <- LV6[-robs4[,2],]
robs5 <- removeobs[329:391,]
LV6 <- LV6[-robs5[,2],]
robs6 <- removeobs[392:447,]
LV6 <- LV6[-robs6[,2],]
robs7 <- removeobs[448:471,]
LV6 <- LV6[-robs7[,2],]
robs8 <- removeobs[472:483,]
LV6 <- LV6[-robs8[,2],]
dim(LV6)
#***********************************************
xI <- LV6[,-12]
xII <- LV4[,-12]
p <- dim(xI)[2]
h = 0.8;Nruns = 20000;tau = 10000
setwd("C:\\XuF\\1.Ph.D\\SPD\\SPD_code\\Approach")
source("4.BDC.R")
library(mvtnorm)
RES <- BDC_diag(x0=xI,x=xII,h,Nruns,tau,ifwhitening = TRUE)
round(RES,2)
round(diag(RES),3)
#*******************************************************************
# plot Chart of Probability of Delta Equal to One ----
#*******************************************************************
pdf("chap3_wine.pdf",width = 9,height = 6.15)
Components <- c("FA","VA","CA","RS","CH","FSD","TSD","DE","PH","SU","AL","quality")
Components <- Components[12:1]
# pro1 <- diag(RES)
pro1 = c(0.999, 1.000, 1.000, 0.190, 1.000, 1.000, 0.988, 1.000, 0.552, 0.740, 0.268)
pro1 <- pro1[11:1]
probability <- t(cbind(pro1,1-pro1))
op = par(no.readonly = TRUE)
par(mar=c(3.5,3.5,0.1,5)+0.1,mgp=c(2.2,0.1,0),tck=0.01)
# xy <- par("usr");print(xy)
xy = c(-0.01,  1.00,  0.46, 28.54)
shade = 1:22
angle1 <- rep(c(45,0),11)[shade]
angle2 <- rep(c(135,0),11)[shade]
dens <- rep(c(20,0),11)[shade]
xtick <- c('0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0')
barplot(probability,horiz=TRUE,beside = FALSE,space = 1.5,names.arg = Components[-1],axis.lty = 1,
        xlab = list(expression(paste('Probability of ', tilde(eta))),cex=1.2), las=2,
        ylab = list('Components',cex=1.2),axes =F,
        angle = angle1, density = dens, col = "blue",
        legend.text = c(expression(paste(tilde(eta),'=1  ')),expression(paste(tilde(eta),'=0  '))),
        args.legend=list(x=xy[2]+0.15,y= xy[4],bty='n', cex=1.2),cex.axis = 1.2,cex = 1.2)
barplot(probability,horiz=TRUE,beside = FALSE,space = 1.5,names.arg = Components[-1],axis.lty = 1,
        xlab = list(expression(paste('Probability of ', tilde(eta))),cex=1.2), las=2,
        ylab = list('Components',cex=1.2),axes =F,
        angle = angle2, density = dens, col = "blue",
        legend.text = c(expression(paste(tilde(eta),'=1  ')),expression(paste(tilde(eta),'=0  '))),
        args.legend=list(x=xy[2]+0.15,y= xy[4],bty='n', cex=1.2),cex.axis = 1.2,cex = 1.2,add =TRUE)
axis(side=1,at=seq(0,1,by=0.1),labels=xtick,cex.axis=1.2)
grid()
graphics::box()
dev.off()
#*******************************************************************
# plot point-and-figure chart of standardized observations ----
#*******************************************************************
N <- dim(xI)[1];N
n <- dim(xII)[1];n
p <- dim(xII)[2];p
# plot sactter
muI <- apply(xI,2,mean)
SigmaI <- cov(xI)
# eigenvalue and eigenvector
eigenSigmaI <- eigen(SigmaI)  
Lambda <- eigenSigmaI$values
TV <- eigenSigmaI$vectors
# standardizied phase I data
x0 = t(apply(xI,1,function(xI) TV%*%solve(diag(sqrt(Lambda)))%*%t(TV)%*%(xI-muI)))
barx0 <-  round(apply(x0,2,mean),3)
hatSigma0 <- round(cov(x0),3)
invhatSigma0 <- solve(hatSigma0)
# standardizied phase II data
muII <- apply(xII,2,mean)
x = t(apply(xII,1,function(xII) TV%*%solve(diag(sqrt(Lambda)))%*%t(TV)%*%(xII-muII)))
hatSigma <- round(cov(x),1)
barx <- round(apply(x,2,mean),3)
#***********************************************
dim(x0);dim(x)
nameC <- c("Fixed Acidity", "Volatile Acidity","Citric Acid","Residual Sugar",
           "Chlorides","Free Sulfur Dioxide","Total Sulfur Dioxide","Density",
           "pH","Sulphates","Alcohol")  

for(i in 1:11){
  components = i
  name2 <- c("FA","VA","CA","RS","CH","FSD","TSD","DE","PH","SU","AL")
  pdf(paste0("chap3_",name2[components],".pdf"),width = 7.5,height = 6.755)
  par(mar=c(2.5,2.5,2,0.1)+0.1,mgp=c(1.5,0.1,0),tck=0.01)
  phaseI <- x0[1553:1715,components]
  phaseII <- x[,components]
  plot(c(phaseI,phaseII),pch=c(rep(1,163),rep(8,163)),xlab="Observations",type="o",col="blue",
       ylab=nameC[components],cex=1.2,cex.axis=1.2,cex.lab=1.2)
  abline(v=163,lwd=1,col="red",lty=2)
  # plot(c(LV6[,1],LV4[,1]))
  # abline(v=2198,lwd=0.1,col="blue",lty=2)
  Imusd <- round(c(mean(x0[,components]),var(x0[,components])),2)
  IImusd <- round(c(mean(x[,components]),var(x[,components])),2)
  xy <- par("usr")
  legend(x=0,y=xy[4]+yinch(0.5),xpd = TRUE,bty = 'n',cex = 1.2,col="blue",
         legend = substitute(paste("phase I: ",mu,"=",mu1,", ",sigma^2,"=",sigma1),
                             list(mu1=Imusd[1],sigma1=format(Imusd[2],nsmall=2))),pch=1)
  legend(x=170,y=xy[4]+yinch(0.5),xpd = TRUE,bty = 'n',cex = 1.2,col="blue",
         legend = substitute(paste("phase II: ",mu,"=",mu1,", ",sigma^2,"=",sigma1),
                             list(mu1=IImusd[1],sigma1=format(IImusd[2],nsmall=2))),pch=8)
  grid()
  dev.off()
}

# legend(x=0,y=xy[4]+yinch(0.5),xpd = TRUE,bty = 'n',cex = 1.2,
#        legend = expression(paste("phase I: ",mu,"=","   , ",sigma^2,"=          ")),pch=1)
# legend(x=55,y=xy[4]+yinch(0.45),xpd = TRUE,bty = 'n',cex = 1.2,
#        legend = paste0(Imusd[1],"         ",format(Imusd[2],nsmall=2)))
# 
# legend(x=170,y=xy[4]+yinch(0.5),xpd = TRUE,bty = 'n',cex = 1.2,
#        legend = expression(paste("phase II: ",mu,"=","   , ",sigma^2,"=            ")),pch=8)
# legend(x=228,y=xy[4]+yinch(0.45),xpd = TRUE,bty = 'n',cex = 1.2,
#        legend = paste0(IImusd[1],"         ",format(IImusd[2],nsmall=2)))
# par(op)
