rm(list=ls())
options(digits = 3, scipen=200, max.print=1000) 
library(mvtnorm)
library(MCMCpack)
library(doParallel)

setwd("E:\\XuF\\XuF\\chapter3")
source("0.Estimate_Covariance_of_Delta.R")
source("1.LEB.R")
source("2.Step_down.R")
source("4.BDC.R")

setwd("E:\\XuF\\XuF\\chapter3")
source("3.Chapter3_function.R")

#***********************************************
# section 4.1 Effect of h ----
#***********************************************

N = 200;n = 200;p = 6;ns = 10;h = c(0.28,0.35,0.45);
Delta = seq(0.5,1.5,by=0.25);
hNum = length(h)
resColnames <- c('Delta_0.00',paste0("h",rep(1:hNum,rep(2,hNum)),c("_I","_II")))
write.table(t(resColnames),"N200n200_e.txt",quote = FALSE,col.names = FALSE,row.names=FALSE)
Delta0 = c('0.50','0.75','1.00','1.25','1.50')
# detectCores()       
cl = makeCluster(5) 
registerDoParallel(cl)
startTime <- Sys.time()

# for(OC in 1:7){
#   for(jj in 1:5){
    OC = 1;jj =  5
    Res <- foreach(p = rep(p,ns), .combine='rbind',.multicombine = TRUE,.maxcombine=ns,
                   .inorder=FALSE,.packages=c("mvtnorm",'MCMCpack','doParallel')) %dopar%
      simBDC(h = h, p = p, N = N, n = n, Delta = Delta[jj], OC = OC,index = "ErrorRate")
    ErrorRates0 = apply(Res,2,sum)/ns
    names(ErrorRates0) <- paste0("h",rep(1:hNum,rep(2,hNum)),c("_I","_II"))
    write.table(t(round(ErrorRates0,2)),"N200n200_e.txt",append = TRUE,quote = FALSE,
                col.names = FALSE,row.names = paste0('Delta',OC,'_',Delta0[jj]))
#   }
# }

runTime(startTime)
stopCluster(cl)

#***********************************************
# section 4.2 Comparison with LEB and Step-down procedures ----
#***********************************************
p = 5;n = 100;N = 2*n;
Delta = seq(0.5,1.5,by=0.25);
nsBDC = 200;nsLEB=10^4
resColnames <- c('Delta_0.00','LEBC','LEBe','SD5C','SD5e','SD1C','SD1e','SD2C','SD2e',
                 'BD6C','BD6e','BD7C','BD7e')
# write.table(t(resColnames),"N300n150p5_c.txt",quote = FALSE,col.names = FALSE,row.names=FALSE)
# Delta0 = c('0.4','0.8','1.2','1.6','2.0')
Delta0 = c('0.50','0.75','1.00','1.25','1.50')
cl = makeCluster(5) 
registerDoParallel(cl)
startTime <- Sys.time()

for(OC in 1:7){
  for(jj in 1:5){
    BDC <- foreach(p = rep(p,nsBDC), .combine='rbind',.multicombine = TRUE,.maxcombine=nsBDC,
                   .inorder=FALSE,.packages=c("mvtnorm",'MCMCpack','doParallel')) %dopar%
      simBDC(h = c(0.6,0.7), p = p, N = N, n = n, Delta = Delta[jj], OC = OC,index = "CENE")
    CENE_BDC <- apply(BDC,2,sum)/nsBDC
    # 
    LEB_Stepdown <- foreach(p = rep(p,nsLEB), .combine='rbind',.multicombine = TRUE,.maxcombine=nsLEB,
                            .inorder=FALSE,.packages=c("mvtnorm",'MCMCpack','doParallel')) %dopar%
      simLEB_SD(p,N,n,Delta = Delta[jj],OC = OC)
    CENE_LEB_Stepdown <- apply(LEB_Stepdown,2,sum)/(nsLEB)
    
    CENE <- c(CENE_LEB_Stepdown,CENE_BDC)
    write.table(t(round(CENE,8)),"N200n100p5_c.txt",append = TRUE,quote = FALSE,
                col.names = FALSE,row.names = paste0('Delta',OC,'_',Delta0[jj]))
  }
}

stopCluster(cl)
runTime(startTime)


#***********************************************
#  draw figure ----
#***********************************************
library(openxlsx)
table <- read.xlsx("Table1to6.xlsx",sheet = "Tab6", startRow = 3);
setEPS()
postscript(paste0("Table_6.eps"),width = 12, height = 5.98)
CENE.fig(data=table)
dev.off()


