#################### NO-Oultier Case
rm(list=ls())
setwd("/Users/zhenzhen/Desktop")
#install.packages("R2jags")
#install.packages("arm")
#install.packages("BRugs") 
#install.packages("R2OpenBUGS")
#library(MSBVAR)
#library(VGAM)
#library(lars)
#library(glmnet)
#library(R2WinBUGS)
sgn<-function(a){if(a>=0) 1 else -1}

#*******************************************************************************
#************* Set up the true values
n=50
p=20#   p=8#    p=50#   
simulation=50
beta.true<-rep(0,times=p)
#beta.true[1:8]<-c(3,1.5,0,0,2,0,0,0)
#sim00: large cluster, small values
#beta.true[1:12]<-rep(2,times=12)
#sim0: large cluster, large values
#beta.true[1:12]<-rep(-10,times=12)
#sim1: large cluster and symmetric
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,10,10,10,10,10,10)
#sim2: large cluster and non-symmetric
beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,4,4,4,4,4,4)
#sim3: large and small cluster and symmetric
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,10,10)
#sim4: large and small cluster and non-symmetric
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,4,4)
#sim5: small cluster and symmetric
#beta.true[1:12]<-c(-10,-10,-4,-4,-2,-2,2,2,4,4,10,10)
#sim 5a:
#beta.true[1:8]<-c(-10,-10,-2,-2,2,2,10,10)
#sim 5b: 
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,2,2,2,2,2,2)
#sim 5c:
#beta.true[1:12]<-c(-10,-10,-8,-8,-4,-4,4,4,8,8,10,10)
#sim 5d:
#beta.true[1:12]<-c(-8,-8,-6,-6,-4.5,-4.5,4,4,5,5,7,7)
#sim 5e:
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,50,50,50,50,50,50)
#sim 5f:
#beta.true[1:12]<-c(-10,-10,-10,-10,-10,-10,20,20,20,20,20,20)
#sim 5g:
#beta.true[1:12]<-c(-10,-10,-8,-8,-4,-4,12,12,15,15,20,20)
#sim6: small cluster and non-symmetric
#beta.true[1:12]<-c(-10,-10,-4,-4,-2,-2,3,3,6,6,8,8)
#sim6a: small cluster and non-symmetric, 
#beta.true[1:12]<-c(-6,-6,-4,-4,-2,-2,3,3,5,5,7,7)
#sim7: small cluster and non-non symmetric
#beta.true[1:12]<-c(-10,-10,-8,-8,-6,-6,-4,-4,-2,-2,3,3)
#sim8: small cluster and non-non symmetric
#beta.true[1:12]<-c(-10,-10,-8,-8,-6,-6,-4,-4,-2,-2,2,2)
#---------------------------------------------------------
#sim9:
#beta.true[1:12]<-c(-10,-10,-8,-8,-4,-4,12,12,12,15,15,15)
#sim9a:
#beta.true[1:5]<-c(0.6,1.2,1.8,2.4,3.0)
#----------------------------------------------------------
#sim regular:
#beta.true[1:8]<-c(3,1.5,0,0,2,0,0,0)

#sim10: significant values, non-symmmetric nozeros, un-equally spaced, large # of clusters, concentrated values
#beta.true[1:12]<-c(-10,-10,-15,12,16,16,16,18,20,20,20,20)
                                                            
sigma2.true<-1 
#eta.true=1.8# 
eta.true=0.5#       eta.true=1#     
sigmax<-array(0, dim=c(p,p))
#sigmax[row(sigmax)==col(sigmax)]<-1
for(i in 1:p){
  for(j in 1:p){
    sigmax[i,j]<-0.5^abs(i-j)
  }
}


################################################################################
################################################################################
# copies of simulations 
mode.parameter<-array(0,dim=c(simulation,p+2))
median.parameter<-array(0,dim=c(simulation, p+2))
tmse.mode<-rep(0,times=simulation)
tmse.beta<-rep(0,times=simulation)
logl<-rep(0,times=simulation)
logl.true<-rep(0,times=simulation)
logl.ratio<-rep(0,times=simulation)

set.seed(100)
for(sim in 1:simulation){
  cat("Simulation=", sim, "\n")
  #*****************************************************************************
  #Create the dataset
  #x.raw<-rmultnorm(n, mux, vmat=sigmax)
  x.raw<- matrix(rnorm(n*p),n,p) %*% t(chol(sigmax))
  #x<-scale(x.raw, center=TRUE, scale=FALSE) 
  x<-scale(x.raw, center=TRUE, scale=TRUE) 
  med<-rep(0,times=n)
  e<-rep(0,times=n)   
  gmedian<-rep(0,times=n)
  gy<-rep(0,times=n)
  y<-rep(0,times=n)
  for(i in 1:n){
    med[i]<-x[i,]%*%beta.true
    gmedian[i]<-(sgn(med[i])*abs(med[i])^eta.true-1)/eta.true
    e[i]<-rnorm(1,mean=gmedian[i],sd=sqrt(sigma2.true))
    gy[i]<-e[i]   
    y[i]<-sgn(eta.true*gy[i]+1)*abs(gy[i]*eta.true+1)^(1/eta.true)
  }
   
  beta.ols<-lm(y~x)$coefficients
  tau.l2.mean<-1/var(beta.ols)
  #tau.l2.mean<-2
  #*****************************************************************************
  #***************** JAGS loops to estimate parameters
  library(rjags)
  inits=list(gbeta.l2=rep(1,times=p), betau=runif(p), scale.eta=0.2, i_sigma2=0.3)
  data=list('n'=n, 'p'=p,'y'=y, 'x'=x, 'pi'=pi,'tau.l2.mean'=tau.l2.mean)
  parameter=c('beta.lasso', 'eta', 'sigma2')
  jags<-jags.model('pgaussian_tbs_jags.r', data, inits, n.chains = 1, n.adapt = 5000)
  coda<-coda.samples(jags, variable.names=parameter, n.iter=1000)
  latent<-as.matrix(coda[[1]])
  coda.matrix<-array(0,dim=c(nrow(latent),(p+2)))
  coda.matrix<-latent[,1:(p+2)]
  
  
  
    
  #Find the mode of beta parameter  
  #Find the median of beta parameter
  for(j in 1:(p+2)){
    median.parameter[sim,j]<-median(coda.matrix[,j])
  }
  tmse.beta[sim]<-sum((median.parameter[sim,1:p]-beta.true)^2)
  gxbeta<-rep(0,n)
  for(i in 1:n){
    gxbeta[i]<-(sgn(x[i,]%*%median.parameter[sim,1:p])*abs(x[i,]%*%median.parameter[sim,1:p])^eta.true-1)/eta.true
  }
  logl[sim]<- -sum((gy-gxbeta)^2)/(2*sigma2.true)-0.5*n*log(2*pi*sigma2.true)+(eta.true-1)*sum(log(abs(y)))  
  logl.true[sim]<- -sum((gy-gmedian)^2)/(2*sigma2.true)-0.5*n*log(2*pi*sigma2.true)+(eta.true-1)*sum(log(abs(y)))  
  logl.ratio[sim]<-logl[sim]/logl.true[sim]-1
  
}
################################################################################
################################################################################
# estimate of parameters 
round(colMeans(median.parameter),digits=3)     
#calculate total mse
round(median(tmse.beta),digits=5)
round(sd(tmse.beta),digits=5)
summary(logl.ratio)


#*******************************************************************************
#*********** Mask, Swamp, JD
nzerosindex<-which(abs(beta.true)>0)
zerosindex<-which(abs(beta.true)==0)
nzeros.median<-rep(0,times=simulation)
mask.median<-rep(0,times=simulation)
swamp.median<-rep(0,times=simulation)

for(sim in 1:simulation){
  nzeros.median[sim]<-sum(abs(median.parameter[sim,1:p])>0)
  mask.median[sim]<-sum(median.parameter[sim, nzerosindex]==0)/sum(abs(beta.true)>0)
  swamp.median[sim]<-sum(abs(median.parameter[sim, zerosindex])>0)/sum(beta.true==0)
}
jd.median<-sum(mask.median==0)/simulation

mean(nzeros.median)
mean(mask.median)
mean(swamp.median)
jd.median




pdf('coda_pgaussian.pdf')
plot(coda)
dev.off()




