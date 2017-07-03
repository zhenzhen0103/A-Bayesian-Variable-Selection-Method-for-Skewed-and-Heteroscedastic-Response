
rm(list=ls())
install.packages("R2jags")

sgn<-function(a){if(a>=0) 1 else -1}

#*******************************************************************************
#************* Set up the true values
n=50
p=20
simulation=50
beta.true<-rep(0,times=p)
beta.true[1:8]<-c(3,1.5,0,0,2,0,0,0)
                      
sigma2.true<-1 
#eta.true=1.8# 
eta.true=0.5#      
sigmax<-array(0, dim=c(p,p))
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

