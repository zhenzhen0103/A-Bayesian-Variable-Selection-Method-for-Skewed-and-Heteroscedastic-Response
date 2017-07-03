data {
  C <- 1000000    # this just has to be large enough to ensure all phi's < 1
                  # Note that JAGS will NOT WARN you if C is too small, but will sample from a truncated distribution
  for (i in 1:n) {
      ones[i] <- 1
      #zeros[i]<-0
  }
}

model{
  for(i in 1:n){
    gy[i]<- (step(y[i])*pow(abs(y[i]),eta)-(1-step(y[i]))*pow(abs(y[i]),eta)-1)/eta
#    gy[i]<-y[i]
    med[i]<-inprod(beta.lasso[1:p],x[i,1:p])
    gmedian[i]<-(step(med[i])*pow(abs(med[i]),eta)-(1-step(med[i]))*pow(abs(med[i]),eta)-1)/eta
#    gmedian[i]<-med[i]
    phi[i]<-(gy[i]-gmedian[i])*(gy[i]-gmedian[i])/(2*sigma2)+0.5*log(2*pi*sigma2)-(eta-1)*log(abs(y[i]))
    pones[i]<-exp(-phi[i])/C
    #pones[i]<- dnorm((gy[i]-gmedian[i]), 0, i_sigma2)*pow(abs(y[i]),(eta-1))/C
    ones[i]~dbern(pones[i])
  }
  #prior for beta
  for(j in 1:p){
    beta.lasso[j]<- beta.l2[j]*step(betau[j]-betau.pi)
    betau[j]~dunif(0,1)
   beta.l2[j]<-step(eta*gbeta.l2[j]+1)*pow(abs(eta*gbeta.l2[j]+1),1/eta)-(1-step(eta*gbeta.l2[j]+1))*pow(abs(eta*gbeta.l2[j]+1),1/eta)
#    beta.l2[j]<-gbeta.l2[j]
    gbeta.l2[j]~dnorm(0,tau.l2)
  }
  
  # the grid is the distribution probability for point 0, negative gamma and positive gamma
  betau.pi~dbeta(1.5,1.5) #betau.pi~dbeta(1,1.5) #
  tau.l2~dgamma(tau.l2.mean,1)
  sigma2<-1/i_sigma2
  i_sigma2~dgamma(1.5,1.5)
  eta<-2*scale.eta
  scale.eta~dbeta(1,1)
}
