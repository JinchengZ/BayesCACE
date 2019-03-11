
if(!is.null(trace)){
  cat("Start saving MCMC samples for trace plots...\n")
  for (i in 1:length(fullparam)){
    if(is.element(fullparam[i],trace)) {
      for(j in 1:n.chains){
        out[[paste(fullparam[i], "_chain_", j, sep="")]]<-as.vector(jags.out$samples[[j]][,fullparam[i]])
      }
    }
  } 
}



if(is.element("CACE",trace)){
  mcmc
  par(mfrow=c(n.chains,1))
  for(j in 1:n.chains){
    mcmc<-as.vector(jags.out$samples[[j]][,"CACE"])
    plot(temp,type="l",col="red",ylab="CACE",xlab="Iterations",main=paste("Chain",j))
  }
}

if(postdens){
  cat("Start creating posterior density plot for CACE...\n")
  mcmc<-NULL
  for(j in 1:n.chains){
    mcmc<-c(mcmc,as.vector(jags.out$samples[[j]][,"CACE"]))
  }
  out$mcmc <- mcmc
}