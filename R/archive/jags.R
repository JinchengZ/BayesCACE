CACEjags <-
  function(data, 
           param=c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia", "alpha.n", "alpha.a"),
           model="full",prior.type="default", #a=0.001,b=0.001,c=10,
           delta.n=TRUE, delta.a=TRUE, delta.u=TRUE, delta.v=TRUE, delta.s=TRUE, delta.b=TRUE, 
           higher.better=FALSE,digits=4,n.adapt=5000,n.iter=100000,
           n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),
           conv.diag=FALSE,trace=NULL,dic=TRUE,postdens=FALSE,mcmc.samples=FALSE)    {
    ## check the input parameters
    options(warn=1)
    
    if(missing(data)) stop("need to specify data")
    if(!missing(data) ){
      study.id <- data$study.id[complete.cases(data)]
      n000<-data$n000[complete.cases(data)]
      n001<-data$n001[complete.cases(data)]
      n010<-data$n010[complete.cases(data)]
      n011<-data$n011[complete.cases(data)]
      n100<-data$n100[complete.cases(data)]
      n101<-data$n101[complete.cases(data)]
      n110<-data$n110[complete.cases(data)]
      n111<-data$n111[complete.cases(data)]
      cat("NA is not allowed in the input data set;\n")
      cat("the rows containing NA are removed.\n")
    }

    if(length(study.id)!=length(n000) | length(n000)!=length(n001) | length(n001)!=length(n010) | 
       length(n010)!=length(n011) | length(n011)!=length(n100) | length(n100)!=length(n101) |
       length(n101)!=length(n110) | length(n110)!=length(n111) ){
      stop("study.id, n000, n001, n010, n011, n100, n101, n110, and n111 have different lengths.")
    }
    

    ## jags model
    if(model=="full"){
      modelstring<-model.full(prior.type)
    }
    Ind <- rep(1, 6)
    if (!delta.n) Ind[1]=0
    if (!delta.a) Ind[2]=0
    if (!delta.u) Ind[3]=0
    if (!delta.v) Ind[4]=0
    if (!delta.s) Ind[5]=0
    if (!delta.b) Ind[6]=0

    ## jags data
    if(model == "full"){
      if(prior.type == "default"){
        Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
        N0 <- n000+n001+n010+n011
        N1 <- n100+n101+n110+n111
        R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
        I <- length(Ntol)
        r <- (n100+n101+n110+n111)/Ntol
        Ind <- Ind
        data.jags <- list(N0=N0, N1=N1, R=R, I=I, Ind=Ind)
      }
    }

    
    ## jags initial value
    rng.seeds<-sample(1000000,n.chains)
    init.jags <- vector("list", n.chains)
    for(ii in 1:n.chains){
      init.jags[[ii]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
    }
    
    ## parameters to be monitored in jags
    if(!is.element("CACE",param)) param<-c("CACE",param)
    if(!is.null(trace)){
      if(!any(is.element(trace, param))) stop("at least one effect size in argument trace is not specified in argument param.")
    }
    fullparam <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia", "sigma.n", "sigma.a", 
                   "sigma.b", "sigma.v", "sigma.u", "sigma.s", 
                   "delta.n", "delta.a", 
                   "delta.b", "delta.v", "delta.u", "delta.s", 
                   "alpha.n", "alpha.a, rho")
    if(!any(is.element(param, fullparam))) stop("parameters must be specified from the following:  CACE, u1out, v1out, s1out, b1out, 
                   pic, pin, pia, sigma.n, sigma.a, 
                   sigma.b, sigma.v, sigma.u, sigma.s, 
                   delta.n, delta.a, 
                   delta.b, delta.v, delta.u, delta.s, 
                   alpha.n, alpha.a, rho.")
    monitor<-param[is.element(param,fullparam)]
    
    ## run jags
    cat("Start running MCMC...\n")
    jags.m<-jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,
                       n.chains=n.chains,n.adapt=n.adapt)
    update(jags.m,n.iter=n.burnin)
    jags.out<-coda.samples.dic(model=jags.m,variable.names=param,n.iter=n.iter,thin=n.thin)
    smry<-summary(jags.out$samples)
    smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")])
    smry<-signif(smry,digits=digits)
    
    out<-NULL
    out$model<-paste(model, "model", sep=" ")

    #dic
    dev<-jags.out$dic[[1]] # mean deviance
    pen<-jags.out$dic[[2]] # pD
    pen.dev<-dev+pen # DIC
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("D.bar","pD","DIC")
    colnames(dic.stat)<-""
    out$DIC<-dic.stat

    out$CACE<-smry[c("CACE"), ]

    
    if(conv.diag){
      cat("Start calculating MCMC convergence diagnostic statistics...\n")
      conv.out<-gelman.diag(jags.out$samples,multivariate=FALSE)
      out$conv.out<-conv.out$psrf
    }

    
    if(mcmc.samples){
      out$mcmc.samples<-jags.out$samples
    }
    
    if(!is.null(trace)){
      cat("Start creating trace plots...\n")
    }
    
    if(is.element("CACE",trace)){
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(jags.out$samples[[j]][,"CACE"])
        plot(temp,type="l",col="red",ylab="CACE",xlab="Iterations",main=paste("Chain",j))
      }
    }

    if(postdens){
      cat("Start creating posterior density plot for XAXE...\n")
      mcmc<-NULL
      for(j in 1:n.chains){
        mcmc<-c(mcmc,as.vector(jags.out$samples[[j]][,"CACE"]))
      }
      out$mcmc <- mcmc
    }
    class(out)<-"CACEjags"
    return(out)
    options(warn=0)   
}


