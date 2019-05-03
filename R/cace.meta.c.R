cace.meta.c <-
  function(data, 
           param = c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia"),
           prior.type = "default", 
           delta.n = TRUE, delta.a = TRUE, delta.u = TRUE, delta.v = TRUE, 
           delta.s = TRUE, delta.b = TRUE, cor = TRUE, 
           digits = 3, n.adapt = 1000, n.iter = 100000,
           n.burnin = floor(n.iter/2), n.chains = 3, n.thin = max(1,floor((n.iter-n.burnin)/100000)),
           conv.diag = FALSE, mcmc.samples = FALSE, study.specific = FALSE)    {
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
       length(n101)!=length(n110) | length(n110)!=length(n111) )
      stop("study.id, n000, n001, n010, n011, n100, n101, n110, and n111 have different lengths. \n")
    
    if ((!(delta.n & delta.a)) & cor){
      warning("'cor' can be assigned as TRUE only if both delta.n and delta.a are TRUE.\n
              the model is continued by forcing delta.n=TRUE and delta.a=TRUE")
      delta.n=TRUE
      delta.a=TRUE
    }
    
    if ((!(delta.u|delta.v)) & study.specific){
      warning("no random effect is assigned to the response rate u1 or v1, \n
           study-specific CACE is the same across studies. \n
           the model is continued by making 'study.specific=FALSE'. \n
           to make a CACE forestplot, please run 'cace.study' to estimate study level CACEs. \n")
      study.specific=FALSE
    }
      
    
    Ind <- rep(1, 7)
    if (!delta.n) Ind[1]=0
    if (!delta.a) Ind[2]=0
    if (!delta.u) Ind[3]=0
    if (!delta.v) Ind[4]=0
    if (!delta.s) Ind[5]=0
    if (!delta.b) Ind[6]=0
    if (!cor) Ind[7]=0
    
    ## jags model
    modelstring<-model.meta.c(prior.type, Ind)

    ## jags data
    if(prior.type == "default"){
      Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
      N0 <- n000+n001+n010+n011
      N1 <- n100+n101+n110+n111
      R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
      I <- length(Ntol)
      pi <- pi
      data.jags <- list(N0=N0, N1=N1, R=R, I=I, Ind=Ind, pi=pi)
    }

    
    ## jags initial value
    rng.seeds<-sample(1000000,n.chains)
    init.jags <- vector("list", n.chains)
    for(ii in 1:n.chains){
      init.jags[[ii]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
    }
    
    ## parameters to be paramed in jags
    if(!is.element("CACE",param)) param<-c("CACE",param)

    fullparam <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia", "alpha.n", "alpha.a", 
                   "alpha.u", "alpha.v", "alpha.s", "alpha.b", 
                   "sigma.n", "sigma.a", "rho", "Sigma.rho",
                   "sigma.u", "sigma.v", "sigma.s", "sigma.b")
    if(!any(is.element(param, fullparam))) stop("parameters must be specified from the following:  
                  CACE, u1out, v1out, s1out, b1out, 
                   pic, pin, pia, alpha.n, alpha.a, 
                alpha.u, alpha.v, alpha.s, alpha.b, 
                sigma.n, sigma.a, rho, Sigma.rho,
                sigma.u, sigma.v, sigma.s, sigma.b")
    if(study.specific) param<-c("cacei",param)
    
    ## run jags
    cat("Start running MCMC...\n")
    jags.m<-jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,
                       n.chains=n.chains,n.adapt=n.adapt)
    update(jags.m,n.iter=n.burnin)
    jags.out<-coda.samples.dic(model=jags.m,variable.names=param,n.iter=n.iter,thin=n.thin)
    
    out<-NULL
    out$Ind <- Ind
    out$model<-"cace.meta.c"
    
    smry<-summary(jags.out$samples)
    smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")], 
                smry$statistics[,c("Naive SE","Time-series SE")])
    smry<-signif(smry,digits=digits)
    out$smry <- smry

    #dic
    dev<-jags.out$dic[[1]] # mean deviance
    pen<-jags.out$dic[[2]] # pD
    pen.dev<-dev+pen # DIC
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("D.bar","pD","DIC")
    colnames(dic.stat)<-""
    out$DIC<-dic.stat

    for (i in 1:length(fullparam)){
      if(is.element(fullparam[i],param)) {out[[fullparam[i] ]]<-smry[c(fullparam[i]), ]}
    } 
    
    if(conv.diag){
      cat("MCMC convergence diagnostic statistics are calculated and saved in conv.out\n")
      conv.out<-gelman.diag(jags.out$samples,multivariate=FALSE)
      out$conv.out<-conv.out$psrf
    }
    
    if(mcmc.samples){
      out$mcmc.samples<-jags.out$samples
    }

    class(out)<-"cace.Bayes"
    return(out)
    options(warn=0)   
}


