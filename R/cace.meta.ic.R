#' @export
#' 
cace.meta.ic <-
  function(data, param = c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia"),
           prior.type = "default", 
           delta.n = TRUE, delta.a = TRUE, delta.u = TRUE, delta.v = TRUE, 
           delta.s = TRUE, delta.b = TRUE, cor = TRUE, 
           digits = 3, n.adapt = 1000, n.iter = 100000,
           n.burnin = floor(n.iter/2), n.chains = 3, 
           n.thin = max(1,floor((n.iter-n.burnin)/100000)),
           conv.diag = FALSE, mcmc.samples = FALSE, study.specific = FALSE)    {
    ## check the input parameters
    options(warn=1)
    
    if(missing(data)) stop("need to specify data")
    if(!missing(data) ){
      data$miss.r0 <- ifelse((data$n000==0 & data$n001==0 & data$n010==0 & data$n011==0), 1, 0)
      data$miss.r1 <- ifelse((data$n100==0 & data$n101==0 & data$n110==0 & data$n111==0), 1, 0)
      data$miss <- ifelse((data$miss.r0==1|data$miss.r1==1), 1, 0)
      temp <- data[order(data$miss.r0, data$miss.r1),]
      study.id <- temp$study.id[complete.cases(temp)]
      n000<-temp$n000[complete.cases(temp)]
      n001<-temp$n001[complete.cases(temp)]
      n010<-temp$n010[complete.cases(temp)]
      n011<-temp$n011[complete.cases(temp)]
      n100<-temp$n100[complete.cases(temp)]
      n101<-temp$n101[complete.cases(temp)]
      n110<-temp$n110[complete.cases(temp)]
      n111<-temp$n111[complete.cases(temp)]
      n0s0<-temp$n0s0[complete.cases(temp)]
      n0s1<-temp$n0s1[complete.cases(temp)]
      n1s0<-temp$n1s0[complete.cases(temp)]
      n1s1<-temp$n1s1[complete.cases(temp)]
      miss.r0<-temp$miss.r0[complete.cases(temp)]
      miss.r1<-temp$miss.r1[complete.cases(temp)]
      miss<-temp$miss[complete.cases(temp)]
      cat("NA is not allowed in the input data set;\n")
      cat("the rows containing NA are removed.\n")
    }
    
    if(length(study.id)!=length(n000) | length(n000)!=length(n001) | length(n001)!=length(n010) | 
       length(n010)!=length(n011) | length(n011)!=length(n100) | length(n100)!=length(n101) |
       length(n101)!=length(n110) | length(n110)!=length(n111) | length(n111)!=length(n0s0) | 
       length(n0s0)!=length(n0s1) | length(n0s1)!=length(n1s0) | length(n1s0)!=length(n1s1) )
      stop("study.id, n000, n001, n010, n011, n100, n101, n110, n111, 
           n0s0, n0s1, n1s0, and n1s1 have different lengths. \n")
    
    if ((!(delta.n & delta.a)) & cor){
      warning("'cor' can be assigned as TRUE only if both delta.n and delta.a are TRUE.\n
              the model is continued by forcing delta.n=TRUE and delta.a=TRUE")
      delta.n=TRUE
      delta.a=TRUE
    }
    
    if (!(delta.u|delta.v)){
      warning("no random effect is assigned to the response rate u1 or v1, \n
              study-specific CACE is the same across studies. \n
              a CACE forestplot cannot be made. \n")
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
    modelstring<-model.meta.ic(prior.type, Ind)
    
    ## jags data
    if(prior.type == "default"){
    
    n1 <- sum(miss.r0==0 & miss.r1==0) #  4+4
    n2 <- sum(miss.r0==0 & miss.r1==1) #  4+2
    n3 <- sum(miss.r0==1 & miss.r1==0) #  2+4
    n4 <- sum(miss.r0==1 & miss.r1==1) #  4+4
    n <- length(study.id)
    
    N0_4 <- n000+n001+n010+n011
    N0_2 <- n0s1+n0s0
    N1_4 <- n100+n101+n110+n111
    N1_2 <- n1s1+n1s0
    
    R0_4 <- cbind(n000,n001,n010,n011)
    R0_2 <- cbind(n0s0,n0s1)
    R1_4 <- cbind(n100,n101,n110,n111)
    R1_2 <- cbind(n1s0,n1s1)
    
    R1 <- cbind(R0_4, R1_4)
    R2 <- cbind(R0_4, R1_2)
    R3 <- cbind(R0_2, R1_4)
    R4 <- cbind(R0_2, R1_2)
    pi <- pi
    
    data.jags <- list(N0_4=N0_4, N0_2=N0_2, N1_4=N1_4, N1_2=N1_2, 
                 R0_4=R0_4, R0_2=R0_2, R1_4=R1_4, R1_2=R1_2, 
                 n1=n1, n2=n2, n3=n3, n4=n4, Ind=Ind, pi=pi)
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
    out$model<-"cace.meta.ic"
    
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


