CACEjags <-
  function(data, 
           param=c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia", "alpha.n", "alpha.a"),
           model="fixed",prior.type="diffuse", #a=0.001,b=0.001,c=10,
           higher.better=FALSE,digits=4,n.adapt=5000,n.iter=100000,
           n.burnin=floor(n.iter/2),n.chains=3,n.thin=max(1,floor((n.iter-n.burnin)/100000)),
           conv.diag=FALSE,trace=NULL,dic=FALSE,postdens=FALSE,mcmc.samples=FALSE)    {
    ## check the input parameters
    options(warn=1)
    
    if(missing(data)) stop("need to specify data")
    if(!missing(data) ){
      study.id <- study.id[complete.cases(data)]
      n000<-n000[complete.cases(data)]
      n001<-n001[complete.cases(data)]
      n010<-n010[complete.cases(data)]
      n011<-n011[complete.cases(data)]
      n100<-n100[complete.cases(data)]
      n101<-n101[complete.cases(data)]
      n110<-n110[complete.cases(data)]
      n111<-n111[complete.cases(data)]
      cat("NA is not allowed in the input data set;\n")
      cat("the rows containing NA are removed.\n")
    }

    if(length(study.id)!=length(n000) | length(n000)!=length(n001) | length(n001)!=length(n010) | 
       length(n010)!=length(n011) | length(n011)!=length(n100) | length(n100)!=length(n101) |
       length(n101)!=length(n110) | length(n110)!=length(n111) ){
      stop("study.id, n000, n001, n010, n011, n100, n101, n110, and n111 have different lengths.")
    }
    

    ## jags model
    if(model=="fixed"){
      modelstring<-model.full(prior.type)
    }

    ## jags data
    
    if(model == "fixed"){
      if(prior.type == "diffuse"){
        Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
        N0 <- n000+n001+n010+n011
        N1 <- n100+n101+n110+n111
        R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
        I <- length(Ntol)
        r <- (n100+n101+n110+n111)/Ntol
        data.jags <- list(N0=N0, N1=N1, R=R, I=I)
      }
    }

    
    ## jags initial value
    rng.seeds<-sample(1000000,n.chains)
    # mu.init<-numeric(ntrt)
    # for(i in 1:ntrt){
    #   mu.init[i]<-sum(event.n[t.id==t.id[i]])/sum(total.n[t.id==t.id[i]])
    # }
    init.jags <- vector("list", n.chains)
    for(ii in 1:n.chains){
      init.jags[[ii]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
    }
    
    
    ## parameters to be monitored in jags
    if(!is.element("CACE",param)) param<-c("CACE",param)
    if(!is.null(trace)){
      if(!any(is.element(trace, param))) stop("at least one effect size in argument trace is not specified in argument param.")
    }
    monitor<-param[!is.element(param,c("CACE", "u1out", "v1out", "s1out", "b1out", 
                                       "pic", "pin", "pia", "alpha.n", "alpha.a"))]

    
    ## run jags
    cat("Start running MCMC...\n")
    jags.m<-jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,
                       n.chains=n.chains,n.adapt=n.adapt)
    update(jags.m,n.iter=n.burnin)
    jags.out<-coda.samples(model=jags.m,variable.names=param,n.iter=n.iter,thin=n.thin)
    smry<-summary(jags.out)
    smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")])
    smry<-signif(smry,digits=digits)
    
    
    # out<-NULL
    # out$model<-"Full model"
    # CACE.id<-grep("CACE",rownames(smry))
    # CACE.stat<-array(paste(format(round(smry[CACE.id,"Mean"],digits=digits),nsmall=digits)," (",format(round(smry[CACE.id,"SD"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
    # colnames(CACE.stat)<-"Mean (SD)"
    # rownames(CACE.stat)<-trtname
    # CACE.quan<-array(paste(format(round(smry[CACE.id,"50%"],digits=digits),nsmall=digits)," (",format(round(smry[CACE.id,"2.5%"],digits=digits),nsmall=digits),", ",format(round(smry[CACE.id,"97.5%"],digits=digits),nsmall=digits),")",sep=""),dim=c(ntrt,1))
    # colnames(CACE.quan)<-"Median (95% CI)"
    # rownames(CACE.quan)<-trtname
    # out$CACE<-list(Mean_SD=noquote(CACE.stat),Median_CI=noquote(CACE.quan))
    # 
    
    if(conv.diag){
      cat("Start calculating MCMC convergence diagnostic statistics...\n")
      conv.out<-gelman.diag(jags.out,multivariate=FALSE)
      conv.out<-conv.out$psrf
      write.table(conv.out,"ConvergenceDiagnostic.txt",row.names=rownames(conv.out),col.names=TRUE)
    }
    
    if(dic){
      cat("Start calculating deviance information criterion statistics...\n")
      dic.out<-dic.samples(model=jags.m,n.iter=n.iter,thin=n.thin)
      dev<-sum(dic.out$deviance)
      pen<-sum(dic.out$penalty)
      pen.dev<-dev+pen
      dic.stat<-rbind(dev,pen,pen.dev)
      rownames(dic.stat)<-c("D.bar","pD","DIC")
      colnames(dic.stat)<-""
      out$DIC<-dic.stat
    }
    
    if(mcmc.samples){
      out$mcmc.samples<-jags.out
    }
    
    if(!is.null(trace)){
      cat("Start saving trace plots...\n")
    }
    
    if(is.element("CACE",trace)){
      
      png("TracePlot_CACE.png",res=600,height=8.5,width=11,units="in")
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(jags.out[[j]][,"CACE"])
        plot(temp,type="l",col="red",ylab="CACE",xlab="Iterations",main=paste("Chain",j))
      }
      dev.off()
    }

    if(postdens){
      cat("Start saving posterior density plot for XAXE...\n")
      mcmc<-NULL
      dens<-matrix(0,1,3)
      colnames(dens)<-c("ymax","xmin","xmax")
      
        temp<-NULL
        for(j in 1:n.chains){
          temp<-c(temp,as.vector(jags.out[[j]][,"CACE"]))
        }
        mcmc<-temp
        tempdens<-density(temp)
        dens<-c(max(tempdens$y),quantile(temp,0.001),quantile(temp,0.999))
       
      ymax<-dens[1]
      xmin<-dens[2]
      xmax<-dens[3]
      cols<-rainbow(1,s=1,v=0.6)
      pdf("CACEPlot.pdf")
      par(mfrow=c(1,1),mar=c(5.5,5.5,2,2)+0.1)
      plot(density(mcmc),xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="CACE",ylab="Density",
           main="",col=cols[1],lty=1,lwd=2,cex.axis=2,cex.lab=2)
 
        lines(density(mcmc),col=cols[1],lty=1,lwd=2)

      # legend("topright",legend=trtname,col=cols,lty=1,lwd=2,cex=1.5)
      dev.off()
    }
    class(out)<-"CACEjags"
    return(out)
    options(warn=0)   }
