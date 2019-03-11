plot.cacebayes <- 
  function(obj, which = c("trace", "density", "autocorr"), 
           param = c("CACE"),  ...) {
  if(missing(obj)) stop("need to specify obj, the object generated from one of the following functions:\n
         cace.study, cace.meta.c, cace.meta.ic.\n")
  if (!inherits(obj, "cace.Bayes"))
    stop("Use only with 'cace.Bayes' objects, the output generated from one of the following functions:\n
         cace.study, cace.meta.c, cace.meta.ic.\n")
  if(is.null(obj$mcmc.samples))
    stop("'obj$mcmc.samples' is missing. please set mcmc.samples=TRUE when running the function
         cace.study, cace.meta.c or cace.meta.ic.")
  if(!class(obj$mcmc.samples) %in% c("mcmc.list","list") )
    stop("'obj$mcmc.samples' must be a mcmc.list or list generated from one of the following functions:\n
         cace.study, cace.meta.c, cace.meta.ic.")
    if(class(obj$mcmc.samples) =="mcmc.list") {x <- obj$mcmc.samples}
    else if(class(obj$mcmc.samples) =="list") {x <- obj$mcmc.samples[[1]]}
    
  which <- match.arg(which)  
  param <- match.arg(param)
  
  n.chains <- length(x)
  cols<-rainbow(n.chains,s=0.6,v=0.6)
    
  if (which == "trace") {   
    if (n.chains==1) {
      nams <- dimnames(x)[[2]]
      if(is.element(param, nams)) {
        temp<-as.vector(x[,param])
        plot(temp,type="l",col=cols[1],ylab=param,xlab="Iterations", 
             main = paste("Trace Plot of", param))  
      } 
    }
    else if (n.chains >1) {
      nams <- dimnames(x[[1]])[[2]]
      if(is.element(param, nams)) {
        temp<-as.vector(x[[1]][,param])
        plot(temp,type="l",col=cols[1],ylab=param,xlab="Iterations", 
             main = paste("Trace Plot of", param))    
      } 
      for(j in 2:n.chains){
        temp<-as.vector(x[[j]][,param])
        lines(temp,type="l",col=cols[j])
      }
    }
  }  

  else if (which == "density") {   
    if (n.chains==1) {
      nams <- dimnames(x)[[2]]
      if(is.element(param, nams)) {
        temp<-as.vector(x[,param])
        bw <- bw.SJ(temp) * 1.5
        plot(density(temp, bw = bw), xlab = param, 
             main = paste("Density of", param))
      } 
    }
    else if (n.chains >1) {
      nams <- dimnames(x[[1]])[[2]]
      if(is.element(param, nams)) {
        temp<-NULL
        for(j in 1:n.chains){
          temp<-c(temp, as.vector(x[[j]][,param]))
        }
        bw <- bw.SJ(temp) * 1.5
        plot(density(temp, bw = bw), xlab = param, 
             main = paste("Density of", param))
      } 
    }
  }  
  else if (which == "autocorr"){
    if (n.chains==1) {
      nams <- dimnames(x)[[2]]
      if(is.element(param, nams)) {
        temp<-as.vector(x[,param])
        acf(temp, ylab = param, main = paste("Series", param))
      } 
    }
    else if (n.chains >1) {
      nams <- dimnames(x[[1]])[[2]]
      if(is.element(param, nams)) {
        temp<-NULL
        for(j in 1:n.chains){
          temp<-c(temp, as.vector(x[[j]][,param]))
        }
        acf(temp, ylab = param, main = paste("Series", param))
      } 
    }
  }
  invisible()
}

