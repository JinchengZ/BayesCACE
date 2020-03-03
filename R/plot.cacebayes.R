#' @title plot.cacebayes
#' @importFrom grDevices rainbow 
#' @importFrom graphics par plot
#' @importFrom stats acf bw.SJ density
#' @param obj an S3 objective to plot
#' @param param parameters
#' @param trialnumber indicator for the trial number
#' @return
#' @export

plot.cacebayes <- 
  function(obj, which = c("trace", "density", "autocorr"), 
           param = c("CACE"), trialnumber = 1, ...) {
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
    if(class(obj$mcmc.samples) =="list") {
      if(!trialnumber %in% c(1:length(obj$mcmc.samples)))
        stop("'trialnumber' should be the trial number to be plotted \n")
      else {x <- obj$mcmc.samples[[trialnumber]]}
    }
  
  for (i in 1:length(which)) {
    whichi <- which[i]
    param <- param
    
    n.chains <- length(x)
    cols<-rainbow(n.chains,s=0.6,v=0.6)
    
    if (whichi == "trace") {   
      # if (n.chains==1) {
      #   nams <- dimnames(x)[[2]]
      #   if(is.element(param, nams)) {
      #     temp<-as.vector(x[,param])
      #     plot(temp,type="l",col=cols[1],ylab=param,xlab="Iterations", 
      #          main = paste("Trace Plot of", param))  
      #   } 
      # }
      # else if (n.chains >1) {
      #   nams <- dimnames(x[[1]])[[2]]
      #   if(is.element(param, nams)) {
      #     temp<-as.vector(x[[1]][,param])
      #     plot(temp,type="l",col=cols[1],ylab=param,xlab="Iterations", 
      #          main = paste("Trace Plot of", param))    
      #   } 
      #   for(j in 2:n.chains){
      #     temp<-as.vector(x[[j]][,param])
      #     lines(temp,type="l",col=cols[j])
      #   }
      # }
      
      par(mfrow=c(n.chains,1))
      for(j in 1:n.chains){
        temp<-as.vector(x[[j]][,param])
        plot(temp,type="l",col="red",ylab=param,xlab="Iterations",
            main = paste("Trace Plot of", param,", Chain",j), 
            lwd=1,cex.axis=1.2,cex.lab=1.4)
      }
      par(mfrow=c(1,1))
    }  
    
    else if (whichi == "density") {   
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
    else if (whichi == "autocorr"){
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
  }
  
  invisible()
}

