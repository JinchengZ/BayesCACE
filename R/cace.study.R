#' This function performs CACE analysis for a single study using the 
#' likelihood and model specified in the paper Section 2.1, or a two-step 
#' approach for meta-analysis with complete compliance information as 
#' described in the paper Section 2.2 "The two-step approach".
#' @importFrom stats update complete.cases
#' @import rjags
#' @import coda
#' @import metafor
#' @import Rdpack
#' @export
#' @title CACE analysis for a single study, or a two-step approach for meta-analysis 
#' with complete complice information
#' @param data a input dataset the same structure as the example data `epidural_c`, 
#' containing either one row of observations for a single study, or multiple rows referring 
#' to multiple studies in a meta-analysis. 
#' @param param a character string vector indicating the parameters to be tracked and estimated. 
#' By default all parameters in the model (see \code{details}) are included: \eqn{\theta^{CACE}} 
#' (\code{CACE}), \eqn{u_1} (\code{u1}), \eqn{v_1} (\code{v1}), \eqn{s_1} (\code{s1}), \eqn{b_1} (\code{b1}), 
#' \eqn{\pi_a} (\code{pi.a}), \eqn{\pi_n} (\code{pi.n}), and \eqn{\pi_c=1-\pi_a-\pi_n} (\code{pi.c}). 
#' Users can modify the string vector to only include parameters of interest besides 
#' \eqn{\theta^{CACE}}. 
#' @param prior.type the default priors are used by the default assignment `prior.type="default"`.
#' They are assigned to the transformed scale of the following parameters:
#' \eqn{\pi_{n}=\frac{\exp(n)}{1+\exp(n)+\exp(a)}}, \eqn{\pi_{a}=\frac{\exp(a)}{1+\exp(n)+\exp(a)}}, 
#' \eqn{{logit}(s_1)=\alpha_s}, \eqn{{logit}(b_1)=\alpha_b}, \eqn{{probit}(u_1)=\alpha_u},
#' and \eqn{{probit}(v_1)=\alpha_v}, where \eqn{n, a \sim N(0, 2.5^2)} and \eqn{\alpha_s, \alpha_b,
#' \alpha_u, \alpha_v \sim N(0, 2^2)}. 
#' Alternatively, users can specify their own prior distributions for all parameters, 
#' and save them as a file \code{prior.study.R} under the same directory with the model 
#' function. By assigning `prior.type = "custom"`, the function calls the user-defined 
#' text string as the priors. See example in `Details`. 
#' Note that if users choose the customized priors, the pre-defined \code{prior.study.R} 
#' must include distributions for all parameters. The function cannot combine the default 
#' priors with partial user-defined prior distributions. 
#' @param digits a positive integer specifying the digits after the decimal point for 
#' the effect size estimates. The default is \code{3}.
#' @param n.adapt the number of iterations for adaptation in Markov chain Monte Carlo (MCMC) algorithm; 
#' it is used to maximize the sampling efficiency. 
#' The default is 1,000. If a warning "adaptation incomplete" appears, users may increase 
#' \code{n.adapt}. This argument and the following \code{n.iter}, \code{n.burnin}, \code{n.chains},
#' \code{n.thin} are passed to the functions in R package \code{rjags}. 
#' @param n.iter the number of iterations of each MCMC chain. 
#' The default is \code{100,000}. 
#' @param n.burnin the number of iterations for burn-in period. The default is 
#' the largest integer not greater than \code{n.iter/2}.
#' @param n.chains the number of MCMC chains. The default is \code{3}. 
#' @param n.thin a positive integer indicating thinning rate for MCMC chains, which is used to 
#' avoid potential high auto-correlation and to save computer memory when \code{n.iter} is 
#' large. The default is set as \code{1} or the largest integer not greater than 
#' \code{((n.iter - n.burnin)/1e+05)}, whichever is larger. 
#' @param conv.diag a logical value indicating whether to compute the Gelman and Rubin 
#' convergence statistic (\eqn{\hat{R}}) of each parameter as a convergence diagnostic.
#' It is considered the chains are well mixed and have converged to the target distribution 
#' if \eqn{\hat{R} \le 1.1}. The default is `FALSE`. If `TRUE`, \code{n.chains} must be greater than 1, 
#' and the function saves each chain's MCMC samples for all parameters, which can be used 
#' to produce trace, posterior density, and auto-correlation plots by calling the function 
#' \code{plot.cacebayes}. 
#' @param mcmc.samples a logical value indicating whether to save MCMC posterior samples
#' in the output object. The default is `FALSE`. If `TRUE`, the output object list 
#' includes each chain's MCMC samples for all parameters. They can be used in the function 
#' \code{plot.cacebayes} to generate the trace, posterior density, and auto-correlation plots 
#' for further model diagnostics. 
#' @param two.step a logical value indicating whether to conduct a two-step meta-analysis. 
#' If `two.step = TRUE`, the posterior mean and standard deviation of study-specific 
#' \eqn{\theta^{CACE}_i} are used to perform a standard meta-analysis, using the R package \code{metafor}. 
#' @param method the method used in meta-analysis if `two.step = TRUE`. The default estimation 
#' method is the REML (restricted maximum-likelihood estimator) method for the random-effects 
#' model. Users can change the argument \code{method} to obtain different meta-analysis 
#' estimators from either a random-effects model or a fixed-effect model, e.g., 
#' `method = "DL"` refers to the DerSimonian--Laird estimator, 
#' `method = "HE"` returns the Hedges estimator, and `method = "HS"` gives the Hunter--Schmidt 
#' estimator.  More details are available from the documentation of the function \code{metafor::rma}. 
#' If the input data include only one study, the meta-analysis result is just the same as 
#' the result from the single study. 
#' 
#' @details  
#' The likelihood \deqn{\log L({\boldsymbol{\beta}}) = N_{000}\log\{\pi_{c}(1-v_1)+\pi_{n}(1-s_1)\}+N_{001}
#' \log(\pi_{c}v_1+\pi_{n}s_1)+N_{010}\log\{{\pi}_{a}(1-b_1)\} + N_{011}\log\{\pi_{a}b_1\}+ N_{100}
#' \log\{\pi_{n}(1-s_1)\}+N_{101}\log({\pi}_{n}s_1) + N_{110}\log\{(\pi_{c}(1-u_1)+\pi_{a}(1-b_1)\}+
#' {N_{111}\log(\pi_{c}u_1+\pi_{a}b_1)} + constant}. If the input \code{data} includes more than one study, the study-specific CACEs will be 
#' estimated by retrieving data row by row.
#' One exmaple of the `prior.study.R` file if using `prior.type = "custom"`:
#' \preformatted{prior.study <- function(prior.type = "custom") {
#' string2 <- "n ~ dnorm(0, 0.01) 
#'   a ~ dnorm(0, 0.01) 
#'   alpha.s ~ dnorm(0, 0.01) 
#'   alpha.b ~ dnorm(0, 0.01) 
#'   alpha.u ~ dnorm(0, 0.01) 
#'   alpha.v ~ dnorm(0, 0.01)"
#' return(string2)}
#' }
#' By default, the function \code{cace.study()} returns a list  
#' including posterior estimates (posterior mean, standard deviation, median, and a 95\% 
#' credible interval (CI) with 2.5\% and 97.5\% quantiles as the lower and upper bounds), 
#' and the deviance information criterion (DIC) statistic for each study.
#' 
#' @examples
#' \dontrun{
#' data("epidural_c", package = "BayesCACE")
#' set.seed(123)
#' out.study <- cace.study(data = epidural_c, conv.diag = TRUE, 
#' mcmc.samples = TRUE, two.step = TRUE) 
#' # Show the estimates of \theta^{CACE} for each single study (posterior mean and 
#' # standard deviation, posterior median, 95% credible interval, and time-series 
#' # standard error):
#' out.study$CACE
#' # If the argument conv.diag is specified as `TRUE`, the output list contains 
#' # a sub-list conv.out, which outputs the Gelman and Rubin convergence statistic,
#' labelled Point est.) calculated for each parameter from each single study, and 
#' their upper confidence limits (labelled Upper C.I.).
#' out.study$conv.out[[1]]
#' }
#' @seealso \code{\link[BayesCACE]{cace.meta.c}}, \code{\link[BayesCACE]{cace.meta.ic}}
#' @references {
#' \insertRef{zhou2020software}{BayesCACE}
#' } 
#' 
cace.study <-
  function(data, param = c("CACE", "u1", "v1", "s1", "b1", "pi.c", "pi.n", 
          "pi.a"), prior.type = "default", digits = 3, n.adapt = 1000, 
           n.iter = 100000, n.burnin = floor(n.iter/2), n.chains = 3, n.thin =  
          max(1,floor((n.iter-n.burnin)/1e+05)), conv.diag = FALSE, mcmc.samples
           = FALSE, two.step = FALSE, method = "REML")    {
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
    
    ## jags model
    modelstring<-model.study(prior.type)
    
    ## data prep
    if(prior.type == "default"){
      Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
      N0 <- n000+n001+n010+n011
      N1 <- n100+n101+n110+n111
      R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
      I <- length(Ntol)
    }
    
    ## parameters to be paramed in jags
    if(!is.element("CACE",param)) param<-c("CACE",param)
    fullparam <- c("CACE", "u1", "v1", "s1", "b1", "pi.c", "pi.n", "pi.a")
    if(!any(is.element(param, fullparam))) stop("parameters must be specified from the following:  
                                                CACE, u1, v1, s1, b1, 
                                                pi.c, pi.n, pi.a \n")
    out<-NULL
    out$model<-"cace.single"    
    
  for (i in 1:I) {
    data.jags <- list(N0=N0[i], N1=N1[i], R=R[i,])
    
    ## jags initial value
    rng.seeds<-sample(1000000,n.chains)
    init.jags <- vector("list", n.chains)
    for(ii in 1:n.chains){
      init.jags[[ii]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
    }
    
    ## run jags
    jags.m<-jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,
                       n.chains=n.chains,n.adapt=n.adapt)
    update(jags.m,n.iter=n.burnin)
    jags.out<-coda.samples.dic(model=jags.m,variable.names=param,n.iter=n.iter,thin=n.thin)
    
    smry<-summary(jags.out$samples)
    smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")], 
                smry$statistics[,c("Naive SE","Time-series SE")])
    smry<-signif(smry,digits=digits)
    out$smry[[i]] <- smry
    
    for (j in 1:length(fullparam)){
      if(is.element(fullparam[j],param)) {
        out[[fullparam[j] ]]<-rbind(out[[fullparam[j] ]], smry[c(fullparam[j]), ])
        }
    }
    
    #dic
    dev<-jags.out$dic[[1]] # mean deviance
    pen<-jags.out$dic[[2]] # pD
    pen.dev<-dev+pen # DIC
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("D.bar","pD","DIC")
    colnames(dic.stat)<-""
    out$DIC[[i]]<-dic.stat
    
    
    if(conv.diag){
      cat("MCMC convergence diagnostic statistics are calculated and saved in conv.out\n")
      conv.out<-gelman.diag(jags.out$samples,multivariate=FALSE)
      out$conv.out[[i]]<-conv.out$psrf
    }
    
    if(mcmc.samples){
      out$mcmc.samples[[i]]<-jags.out$samples
    }
  }
    
  if (two.step) {
    out$meta <- rma(yi=Mean, sei=SD, data=out$CACE, method = method)
  }
  class(out)<-"cace.Bayes"
  return(out)
  options(warn=0)   
}


