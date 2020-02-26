#' @return
#' @export
#' 
plot.noncomp <- 
  function(data, overall=TRUE) {
  if(missing(data)) stop("need to specify a data set with complete compliance information. \n")
  
  if (any((data$n000==0 & data$n001==0 & data$n010==0 & data$n011==0) |
      (data$n100==0 & data$n101==0 & data$n110==0 & data$n111==0)) ) {
    warning ("noncompliance plot can only be made for studies with complete data. \n
             only trials with full compliance information are kept. \n")
    data <- data[ !((data$n000==0 & data$n001==0 & data$n010==0 & data$n011==0) |
                    (data$n100==0 & data$n101==0 & data$n110==0 & data$n111==0)),]
  }
    

  data$n00s <- data$n000+data$n001
  data$n01s <- data$n010+data$n011
  data$n10s <- data$n100+data$n101
  data$n11s <- data$n110+data$n111
  data$n0ss <- data$n00s+data$n01s
  data$n1ss <- data$n10s+data$n11s                 
  data$p10 <- data$n01s/data$n0ss   # P(T=1|R=0)
  data$p01 <- data$n10s/data$n1ss   # P(T=0|R=1)

  tmp.b <- t(sapply(split(data, data$study.id), function(x) binom.test(x$n01s, x$n0ss)$conf.int))
  tmp.a <- t(sapply(split(data, data$study.id), function(x) binom.test(x$n10s, x$n1ss)$conf.int))
  data$lower.p10 <- tmp.b[,1]
  data$upper.p10 <- tmp.b[,2]
  data$lower.p01 <- tmp.a[,1]
  data$upper.p01 <- tmp.a[,2]
  
  expit <- function(x){y=exp(x)/(1+exp(x)) 
    return(y)}  
  C <- 16*sqrt(3)/(15*pi)
  
  m10 <- glmer(cbind(n01s, n00s) ~  (1 | study.id),
              data = data, family = binomial)
  mu10 <- c(coef(summary(m10))[1], confint(m10, method="Wald")["(Intercept)", c("2.5 %", "97.5 %")])
  sd10 <- as.data.frame(VarCorr(m10))[1, "sdcor"]
  outp10 <- expit(mu10/sqrt(1+C^2*sd10^2))
  
  p10_out <- paste(formatC(round(outp10[1], 3), format='f', digits=3 ), " (", 
                   formatC(round(outp10[2], 3), format='f', digits=3 ), ",",
                   formatC(round(outp10[3], 3), format='f', digits=3 ), ")", sep="")
  

  m01 <- glmer(cbind(n10s, n11s) ~  (1 | study.id),
              data = data, family = binomial)
  mu01 <- c(coef(summary(m01))[1], confint(m01, method="Wald")["(Intercept)", c("2.5 %", "97.5 %")])
  sd01 <- as.data.frame(VarCorr(m01))[1, "sdcor"]
  outp01 <- expit(mu01/sqrt(1+C^2*sd01^2))
  
  p01_out <- paste(formatC(round(outp01[1], 3), format='f', digits=3 ), " (", 
                   formatC(round(outp01[2], 3), format='f', digits=3 ), ",",
                   formatC(round(outp01[3], 3), format='f', digits=3 ), ")", sep="")

  # outp10 <- rma(xi=n01s, ni=n0ss, measure="PLO", method="ML", data=data)  #p10
  # p10_out <- paste(formatC(round(expit(outp10$beta), 3), format='f', digits=3 ), " (",
  #                  formatC(round(expit(outp10$ci.lb), 3), format='f', digits=3 ),",", 
  #                  formatC(round(expit(outp10$ci.ub), 3), format='f', digits=3 ), ")", sep="")
  # 
  # outp01 <- rma(xi=n10s, ni=n1ss, measure="PLO", method="ML", data=data)  #p01
  # p01_out <- paste(formatC(round(expit(outp01$beta), 3), format='f', digits=3 ), " (",
  #                  formatC(round(expit(outp01$ci.lb), 3), format='f', digits=3 ),",", 
  #                  formatC(round(expit(outp01$ci.ub), 3), format='f', digits=3 ), ")", sep="")
  
if (overall) {
  tabletext<-cbind(
    c("Study (Author, Year)", paste(data$study.name), "Overall"),
    c("P(T=0|R=1)", 
      paste(formatC(round(data$p01, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p01, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p01, 3), format='f', digits=3 ), ")", sep=""), 
      p01_out  
    ),
    c("P(T=1|R=0)", 
      paste(formatC(round(data$p10, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p10, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p10, 3), digits=3 ), ")", sep=""), 
      p10_out  
    )
  )
  
  # tabletext<-cbind(c("Study (Author, Year)", paste(data$study.name), "Overall"))
  
  xticks <- seq(from = 0, to = max(data$upper.p01, data$upper.p10), by = 0.1)
  xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
  attr(xticks, "labels") <- xtlab
  own.f <- fpTxtGp(ticks = gpar(cex=0.85), xlab  = gpar(cex = 0.95))

  forestplot(tabletext,
     graph.pos = 3,
     hrzl_lines = gpar(lwd=1, col="#444444"),
     legend_args = fpLegend(pos = list("top", "inset"=.03, "align"="horizontal")),
     legend = c("P(T=0|R=1)", "P(T=1|R=0)"),
     fn.ci_norm = c(fpDrawCircleCI, fpDrawNormalCI),
     boxsize = .15, # We set the box size to better visualize the type
     line.margin = .36, # We need to add this to avoid crowding
     is.summary=c(TRUE, rep(FALSE, nrow(data)), TRUE),
     mean = cbind(c(NA, data$p01, round(outp01[1], 3)), 
                  c(NA, data$p10, round(outp10[1], 3))),
     lower = cbind(c(NA, data$lower.p01, round(outp01[2], 3)), 
                   c(NA, data$lower.p10, round(outp10[2], 3))),
     upper = cbind(c(NA, data$upper.p01, round(outp01[3], 3)), 
                   c(NA, data$upper.p10, round(outp10[3], 3))),
     clip =c(0, 1),
     col=fpColors(box=c("darkred", "blue"), summary=c("darkred", "blue")),
     grid = structure(c(0, 0.5), gp = gpar(lty = 2, col = "#CCCCFF")), 
     xticks = xticks,
     txt_gp = own.f,
     xlab="Nocompliance Rates")
}
  
else {
  tabletext<-cbind(
    c("Study (Author, Year)", paste(data$study.name)),
    c("P(T=0|R=1)", 
      paste(formatC(round(data$p01, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p01, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p01, 3), format='f', digits=3 ), ")", sep="")
    ),
    c("P(T=1|R=0)", 
      paste(formatC(round(data$p10, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p10, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p10, 3), digits=3 ), ")", sep="")
    )
  )
  
  xticks <- seq(from = 0, to = max(data$upper.p01, data$upper.p10), by = 0.1)
  xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
  attr(xticks, "labels") <- xtlab
  own.f <- fpTxtGp(ticks = gpar(cex=0.85), xlab  = gpar(cex = 0.95))
  
  forestplot(tabletext,
             graph.pos = 3,
             hrzl_lines = gpar(lwd=1, col="#444444"),
             legend_args = fpLegend(pos = list("top", "inset"=.03, "align"="horizontal")),
             legend = c("P(T=0|R=1)", "P(T=1|R=0)"),
             fn.ci_norm = c(fpDrawCircleCI, fpDrawNormalCI),
             boxsize = .15, # We set the box size to better visualize the type
             line.margin = .36, # We need to add this to avoid crowding
             is.summary=c(TRUE, rep(FALSE, nrow(data))),
             mean = cbind(c(NA, data$p01), 
                          c(NA, data$p10)),
             lower = cbind(c(NA, data$lower.p01), 
                           c(NA, data$lower.p10)),
             upper = cbind(c(NA, data$upper.p01), 
                           c(NA, data$upper.p10)),
             clip =c(0, 1),
             col=fpColors(box=c("darkred", "blue"), summary=c("darkred", "blue")),
             grid = structure(c(0, 0.5), gp = gpar(lty = 2, col = "#CCCCFF")), 
             xticks = xticks,
             txt_gp = own.f,
             xlab="Nocompliance Rates")
}
  
  invisible()
  }

