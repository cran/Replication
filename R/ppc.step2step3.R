ppc.step2step3 <- function(step1, y.r, model=model, H0,
                           s.i=NULL,
                           ordered = NULL, sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
                           group = NULL, cluster = NULL, constraints = "", WLS.V = NULL, NACOV = NULL,
                           bayes=FALSE,dp=NULL,convergence="manual",nchains=2){

  #compute likelihood ratio (i.e., D) for all y.s and for y.r

  #f(D_y[s])
  #limited number of posterior samples (saves comp. time), replace=TRUE (!)
  #sample.r <- sample(length(y.s),n.sample.r,replace=TRUE)

  y.s <- step1$y.s
  pT1 <- step1$pT

  pT1 <- pT1[which(pT1$free!=0),]
  pT1 <- pT1[!(duplicated(pT1$label))|pT1$label=="",]

  vars <- pT1$plabel             #var names for estimated parameters
  mat <- create_matrices(varnames=c(vars),hyp=H0)
  R <- mat$R
  r <- mat$r
  E <- mat$E

  llratio.s <- list()
  print("Calculating likelihood ratio for each y.s",quote=FALSE)
  pb = txtProgressBar(min = 0, max = length(y.s), initial = 0, style=3, width = 50)

  if(bayes==TRUE){
    for (i in 1:length(y.s)){
      setTxtProgressBar(pb,i)
      fit_l <- bsem(model, data=y.s[[i]],
                    n.chains=nchains, dp=dp,
                    sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                    group = group, constraints = "", WLS.V = WLS.V, NACOV = NACOV, convergence=convergence)
      pT <- parameterTable(fit_l)
      free.i <- which(pT$free!=0)
      BKcov <- lavInspect(fit_l,"vcov")

      if (sum(duplicated(rownames(BKcov)))>0){
        pT <- pT[free.i,][-which(duplicated(rownames(BKcov))),]
        Q <- pT$est
        BKcov <- BKcov[!duplicated(rownames(BKcov)), !duplicated(colnames(BKcov))]
      }else{
        Q <- pT$est[free.i]
      }

      if(is.null(s.i)==FALSE){
        s <- vector()
        for (j in 1:length(r)){s[j] <- pT$est[ pT$id == s.i[j] ]}
        r.e <- r*s

        llratio.s[[i]] <-tryCatch(llratio.f(BKcov=BKcov,Q=Q,R=R,r=r.e,E=E),
                                  error=function(e) NA)

        llratio.s <- na.omit(unlist(llratio.s))
      }else{

        llratio.s[[i]] <-tryCatch(llratio.f(BKcov=BKcov,Q=Q,R=R,r=r,E=E),
                                  error=function(e) NA)

        llratio.s <- na.omit(unlist(llratio.s))}
    }
  }else{
    for (i in 1:length(y.s)){
      setTxtProgressBar(pb,i)
      fit_l <- sem(model, data=y.s[[i]],
                   ordered = ordered, sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                   group = group, cluster = cluster, constraints = "", WLS.V = WLS.V, NACOV = NACOV)
      pT.s <- parameterTable(fit_l)
      free.i <- which(pT.s$free!=0)
      BKcov <- lavInspect(fit_l,"vcov")

      if (sum(duplicated(rownames(BKcov)))>0){
        pT <- pT.s[free.i,][-which(duplicated(rownames(BKcov))),]
        Q <- pT$est
        BKcov <- BKcov[!duplicated(rownames(BKcov)), !duplicated(colnames(BKcov))]
      }else{
        Q <- pT.s$est[free.i]
      }

      if(is.null(s.i)==FALSE){
        s <- vector()
        for (j in 1:length(r)){s[j] <- pT.s$est[ pT.s$id == s.i[j] ]}
        r.e <- r*s

        llratio.s[[i]] <-tryCatch(llratio.f(BKcov=BKcov,Q=Q,R=R,r=r.e,E=E),
                                  error=function(e) NA)
      }else{

        llratio.s[[i]] <-tryCatch(llratio.f(BKcov=BKcov,Q=Q,R=R,r=r,E=E),
                                  error=function(e) NA)
      }

    }

    llratio.s <- na.omit(unlist(llratio.s))

    if(length(attr(llratio.s,"na.action"))!=0){
      print(paste(length(attr(llratio.s,"na.action")),"datasets could not be analyzed properly, this may relate to non-positive definite variance-covariance matrices."))}

    pT.s <- pT.s[free.i,]
    pT.s <- pT.s[!(duplicated(pT.s$label))|pT.s$label=="",]
    if(identical(pT.s[,1:4],pT1[,1:4])==FALSE){
      print("Warning: the Bayesian parameter table of step1 is not equal to that of step2step3. Check pT.1 and pT.s in the results to see if this affects parameter labels for parameters in H0. If so, specify all model parameters in the model syntax.")
    }
  }

  if(is.null(y.r)==FALSE){
    #for y.r
    fit_r <- sem(model, data=y.r,
                 ordered = ordered, sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                 group = group, cluster = cluster, constraints = "", WLS.V = WLS.V, NACOV = NACOV)
    BKcov.r <- lavInspect(fit_r,"vcov")
    pT.r <- parameterTable(fit_r)
    free.i <- which(pT.r$free!=0)

    if (sum(duplicated(rownames(BKcov.r)))>0){
      pT.r <- pT.r[free.i,][-which(duplicated(rownames(BKcov.r))),]
      Q.r <- pT.r$est
      BKcov.r <- BKcov.r[!duplicated(rownames(BKcov.r)), !duplicated(colnames(BKcov.r))]
    }else{
      Q.r <- pT.r$est[free.i]
    }

    if(is.null(s.i)==FALSE){
      s <- vector()
      for (j in 1:length(r)){s[j] <- pT$est[ pT$id == s.i[j] ]}
      r.e <- r*s
      llratio.r <- llratio.f(BKcov=BKcov.r,Q=Q.r,R=R,r=r.e,E=E)
    }else{
      llratio.r <- llratio.f(BKcov=BKcov.r,Q=Q.r,R=R,r=r,E=E)}

    #plot results
    llratio.s <<- llratio.s
    ppc.plot(llratio.s,llratio.r)

    #prior predictive p
    p <- sum((llratio.s)>=llratio.r)/length(llratio.s) #prior predictive p-value
    results <- list("llratio.r"=llratio.r,"p-value"=p,"llratio.s"=llratio.s,"H0 matrices"=mat,
                    "pT.s"=pT.s[,c(1:4,12)],"pT.1"=pT1[,c(1:4,12)])
  }else{
    results <- list("llratio.s"=llratio.s,"H0 matrices"=mat,
                    "pT.s"=pT.s[,c(1:4,12)],"pT.1"=pT1[,c(1:4,12)])
  }
  cat("\n")
  return(results)
}
