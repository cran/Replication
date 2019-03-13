#obtain posterior ####
posterior.step1 <- function(y.o, model,
                            sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, group = NULL,
                            constraints = "", WLS.V = NULL, NACOV = NULL,
                            nchains=2, nadapt, nburnin, nsample, dp=NULL,convergence="manual",
                            imp=imp){

  #if missing data, impute and combine posterior distributions
  if(is.null(imp) == FALSE){
    y.o.imp <- list()
    post.list <- list()
    pT <- list()
    for (i in 1:imp$m){
      y.o.imp[[i]] <- complete(imp,i)
      y.o.imp[[i]] <- data.frame(sapply(y.o.imp[[i]] , function(x) as.numeric(as.character(x))))
      fit <- bsem(model, data=y.o.imp[[i]], dp=dp,
                  sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                  group = group, constraints = "", WLS.V = WLS.V, NACOV = NACOV,
                  n.chains=nchains, adapt = nadapt, burnin=nburnin, sample=nsample, convergence=convergence)

      pT[[i]] <- parameterTable(fit)
      free.i <- which(pT[[1]]$free!=0)
      fix.i <- which(pT[[1]]$free==0)
      est.fix<-pT[[1]]$est[fix.i]
      traceplot <- plot(fit, pars=1:length(free.i), plot.type="trace")

      post.y.o <- blavInspect(fit,"mcmc")
      post <- list()
      for(j in 1:length(free.i)){
        post[[j]] <- unlist(post.y.o[,j])} #posterior per parameter, combine chains

      post.list[[i]] <- post #posterior parameter samples per imputation
    }

    #save average score for imputation varying values in pT
    pT.m.start <- rowMeans(sapply(pT,'[[',"start"))
    pT.m.est <- rowMeans(sapply(pT,'[[',"est")); pT.m.se <- rowMeans(sapply(pT,'[[',"se"))
    pT.m.psrf <- rowMeans(sapply(pT,'[[',"psrf")); pT.m.logBF <- rowMeans(sapply(pT,'[[',"logBF"))
    pT <- pT[[1]]
    pT$start <- pT.m.start
    pT$est <- pT.m.est; pT$se <- pT.m.se
    pT$psrf <- pT.m.psrf; pT$logBF <- pT.m.logBF

    #save posterior
    post <- list()
    for(j in 1:length(free.i)){
      post[[j]] <- as.vector(sapply(post.list,'[[',j))} #posterior per parameter, combine imputations

  }

  if(is.null(imp)==TRUE){
    fit <- bsem(model, data=y.o, dp=dp,
                sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs,
                group = group, constraints = "", WLS.V = WLS.V, NACOV = NACOV,
                n.chains=nchains, adapt=nadapt, burnin=nburnin, sample=nsample, convergence=convergence)
    #store fixed and free output
    pT <- parameterTable(fit)
    free.i <- which(pT$free!=0)
    traceplot <- plot(fit, pars=1:length(free.i), plot.type="trace")

    #obtain samples posterior
    post.y.o <- blavInspect(fit,"mcmc")
    post <<- list()
    for(i in 1:length(free.i)){post[[i]] <<- unlist(post.y.o[,i])}
  }

  #results
  results <- list(post=post, pT=pT,free.i=free.i, traceplot)
  return(results)

}
