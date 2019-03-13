\name{llratio.imp}
\alias{llratio.imp}
\title{Log likelihood ratio and p-value for data imputed with mice}
\description{
Computes the log likelihood ratio and p-value for data that is imputed with mice}
\usage{
llratio.imp(step2step3,imp,model,effectsize=FALSE,s.i,
            sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
            group = NULL, cluster = NULL, constraints = "", WLS.V = NULL, NACOV = NULL,
            bayes=FALSE,dp=NULL,nchains=2)
}
\arguments{
  \item{step2step3}{The likelihood ratio values for simulated data as obtained with the ppc.step2step3 function.}
  \item{imp}{A mids object created with the R-package mice}
  \item{model}{The lavaan model that is to be applied to the data}
  \item{effectsize}{Logic; if TRUE, the constraints concern effectsizes.}
  \item{s.i}{A vector of length p holding with indices for the (pooled) standard deviation parameters with which the effect sizes should be computed}
  \item{sample.cov}{Numeric matrix. A sample variance-covariance matrix. The rownames and/or colnames must contain the observed variable names. For a multiple group analysis, a list with a variance-covariance matrix for each group. Note that if maximum likelihood estimation is used and likelihood="normal", the user provided covariance matrix is internally rescaled by multiplying it with a factor (N-1)/N, to ensure that the covariance matrix has been divided by N. This can be turned off by setting the sample.cov.rescale argument to FALSE.
}
  \item{sample.mean}{A sample mean vector. For a multiple group analysis, a list with a mean vector for each group.}
  \item{sample.nobs}{
Number of observations if the full data frame is missing and only sample moments are given. For a multiple group analysis, a list or a vector with the number of observations for each group.}
  \item{group}{
A variable name in the data frame defining the groups in a multiple group analysis.}
  \item{cluster}{The cluster variable for multilevel data (beta!).}
  \item{constraints}{
Additional (in)equality constraints not yet included in the model syntax. See model.syntax for more information. Note that the replication hypothesis should not be specified here!}
  \item{WLS.V}{
A user provided weight matrix to be used by estimator "WLS"; if the estimator is "DWLS", only the diagonal of this matrix will be used. For a multiple group analysis, a list with a weight matrix for each group. The elements of the weight matrix should be in the following order (if all data is continuous): first the means (if a meanstructure is involved), then the lower triangular elements of the covariance matrix including the diagonal, ordered column by column. In the categorical case: first the thresholds (including the means for continuous variables), then the slopes (if any), the variances of continuous variables (if any), and finally the lower triangular elements of the correlation/covariance matrix excluding the diagonal, ordered column by column.}
  \item{NACOV}{
A user provided matrix containing the elements of (N times) the asymptotic variance-covariance matrix of the sample statistics. For a multiple group analysis, a list with an asymptotic variance-covariance matrix for each group. See the WLS.V argument for information about the order of the elements.}
  \item{bayes}{Logic; if TRUE, a Bayesian estimator is used. }
  \item{dp}{blavaan default prior distributions on different types of parameters, typically the result of a call to dpriors(). See the dpriors() help file for more information.}
  \item{nchains}{A scalar indicating the number of chains to be used in the Bayesian analysis. Default value = 2.}
}
\value{
\item{pT}{The average parameter table}
\item{llratio.i}{The likelihood ratio values for each of the imputed datasets}
\item{pvals}{The prior-predictive p-values for each of the imputed datasets}
}
\author{M. A. J. Zondervan-Zwijnenburg}

\seealso{
\code{\link{llratio.f}}
}
\examples{
\dontshow{
#create data
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

# simple regression -------------------------------------------------------

set.seed(9)
#step 1 input
#create/load data
n.o=30 #sample size original data
y.o <- data.frame(y=rnorm2(n.o,0,1),x=rnorm2(n.o,3,1))
#y.o <- correlate(as.matrix(y.o), corm=.70); y.o <- data.frame(y=y.o[,1],x=y.o[,2])
n.r=50 #sample size new data
y.r <- data.frame(y=rnorm2(n.r,0.5,1),x=rnorm2(n.r,3,1))

#blavaan model
model <- '
y ~ x     #regression
y ~1      #intercept not default in lavaan (but is in blavaan)
'

#Warning: This is a minimal example;
step1.reg <- ppc.step1(y.o=y.o,model=model,nchains=2,n.r=10,nadapt=10,nburnin=10,nsample=10,nsim=10)

step23.reg <- ppc.step2step3(step1=step1.reg,y.r=NULL,model=model,H0=".p2. = -0.878&.p1. > 0.301")

y.r$y[runif(5,1,n.r)] <- NA #random missing data
imp <- mice(data=y.r,maxit=2,m=3)
llratio.imp(step2step3=step23.reg,imp=imp,model=model)
}\donttest{
#the following example can be used, but takes >10 seconds

#create data
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

#step 1 input
#create/load data
n.o=30 #sample size original data
y.o <- data.frame(y=rnorm2(n.o,0,1),x=rnorm2(n.o,3,1))
n.r=80 #sample size new data
y.r <- data.frame(y=rnorm2(n.r,0.5,1),x=rnorm2(n.r,3,1))
y.r$y[runif(5,1,n.r)] <- NA #random missing data

#blavaan model
model <- '
y ~ x     #regression
y ~1      #intercept not default in lavaan (but is in blavaan)
'

step1.reg <- ppc.step1(y.o=y.o,model=model,n.r=n.r)

#H0: #reg > est, int = est       =       B1>0.302 & B0= -0.878
pT <- step1.reg$pT #parameter table
int.id <- which(pT$lhs=="y"&pT$op=="~1"&pT$rhs=="") #identify B0
reg.id <- which(pT$lhs=="y"&pT$op=="~"&pT$rhs=="x") #identify B1
hyp <- cbind(pT[c(int.id,reg.id),"plabel"],c("=",">"),c(pT[c(int.id,reg.id),"est"]))
print(hyp)
H0 <- paste(hyp[,1],hyp[,2],hyp[,3],collapse="&")

step23.reg <- ppc.step2step3(step1=step1.reg,y.r=NULL,model=model,H0)

y.r$y[runif(5,1,n.r)] <- NA #random missing data
imp <- mice(y.r,maxit=10,m=10)
llratio.imp(step2step3=step23.reg,imp=imp,model=model)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
