\name{ppc.step2step3}
\alias{ppc.step2step3}
\title{Prior predictive check step 2 and 3}
\description{
Calculates an approximate likelihood ratio (D) for new data (y.r) and predicted data (y.s) according to the proposed constraints, and generates a prior predictive p-value.}
\usage{
ppc.step2step3(step1, y.r, model = model, H0, s.i,
        ordered = NULL, sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
        group = NULL, cluster = NULL, constraints = "", WLS.V = NULL, NACOV = NULL,
        bayes = FALSE, dp = NULL, convergence = "manual", nchains = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{step1}{An object containing the output of ppc.step1.}
  \item{y.r}{A data.frame with the new data. If y.r = NA, the approximate likelihood ratio will only be computed for the predicted data y.s.}
  \item{model}{The (b)lavaan model that is to be fitted to the data.}
  \item{H0}{The replication hypothesis within quotes "" with the lavaan plabels as parameter names and parts separated with &. For more information on hypothesis specification, see the details section below.}
  \item{s.i}{A vector of length p with indices for the (pooled) standard deviation parameters with which the effect sizes should be computed. Default = NULL.}
  \item{ordered}{Character vector. Only used if the data is in a data.frame. Treat these variables as ordered (ordinal) variables, if they are endogenous in the model. Importantly, all other variables will be treated as numeric (unless they are declared as ordered in the original data.frame}
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
  \item{dp}{If bayes = TRUE, blavaan default prior distributions on different types of parameters, typically the result of a call to dpriors(). See the dpriors() help file for more information.}
    \item{convergence}{If bayes = TRUE, default convergence setting = "manual". If "auto", parameters will be sampled until convergence is achieved (via autorun.jags). In this case, the arguments burnin and sample are passed to autorun.jags as startburnin and startsample, respectively. Otherwise, parameters are sampled as specified by the user (or by the run.jags defaults).}
  \item{nchains}{If bayes = TRUE, A scalar indicating the number of chains to be used in the Bayesian analysis. Default value = 2.}
  }

\details{
The specification of `H0` in `ppc.step2step3`:

`H0` is a character string that specifies which informative
hypothesis has to be evaluated. A simple example is \code{H0 <- ".p1. >
.p2. > .p3. & .p1. = 2"} which specifies a hypothesis using three estimates with
names ".p1.", ".p2.", and ".p3.", respectively.

The hypothesis specified has to adhere to the following rules:
\itemize{
\item When using \code{ppc.step2step3}, the `plabels` of the \code{blavaan} output resulting from \code{ppc.step1} have to be used to indicate which parameters are involved in the informative hypothesis H0. For example `.p1.` and `.p2.` can be the labels of the parameters of interest.

\item Linear combinations of parameters must be specified adhering to the
following rules:
a) Each parameter name is used at most once.
b) Each parameter name may or may not be pre-multiplied with a number.
c) A constant may be added or subtracted from each parameter name.
Examples are: \code{"3 *.p1.+ 5"}; \code{".p1.+ 2 * .p2.+ 3 * .p3.- 2"} and \code{".p1.- .p2."}.

\item (Linear combinations of) parameters can be constrained using \code{<}, \code{>}, and
\code{=}. For example, \code{".p1.> 0"} or
\code{".p1.> .p2.= 0"} or \code{"2 *.p1.< .p2.+ .p3.> 5"}.

\item The ampersand \code{&} can be used to combine different parts of a hypothesis.
For example, \code{".p1.> .p2.& .p2.> .p3."} which is equivalent to \code{".p1.> .p2.> .p3."} or
\code{".p1.> 0 & .p2.> 0 & .p3.> 0"}.

\item Sets of (linear combinations of) parameters subjected to the same
constraints can be specified using (). For
example, \code{".p1.> (.p2.,.p3.)"} which is equivalent to \code{".p1.> .p2.&.p1.> .p3."}.

\item Hypotheses have to be possible. A hypothesis is impossible if estimates in agreement with the hypothesis do not
exist. For example: values for \code{.p1.} in agreement with \code{".p1.= 0 &
.p1. > 2"} do not exist. It is the responsibility of the user to ensure that the
hypotheses specified are possible. If not, \code{ppc.step2step}` will return an
error message: \code{Error in solve.QP(Dmat, dvec = dvec, t(R), r, meq = E, factorized = FALSE) : constraints are inconsistent, no solution!}.}
}

\value{
Generates a histogram of llratio.s in which llratio.r is indicated with a vertical line. The proportion of llratio.s at the right of this line constitutes the prior predictive p-value.
\item{llratio.r}{The likelihood ratio for the new dataset.}
\item{p-value}{The prior predictive p-value.}
\item{llratio.s}{The likelihood ratio's for each of the datasets y.s.}
\item{H0 matrices}{Matrices R, r, and E that constitute the matrix form of H0: R*theta>r or R*theta=r for equality constraints. The value in E specifies the number of equality constraints in the hypothesis.}
\item{pT.s}{Parameter table columns for y.s indicating which parameter received which id and which label in the lavaan analysis of the predicted data. If the id's and labels differ from those of the blavaan analysis used in step1, the function will generate a warning for the user to check whether this difference affects parameters in H0.}
}
\author{M. A. J. Zondervan-Zwijnenburg}
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

step23.reg <- ppc.step2step3(step1=step1.reg,y.r=y.r,model=model,H0=".p2. = -0.878&.p1. > 0.301")
}
\donttest{
#the following example can be used, but may take >10 seconds

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
step1.reg <- ppc.step1(y.o=y.o,model=model,nchains=2,n.r=50)

print(step1.reg$pT)
#H0: #reg > est, int = est       =       B1>0.302 & B0= -0.878
pT <- step1.reg$pT #parameter table
int.id <- which(pT$lhs=="y"&pT$op=="~1"&pT$rhs=="") #identify B0
reg.id <- which(pT$lhs=="y"&pT$op=="~"&pT$rhs=="x") #identify B1
hyp <- cbind(pT[c(int.id,reg.id),"plabel"],c("=",">"),c(pT[c(int.id,reg.id),"est"]))
print(hyp)
H0 <- paste(hyp[,1],hyp[,2],hyp[,3],collapse="&")

step23.reg <- ppc.step2step3(step1=step1.reg,y.r=y.r,model=model,H0)
}
}
\keyword{htest }% use one of  RShowDoc("KEYWORDS")
\keyword{models }% __ONLY ONE__ keyword per line
