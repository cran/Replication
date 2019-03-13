\name{llratio.f}
\alias{llratio.f}
\title{Likelihood ratio statistic for constrained versus unconstrained model}
\description{
The function first uses solve.QP to find the best fitting estimates under the imposed
constraints. Subsequently, it calculates the likelihood ratio for the constrained versus unconstrained model.}
\usage{
llratio.f(BKcov, Q, R, r = NULL, E = 0L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BKcov}{BKcov is the variance-covariance matrix of the model parameters.}
  \item{Q}{Q is a vector with the (unconstrained) model estimates.}
  \item{R}{A matrix defining the constraints under which we want to minimize the quadratic function.}
  \item{r}{vector holding the minimimum sum for each row in R.}
  \item{E}{Numeric. The first E constraints are treated as equality constraints, all further as inequality constraints. Default value = 0.}
}
\value{
\item{llratio}{scalar, the value of the resulting likelihood ratio.}
}
\author{M. A. J. Zondervan-Zwijnenburg}
\seealso{\code{\link{solve.QP}}}
\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
aggregate(data$y,by=list(data$g),mean)
reg <- lm(data$y~data$g)
BKcov <- vcov(reg)
Q <- reg$coefficients

#b1=102, b3-b2>20, b3-b4>10, b4>30
R <- rbind(c(1,0,0,0),c(0,-1,1,0),c(0,0,1,-1),c(0,0,0,1))
r = c(102,20,10,30)
E = 1

llratio.f(BKcov=BKcov, Q=Q, R=R, r=r, E = E)
}
\keyword{htest}
