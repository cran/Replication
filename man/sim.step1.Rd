\name{sim.step1}
\alias{sim.step1}
\title{Simulate Data for the Predictive Distribution}
\description{
Simulates data y.s using lavaan simulateData. The data y.s are based on samples from the posterior and represent samples from the predictive distribution. }
\usage{
sim.step1(n.r, nsim, post, pT, free.i, group=NULL, n.groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.r}{Sample size for the new data (y.r)}
  \item{nsim}{Number of datasets y.s to be created.}
  \item{post}{Matrix with samples from the posterior distribution for the original data y.o}
  \item{pT}{The blavaan parameter table. }
  \item{free.i}{A vector with indices for the freely estimated parameters as given in pT.}
  \item{group}{If applicable, the name of the grouping variable.}
  \item{n.groups}{If applicable, the number of groups.}
}
\value{
\item{y.s}{A list containing all y.s.}
}
\author{
M. A. J. Zondervan-Zwijnenburg}

\seealso{
\code{\link{simulateData}}
}
\keyword{models }
