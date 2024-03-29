% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moments.r
\name{cumulant2moment}
\alias{cumulant2moment}
\title{Convert raw cumulants to moments.}
\usage{
cumulant2moment(kappa)
}
\arguments{
\item{kappa}{a vector of the raw cumulants. 
The first element is the first cumulant, which is also the first moment.}
}
\value{
a vector of the raw moments.
The first element of the input shall be the same as the first element of the output.
}
\description{
Conversion of a vector of raw cumulatnts to moments.
}
\details{
The \sQuote{raw} cumulants \eqn{\kappa_i}{kappa_i} are connected
to the \sQuote{raw} (uncentered) moments, \eqn{\mu_i'}{mu'_i} via
the equation
\deqn{\mu_n' = \kappa_n + \sum_{m=1}^{n-1} {n-1 \choose m-1} \kappa_m \mu_{n-m}'}
}
\examples{
# normal distribution, mean 0, variance 1
n.mom <- cumulant2moment(c(0,1,0,0,0,0))
# normal distribution, mean 1, variance 1
n.mom <- cumulant2moment(c(1,1,0,0,0,0))
}
\seealso{
\code{\link{moment2cumulant}}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{distribution}
