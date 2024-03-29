% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gram_charlier.r
\name{dapx_gca}
\alias{dapx_gca}
\alias{papx_gca}
\title{Approximate density and distribution via Gram-Charlier A expansion.}
\usage{
dapx_gca(x, raw.moments, support=NULL, 
 basis=c('normal','gamma','beta','arcsine','wigner'), 
 basepar=NULL, log=FALSE)

papx_gca(q, raw.moments, support=NULL, 
 basis=c('normal','gamma','beta','arcsine','wigner'), 
 basepar=NULL, lower.tail=TRUE, log.p=FALSE)
}
\arguments{
\item{x}{where to evaluate the approximate density.}

\item{raw.moments}{an atomic array of the 1st through kth raw moments
of the probability distribution.}

\item{support}{the support of the density function. It is assumed
that the density is zero on the complement of this open interval.
This defaults to \code{c(-Inf,Inf)} for the normal basis,
\code{c(0,Inf)} for the gamma basis, and
\code{c(0,1)} for the Beta, and 
\code{c(-1,1)} for the arcsine and wigner.}

\item{basis}{the basis under which to perform the approximation. \code{'normal'}
gives the classical \dQuote{A} series expansion around the PDF and CDF of the normal
distribution via Hermite polynomials. \code{'gamma'} expands around a
gamma distribution with parameters \code{basepar$shape} and
\code{basepar$scale}. 
\code{'beta'} expands around a beta distribution with parameters
\code{basepar$shape1} and \code{basepar$shape2}.}

\item{basepar}{the parameters for the base distribution approximation. 
If \code{NULL}, the shape and rate are inferred from the first two moments
and/or from the \code{support} as appropriate.}

\item{log}{logical; if TRUE, densities \eqn{f} are given 
as \eqn{\mbox{log}(f)}{log(f)}.}

\item{q}{where to evaluate the approximate distribution.}

\item{log.p}{logical; if TRUE, probabilities p are given 
as \eqn{\mbox{log}(p)}{log(p)}.}

\item{lower.tail}{whether to compute the lower tail. If false, we approximate the survival function.}
}
\value{
The approximate density at \code{x}, or the approximate CDF at
\code{q}.
}
\description{
Approximate the probability density or cumulative distribution function of a distribution via its raw moments.
}
\details{
Given the raw moments of a probability distribution, we can approximate the probability 
density function, or the cumulative distribution function, via a Gram-Charlier 
expansion on the standardized distribution. This expansion uses
some weighting function, \eqn{w}{w}, typically the density of some \sQuote{parent}
probability distribution, and polynomials, \eqn{p_n}{p_n} which are
orthogonal with respect to that weighting:
\deqn{\int_{-\infty}^{\infty} p_n(x) p_m(x) w(x) \mathrm{d}x = h_n \delta_{mn}.}{integral_-inf^inf p_n(x) p_m(x) w(x) dx = h_n delta_mn.}

Let \eqn{f(x)}{f(x)} be the probability density of some random variable,
with cumulative distribution function \eqn{F(x)}{F(x)}. We express
\deqn{f(x) = \sum_{n \ge 0} c_n p_n(x) w(x)}{f(x) = sum_n c_n p_n(x) w(x)}
The constants \eqn{c_n}{c_n} can be computed from the known moments
of the distribution.

For the Gram Charlier \dQuote{A} series, the weighting function is the PDF of the
normal distribution, and the polynomials are the (probabilist's) Hermite
polynomials. As a weighting function, one can also use the PDF of the gamma
distribution (resulting in generalized Laguerre polynomials), or the
PDF of the Beta distribution (resulting in Jacobi polynomials).
}
\note{
Monotonicity of the CDF is not guaranteed.
}
\examples{
# normal distribution:
xvals <- seq(-2,2,length.out=501)
d1 <- dapx_gca(xvals, c(0,1,0,3,0), basis='normal')
d2 <- dnorm(xvals)
# they should match:
d1 - d2

qvals <- seq(-2,2,length.out=501)
p1 <- papx_gca(qvals, c(0,1,0,3,0))
p2 <- pnorm(qvals)
p1 - p2

xvals <- seq(-6,6,length.out=501)
mu <- 2
sigma <- 3
raw.moments <- c(2,13,62,475,3182)
d1 <- dapx_gca(xvals, raw.moments, basis='normal')
d2 <- dnorm(xvals,mean=mu,sd=sigma)
\donttest{
plot(xvals,d1)
lines(xvals,d2,col='red')
}
p1 <- papx_gca(xvals, raw.moments, basis='normal')
p2 <- pnorm(xvals,mean=mu,sd=sigma)
\donttest{
plot(xvals,p1)
lines(xvals,p2,col='red')
}

# for a one-sided distribution, like the chi-square
chidf <- 30
ords <- seq(1,9)
raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
xvals <- seq(0.3,10,length.out=501)
d1g <- dapx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
d2 <- dchisq(xvals,df=chidf)
\donttest{
plot(xvals,d1g)
lines(xvals,d2,col='red')
}

p1g <- papx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
p2 <- pchisq(xvals,df=chidf)
\donttest{
plot(xvals,p1g)
lines(xvals,p2,col='red')
}

# for a one-sided distribution, like the log-normal
mu <- 2
sigma <- 1
ords <- seq(1,8)
raw.moments <- exp(ords * mu + 0.5 * (sigma*ords)^2)
xvals <- seq(0.5,10,length.out=501)
d1g <- dapx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
d2 <- dnorm(log(xvals),mean=mu,sd=sigma) / xvals
\donttest{
	plot(xvals,d1g)
	lines(xvals,d2,col='red')
}
}
\references{
Jaschke, Stefan R. "The Cornish-Fisher-expansion in the context of 
Delta-Gamma-normal approximations." No. 2001, 54. Discussion Papers, 
Interdisciplinary Research Project 373: Quantification and Simulation of 
Economic Processes, 2001. 
\url{http://www.jaschke-net.de/papers/CoFi.pdf}

S. Blinnikov and R. Moessner. "Expansions for nearly Gaussian
distributions." Astronomy and Astrophysics Supplement 130 (1998): 193-205.
\url{http://arxiv.org/abs/astro-ph/9711239}
}
\seealso{
\code{\link{qapx_cf}}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{distribution}
