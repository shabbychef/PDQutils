#' @details
#'
#' Given the raw moments of a probability distribution, we can approximate the probability 
#' density function, or the cumulative distribution function, via a Gram-Charlier 
#' expansion on the standardized distribution. This expansion uses
#' some weighting function, \eqn{w}{w}, typically the density of some \sQuote{parent}
#' probability distribution, and polynomials, \eqn{p_n}{p_n} which are
#' orthogonal with respect to that weighting:
#' \deqn{\int_{-\infty}^{\infty} p_n(x) p_m(x) w(x) \mathrm{d}x = h_n \delta_{mn}.}{integral_-inf^inf p_n(x) p_m(x) w(x) dx = h_n delta_mn.}
#'
#' Let \eqn{f(x)}{f(x)} be the probability density of some random variable,
#' with cumulative distribution function \eqn{F(x)}{F(x)}. We express
#' \deqn{f(x) = \sum_{n \ge 0} c_n p_n(x) w(x)}{f(x) = sum_n c_n p_n(x) w(x)}
#' The constants \eqn{c_n}{c_n} can be computed from the known moments
#' of the distribution.
#'
#' For the Gram Charlier \dQuote{A} series, the weighting function is the PDF of the
#' normal distribution, and the polynomials are the (probabilist's) Hermite
#' polynomials. As a weighting function, one can also use the PDF of the gamma
#' distribution (resulting in generalized Laguerre polynomials), or the
#' PDF of the Beta distribution (resulting in Jacobi polynomials).
#'
