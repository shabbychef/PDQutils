#' @details
#'
#' Given the raw cumulants of a probability distribution, we can approximate the probability 
#' density function, or the cumulative distribution function, via an Edgeworth
#' expansion on the standardized distribution.
#'
#' Suppose \eqn{f(x)}{f(x)} is the probability density of some random
#' variable, and let \eqn{F(x)}{F(x)} be the cumulative distribution function.
#' Let \eqn{He_j(x)}{He_j(x)} be the \eqn{j}{j}th probabilist's Hermite
#' polynomial. These polynomials form an orthogonal basis, with respect to the
#' function \eqn{w(x)}{w(x)} of the Hilbert space of functions which are square
#' \eqn{w}{w}-weighted integrable. The weighting functimn is 
#' \eqn{w(x) = e^{-x^2/2} = \sqrt{2\pi}\phi(x)}{w(x) = e^{-x^2/2} = sqrt(2pi) phi(x)}.
#' The orthogonality relationship is
#' \deqn{\int_{-\infty}^{\infty} He_i(x) He_j(x) w(x) \mathrm{d}x = \sqrt{2\pi} j! \delta_{ij}.}{integral_-inf^inf He_i(x) He_j(x) w(x) dx = sqrt(2pi)j!dirac_ij.}
#'
#' ... 2FIX: fill in here.
#'
