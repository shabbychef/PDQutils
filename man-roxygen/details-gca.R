#' @details
#'
#' Given the raw moments of a probability distribution, we can approximate the probability 
#' density function, or the cumulative distribution function, via a Gram-Charlier A 
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
#' Expanding the density \eqn{f(x)}{f(x)} in terms of these polynomials in the
#' usual way (abusing orthogonality) one has
#' \deqn{f(x) = \sum_{0\le j} \frac{He_j(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{f(x) = sum_{0 <= j} (He_j(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#' The cumulative distribution function is 'simply' the integral of this
#' expansion. Abusing certain facts regarding the PDF and CDF of the normal
#' distribution and the probabilist's Hermite polynomials, the CDF has
#' the representation
#' \deqn{F(x) = \Phi(x) - \sum_{1\le j} \frac{He_{j-1}(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{F(x) = Phi(x) - sum_{1 <= j} (He_{j-1}(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#'
#' These series contain coefficients defined by the probability distribution 
#' under consideration. They take the form
#' \deqn{c_j = \frac{1}{j!}\int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{c_j = (1/j!) integral_-inf^inf f(z) He_j(z) dz.}
#' Using linearity of the integral, these coefficients are easily computed in
#' terms of the coefficients of the Hermite polynomials and the raw, uncentered
#' moments of the probability distribution under consideration. Note that it may be the
#' case that the computation of these coefficients suffers from bad numerical
#' cancellation for some distributions, and that an alternative formulation
#' may be more numerically robust.
#'
