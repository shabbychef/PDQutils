#' @details
#'
#' Given the raw cumulants of a probability distribution, we can approximate the probability 
#' density function, or the cumulative distribution function, via an Edgeworth
#' expansion on the standardized distribution. The derivation of the Edgeworth
#' expansion is rather more complicated than that of the Gram Charlier
#' approximation, involving the characteristic function and an expression of
#' the higher order derivatives of the composition of functions; see
#' Blinnikov and Moessner for more details. The Edgeworth expansion can
#' be expressed succinctly as
#' \deqn{\sigma f(\sigma x) = \phi(x) + \phi(x)\sum_{1 \le s}\sigma^s \sum_{\{k_m\}} He_{s+2r}(x) c_{k_m},}{sigma f(sigma x) = phi(x) + phi(x)  sum_{1 <= s} sigma^s sum_{k_m} He_{s+2r}(x) c_{k_m},}
#' where the second sum is over some partitions, and the constant \eqn{c} 
#' involves cumulants up to order \eqn{s+2}. Unlike the Gram Charlier
#' expansion, of which it is a rearrangement, the Edgeworth expansion
#' is arranged in increasing powers of the standard deviation
#' \eqn{\sigma}{sigma}.
#'
