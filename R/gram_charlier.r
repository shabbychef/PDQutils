# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of PDQutils.
#
# PDQutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PDQutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PDQutils.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.02.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# for the Hermite Polynomials
require(orthopolynom)
require(moments)

.guess_par <- function(raw.moments,gammapar=NULL) {
	# moments for gamma are kΘ and kΘ² + (kθ)²
	if (is.null(gammapar)) {
		theta <- (raw.moments[2]/raw.moments[1]) - raw.moments[1]
		k <- raw.moments[1] / theta
		gammapar <- list(shape=k,scale=theta)
	}
	gammapar
}

# here's the skinny. let g(x) be the pdf of the gamma RV w/ shape k and scale Θ.
# g(x) = (1 / (Γ(k) Θ)) * (x/Θ)^(k-1) * exp(-(x/Θ))
# Let L_n^(k-1) be the nth generalized Laguerre polynomial with alpha=k-1
# then
# integral_0^\inf L_n^(k-1)(x/Θ) L_m^(k-1)(x/Θ) g(x)dx = δ_n,m Γ(n+k) / (Γ(k) n!)

#' @title Approximate density and distribution via Gram-Charlier A expansion.
#'
#' @description 
#'
#' Approximate the probability density or cumulative distribution function of a distribution via its raw moments.
#'
#' @template details-gca
#'
#' @usage
#'
#' dapx_gca(x, raw.moments, support=c(-Inf,Inf), basis=c('normal','gamma'), gammapar=NULL, log=FALSE)
#'
#' papx_gca(q, raw.moments, support=c(-Inf,Inf), basis=c('normal','gamma'), gammapar=NULL, lower.tail=TRUE, log.p=FALSE)
#'
#' @param x where to evaluate the approximate density.
#' @param q where to evaluate the approximate distribution.
#' @param raw.moments an atomic array of the 1st through kth raw moments
#' of the probability distribution. 
#' @param support the support of the density function. It is assumed
#' that the density is zero on the complement of this open interval.
#' @param basis the basis under which to perform the approximation. \code{'normal'}
#' gives the classical expansion around the PDF and CDF of the normal
#' distribution via Hermite polynomials. \code{'gamma'} expands around a
#' gamma distribution with parameters \code{gammapar$shape} and
#' \code{gammapar$scale}.
#' @param gammapar the parameters for the gamma distribution approximation. Only
#' used under the gamma basis. If \code{NULL}, the shape and rate are inferred
#' from the first two moments.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail whether to compute the lower tail. If false, we approximate the survival function.
#' @return The approximate density at \code{x}, or the approximate CDF at
#' \code{q}.
#'
#' @keywords distribution 
#' @seealso \code{\link{qapx_cf}}
#' @export 
#' @template ref-Jaschke
#' @template ref-Blinnikov
#' @aliases papx_gca 
#' @note 
#'
#' Monotonicity of the CDF is not guaranteed.
#'
#' @examples 
#' # normal distribution:
#' xvals <- seq(-2,2,length.out=501)
#' d1 <- dapx_gca(xvals, c(0,1,0,3,0))
#' d2 <- dnorm(xvals)
#' d1 - d2
#'
#' qvals <- seq(-2,2,length.out=501)
#' p1 <- papx_gca(qvals, c(0,1,0,3,0))
#' p2 <- pnorm(qvals)
#' p1 - p2
#' @template etc
dapx_gca <- function(x,raw.moments,support=c(-Inf,Inf),basis=c('normal','gamma'),gammapar=NULL,log=FALSE) {#FOLDUP
	order.max <- length(raw.moments)
	basis <- match.arg(basis)

	# standardize the distribution 
	mu.central <- moments::raw2central(c(1,raw.moments))
	mu.std <- central2std(mu.central)

	# mean and standard deviation
	mu <- raw.moments[1]
	sigma <- sqrt(mu.central[3])

	if (basis == 'normal') {
		# the studentized variable:
		eta <- (x - mu) / sigma
		hermi <- orthopolynom::hermite.he.polynomials(order.max+1, normalized=FALSE)

		retval <- dnorm(eta)
		phi.eta <- retval
		for (iii in c(1:order.max)) {
			ci <- (sum(coef(hermi[[iii+1]]) * mu.std[1:(iii+1)])) / factorial(iii)
			retval <- retval + ci * phi.eta * as.function(hermi[[iii+1]])(eta)
		}

		# adjust back from standardized
		retval <- retval / sigma
	} else if (basis == 'gamma') {
		gammapar <- .guess_par(raw.moments,gammapar)

		# the studentized variable, for a chi-square
		scalby <- 1/gammapar$scale
		eta <- x * scalby

		alpha <- gammapar$shape - 1
		glag <- orthopolynom::glaguerre.polynomials(order.max+1, alpha, normalized=FALSE)

		retval <- dgamma(eta,shape=gammapar$shape,scale=1)
		phi.eta <- retval
		for (iii in c(1:order.max)) {
# 2FIX: the weighting all bad...
			ci <- (sum(coef(glag[[iii+1]]) * mu.std[1:(iii+1)])) / factorial(iii)
			retval <- retval + ci * phi.eta * as.function(glag[[iii+1]])(x)
		}

		# adjust back from standardized
		retval <- retval * scalby
	} else {
		stop('unknown basis')
	}

	# sanity check; shall I throw a warning?
	retval <- pmax(0,retval)

	# support support
	if (is.finite(min(support))) {
		retval[x <= min(support)] <- 0
	}
	if (is.finite(max(support))) {
		retval[x >= max(support)] <- 0
	}

	# must be a better way to do this ... 
	if (log)
		retval <- log(retval)
	return(retval)
}#UNFOLD
#' @export 
papx_gca <- function(q,raw.moments,support=c(-Inf,Inf),basis=c('normal','gamma'),gammapar=NULL,lower.tail=TRUE,log.p=FALSE) {#FOLDUP
	order.max <- length(raw.moments)
	basis <- match.arg(basis)

	# 2FIX: would it not be better to pass lower.tail to pnorm and dnorm below?
	# or subtract 1 from the end result?
	if (!lower.tail) {
		# transform q and the raw moments
		q <- - q;
		raw.moments <- raw.moments * (-1^(1:order.max))
	}

	mu.central <- moments::raw2central(c(1,raw.moments))
	mu.std <- central2std(mu.central)
	eta <- (q - raw.moments[1]) / sqrt(mu.central[3])
	hermi <- orthopolynom::hermite.he.polynomials(order.max, normalized=FALSE)

	retval <- pnorm(eta)
	phi.eta <- dnorm(eta)
	for (iii in c(1:order.max)) {
		ci <- (sum(coef(hermi[[iii+1]]) * mu.std[1:(iii+1)])) / factorial(iii)
		# n.b. the minus here!
		# and we are using He_{j-1} instead of He_j as in dapx_gca
		retval <- retval - ci * phi.eta * as.function(hermi[[iii]])(eta)
	}
	# sanity check; shall I throw a warning?
	retval <- pmin(1,pmax(0,retval))

	# support support
	if (is.finite(min(support))) {
		retval[q <= min(support)] <- 0
	}
	if (is.finite(max(support))) {
		retval[q >= max(support)] <- 1
	}

	# must be a better way to do this ... 
	if (log.p)
		retval <- log(retval)
	return(retval)
}#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
