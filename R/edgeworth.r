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

# Created: 2015.02.19
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# for the Hermite Polynomials
require(orthopolynom)
require(moments)

# blinnikov and moessner's ADVANCE function
advance <- function(kms) {
	current <- kms$current
	mold <- kms$mold

	n <- length(current)
	ords <- 1:(length(current))
	stopifnot(n == sum(current * ords))
	sumcur <- n
	m <- 1
	is.done <- FALSE
	while (!is.done) {
		sumcur <- sumcur - current[m]*m + (m+1)
		current[m] <- 0
		current[m+1] <- current[m+1] + 1
		m <- m+1
		is.done <- (sumcur <= n) || (m > mold)
	}
	mold <- max(mold,m)
	current[1] <- n - sumcur
	retv <- list(current=current,mold=mold)
	return(retv)
}

#' @title Approximate density and distribution via Edgeworth expansion.
#'
#' @description 
#'
#' Approximate the probability density or cumulative distribution function of a distribution via its raw cumulants.
#'
#' @template details-edgeworth 
#'
#' @usage
#'
#' dapx_edgeworth(x, raw.cumulants, support=c(-Inf,Inf), log=FALSE)
#'
#' papx_edgeworth(q, raw.cumulants, support=c(-Inf,Inf), lower.tail=TRUE, log.p=FALSE)
#'
#' @param x where to evaluate the approximate density.
#' @param q where to evaluate the approximate distribution.
#' @param raw.cumulants an atomic array of the 1st through kth raw cumulants
#' of the probability distribution. The first cumulant is the mean, the
#' second is the variance. The third is \emph{not} the typical unitless skew.
#' @param support the support of the density function. It is assumed
#' that the density is zero on the complement of this open interval.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail whether to compute the lower tail. If false, we approximate the survival function.
#' @return The approximate density at \code{x}, or the approximate CDF at
#' \code{q}.
#'
#' @keywords distribution 
#' @seealso \code{\link{dapx_gca}, \link{papx_gca}}
#' @export 
#' @template ref-Blinnikov
#' @aliases papx_edgeworth 
#' @note 
#'
#' Monotonicity of the CDF is not guaranteed.
#'
#' @examples 
#' # normal distribution, for which this is silly
#' xvals <- seq(-2,2,length.out=501)
#' d1 <- dapx_edgeworth(xvals, c(0,1,0,0,0,0))
#' d2 <- dnorm(xvals)
#' d1 - d2
#'
#' qvals <- seq(-2,2,length.out=501)
#' p1 <- papx_edgeworth(qvals, c(0,1,0,0,0,0))
#' p2 <- pnorm(qvals)
#' p1 - p2
#' @template etc
dapx_edgeworth <- function(x,raw.cumulants,support=c(-Inf,Inf),log=FALSE) {#FOLDUP
	order.max <- length(raw.cumulants)

	mu <- raw.cumulants[1]
	sigma <- sqrt(raw.cumulants[2])

	eta <- (x - mu) / sigma

	phi.eta <- dnorm(eta)
	retval <- phi.eta

	if (order.max > 2) {
		hermi <- orthopolynom::hermite.he.polynomials(3*order.max+2, normalized=FALSE)
		Sn <- raw.cumulants[3:order.max] / (raw.cumulants[2] ^ (2:(order.max-1)))

		for (s in c(1:(order.max-2))) {
			nexterm <- rep(0,length(phi.eta))
			kms <- list(mold=1,current=c(s,rep(0,s-1)))
			while (kms$mold <= s) {
				r <- sum(kms$current)
				coefs <- ((Sn[1:s] / factorial(3:(s+2))) ^ kms$current) / factorial(kms$current)
				coef <- prod(coefs)
				nexterm <- nexterm + coef * as.function(hermi[[s+2*r+1]])(eta)
				kms <- advance(kms)
			}
			retval <- retval + phi.eta * (sigma^s) * nexterm
		}
	}

	# adjust back from standardized
	retval <- retval / sigma

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
papx_edgeworth <- function(q,raw.cumulants,support=c(-Inf,Inf),lower.tail=TRUE,log.p=FALSE) {#FOLDUP
	order.max <- length(raw.cumulants)

	# 2FIX: would it not be better to pass lower.tail to pnorm and dnorm below?
	# or subtract 1 from the end result?
	if (!lower.tail) {
		# transform q and the raw cumulants
		q <- - q;
		raw.cumulants <- raw.cumulants * (-1^(1:order.max))
	}

	mu <- raw.cumulants[1]
	sigma <- sqrt(raw.cumulants[2])

	eta <- (q - mu) / sigma

	phi.eta <- dnorm(eta)
	retval <- pnorm(eta)

	if (order.max > 2) {
		hermi <- orthopolynom::hermite.he.polynomials(3*order.max+2, normalized=FALSE)
		Sn <- raw.cumulants[3:order.max] / (raw.cumulants[2] ^ (2:(order.max-1)))

		for (s in c(1:(order.max-2))) {
			nexterm <- rep(0,length(phi.eta))
			kms <- list(mold=1,current=c(s,rep(0,s-1)))
			while (kms$mold <= s) {
				r <- sum(kms$current)
				coefs <- ((Sn[1:s] / factorial(3:(s+2))) ^ kms$current) / factorial(kms$current)
				coef <- prod(coefs)
				# n.b. the hermite polynomial is one order less than in the dapx
				# and we _subtract_ it
				nexterm <- nexterm - coef * as.function(hermi[[s+2*r]])(eta)
				kms <- advance(kms)
			}
			retval <- retval + phi.eta * (sigma^s) * nexterm
		}
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

## test nakagami#FOLDUP
##q <- seq(-2,2,length.out=1001)
##pv <- papx_edg(q,raw.cumulants=c(0,1,0,1,1,1))
##plot(q,pv)
##lines(q,pnorm(q),col='red')
###plot(x,dnorm(x) - dv)

## compute the first ord.max raw moments of the Nakagami distribution
#nak_moments <- function(m,omega=1,ord.max=10) {
	#ords <- 1:ord.max
	#df <- 2*m
	#moms <- 2 ^ (ords/2.0) * exp(lgamma((df+ords)/2) - lgamma(df/2))
	#moms <- moms * ((omega/df) ^ (ords/2.0))
	#return(moms)
#}
#nak_cumulants <- function(...) {
	#cumulants <- moment2cumulant(nak_moments(...))
	#return(cumulants)
#}


## perhaps the computation of cumulants has
## numerical issues?
#mom <- nak_moments(m=30,omega=150,ord.max=20)
#rac <- moment2cumulant(mom)

#momc <- moments::raw2central(c(1,mom))
#momc <- momc[2:length(momc)]
#racc <- moment2cumulant(momc)
#racc[1] <- mom[1]


#q <- seq(0,20,length.out=1001)
#pv <- papx_edg(q,raw.cumulants=racc,support=c(0,Inf))
#plot(q,pv)
#pv2 <- papx_gca(q,mom,support=c(0,Inf))
#lines(q,pv2,col='red')

##dv <- dapx_edg(q,raw.cumulants=rac,support=c(0,Inf))
##plot(q,dv)
##dv2 <- dapx_gca(q,mom,support=c(0,Inf))
##lines(q,dv2,col='red')

### symmetry trick ...
### does not help...
##mom2 <- rac
##mom2[seq(1,length(mom2),by=2)] <- 0
##rac2 <- moment2cumulant(mom2)
##dv3 <- 0.5 * dapx_edg(q,raw.cumulants=rac2,support=c(9,Inf))
##lines(q,dv3,col='blue')
##UNFOLD

## checking#FOLDUP
## the problem does _not_ seem to be numerical issues
## in converting moments to cumulants.
#lambda <- 0.7
#max.ord <- 15
#ords <- 1:max.ord
#moms <- exp(lfactorial(ords) - ords * log(lambda))
## roundoff issues?
#rac1 <- moment2cumulant(moms)
## the 'exact'
#rac <- exp(lfactorial(ords-1) - ords * log(lambda))
#err <- rac1 - rac
#rerr <- err / rac
##' cumul <- (rate ^ -(1:n)) * factorial(0:(n-1))

#q <- seq(0,20,length.out=1001)
#truv <- pexp(q,rate=lambda)
#plot(q,truv)
#pv <- papx_edg(q,raw.cumulants=rac,support=c(0,Inf))
#lines(q,pv,col='blue')
#pv1 <- papx_edg(q,raw.cumulants=rac1,support=c(0,Inf))
#lines(q,pv1,col='green')
#pv2 <- papx_gca(q,moms,support=c(0,Inf))
#lines(q,pv2,col='red')

#plot(q,pv-pv1)
##UNFOLD

## checking#FOLDUP
## the problem does _not_ seem to be numerical issues
## in converting moments to cumulants.
#lambda <- 0.7
#max.ord <- 6
#ords <- 1:max.ord
#moms <- exp(lfactorial(ords) - ords * log(lambda))
#rac <- exp(lfactorial(ords-1) - ords * log(lambda))

#q <- seq(0,20,length.out=1001)
#truv <- pexp(q,rate=lambda)
#plot(q,truv)
#pv <- papx_edg(q,raw.cumulants=rac,support=c(0,Inf))
#lines(q,pv,col='green')
#pv2 <- papx_gca(q,moms,support=c(0,Inf))
#lines(q,pv2,col='red')
##UNFOLD

## checking#FOLDUP
## the problem does _not_ seem to be numerical issues
## in converting moments to cumulants.
#lambda <- 0.7
#max.ord <- 4
#ords <- 1:max.ord
#moms <- exp(lfactorial(ords) - ords * log(lambda))
#rac <- exp(lfactorial(ords-1) - ords * log(lambda))

#x <- seq(0,20,length.out=1001)
#truv <- dexp(x,rate=lambda)
#plot(x,truv)
#dv <- dapx_edg(x,raw.cumulants=rac,support=c(0,Inf))
#lines(x,dv,col='green')
#dv2 <- dapx_gca(x,moms,support=c(0,Inf))
#lines(x,dv2,col='red')
##UNFOLD

## look at a root exponential distribution#FOLDUP
## with density
## lambda |y| exp(-lambda *y^2)
## this is a density,
## and y^2 is exponential(lambda)

#lambda <- 0.7
#max.ord <- 7
#ords <- 1:max.ord
#moms <- rep(0,2*max.ord)
#moms[seq(2,length(moms),by=2)] <- exp(lfactorial(ords) - ords * log(lambda))
#rac <- moment2cumulant(moms)

#x <- seq(-4,4,length.out=1001)
#truv <- lambda * abs(x) * exp(-lambda * x^2)
#plot(x,truv)

#dv <- dapx_edg(x,raw.cumulants=rac,support=c(-Inf,Inf))
#lines(x,dv,col='green')
#dv2 <- dapx_gca(x,moms,support=c(-Inf,Inf))
#lines(x,dv2,col='red')

## and CDF
#lambda <- 0.7
#max.ord <- 10
#ords <- 1:max.ord
#moms <- rep(0,2*max.ord)
#moms[seq(2,length(moms),by=2)] <- exp(lfactorial(ords) - ords * log(lambda))
#rac <- moment2cumulant(moms)

#q <- seq(-3,3,length.out=1001)
#truv <- rep(0,length(q))
#truv <- 0.5 * (1 + sign(q) * pexp(q^2,rate=lambda))
#plot(q,truv)

#pv <- papx_edg(q,raw.cumulants=rac,support=c(-Inf,Inf))
#lines(q,pv,col='green')
#pv2 <- papx_gca(q,moms,support=c(-Inf,Inf))
#lines(q,pv2,col='red')
##UNFOLD

## look at chisq distribution#FOLDUP

## match blinnikov and moessner
#df <- 5

#lambda <- 0.7
#max.ord <- 12 
#ords <- 1:max.ord
#subords <- ords - 1

#rac <- df * (2^subords) * factorial(subords)
#moms <- cumulant2moment(rac)

#x <- seq(0,25,length.out=1001)
#truv <- dchisq(x,df=df)

#plotx <- (x - df) / sqrt(2*df)
#plot(plotx,truv)

#dv <- dapx_edg(x,raw.cumulants=rac,support=c(0,Inf))
#lines(plotx,dv,col='green')
#dv2 <- dapx_gca(x,moms,support=c(0,Inf))
#lines(plotx,dv2,col='red')

## and CDF
#lambda <- 0.7
#max.ord <- 10
#ords <- 1:max.ord
#moms <- rep(0,2*max.ord)
#moms[seq(2,length(moms),by=2)] <- exp(lfactorial(ords) - ords * log(lambda))
#rac <- moment2cumulant(moms)

#q <- seq(-3,3,length.out=1001)
#truv <- rep(0,length(q))
#truv <- 0.5 * (1 + sign(q) * pexp(q^2,rate=lambda))
#plot(q,truv)

#pv <- papx_edg(q,raw.cumulants=rac,support=c(-Inf,Inf))
#lines(q,pv,col='green')
#pv2 <- papx_gca(q,moms,support=c(-Inf,Inf))
#lines(q,pv2,col='red')
##UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
