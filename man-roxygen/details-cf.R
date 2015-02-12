#' @details
#'
#' The Cornish Fisher approximation is the Legendre
#' inversion of the Edgeworth expansion of a distribution, but ordered
#' in a way that is convenient when used on the mean of a number of
#' independent draws of a random variable. 
#'
#' Suppose \eqn{x_1, x_2, \ldots, x_n}{x_1, x_2, ..., x_n} are \eqn{n} independent 
#' draws from some probability distribution. 
#' Letting 
#' \deqn{X = \frac{1}{\sqrt{n}} \sum_{1 \le i \le n} x_i,}{X = (x_1 + x_2 + ... x_n) / sqrt(n),}
#' the Central Limit Theorem assures us that, assuming finite variance, 
#' \deqn{X \rightarrow \mathcal{N}(\sqrt{n}\mu, \sigma),}{X ~~ N(sqrt(n) mu, sigma),}
#' with convergence in \eqn{n}.
#'
#' The Cornish Fisher approximation gives a more detailed picture of the
#' quantiles of \eqn{X}{X}, one that is arranged in decreasing powers of
#' \eqn{\sqrt{n}}{sqrt(n)}. The quantile function is the function \eqn{q(p)}{q(p)} 
#' such that \eqn{P\left(X \le q(p)\right) = q(p)}{P(x <= q(p)) = p}. The
#' Cornish Fisher expansion is 
#' \deqn{q(p) = \sqrt{n}\mu + \sigma \left(z + \sum_{3 \le j} c_j f_j(z)\right),}{q(p) = sqrt{n}mu + sigma (z + sum_{3 <= j} c_j f_j(z)),}
#' where \eqn{z = \Phi^{-1}(p)}{z = qnorm(p)}, and \eqn{c_j}{c_j} involves
#' standardized cumulants of the distribution of \eqn{x_i}{x_i} of order
#' up to \eqn{j}. Moreover, the \eqn{c_j}{c_j} include decreasing
#' powers of \eqn{\sqrt{n}}{sqrt(n)}, giving some justification for truncation.
#' When \eqn{n=1}{n=1}, however, the ordering is somewhat arbitrary.
#'
