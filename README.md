

# PDQutils

Approximating Probability Distributions.

-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package may be installed from CRAN; the latest version may be
found on [github](https://www.github.com/shabbychef/PDQutils "PDQutils")
via devtools:


```r
if (require(devtools)) {
    # latest greatest
    install_github("shabbychef/PDQutils")
}
```

# Basic Usage

Approximating the distribution of a random variable via the Gram Charlier or Cornish Fisher
expansions is most convenient when the random variable can be decomposed as the sum of a 
small number of independent random variables whose cumulants can be computed. For example, 
suppose $Y = \sum_{1 \le i \le k} \sqrt{X_i / \nu_i}$ where the $X_i$ are independent central 
chi-square random variables with degrees of freedom $\nu_1,\nu_2,...,\nu_k$. I will call this
a 'snak' distribution, for 'sum of Nakagami', since each summand follows a 
[Nakagami distribution](https://en.wikipedia.org/wiki/Nakagami_distribution "Nakagami distribution").
We can easily write code that generates variates from this distribution given a vector
of the degrees of freedom:


```r
rsnak <- function(n, dfs) {
    samples <- Reduce("+", lapply(dfs, function(k) {
        sqrt(rchisq(n, df = k)/k)
    }))
}
```

Let's take one hundred thousand draws from this distribution and see whether it is approximately normal:


```r
n.samp <- 1e+05
dfs <- c(8, 15, 4000, 10000)
set.seed(18181)
# now draw from the distribution
rvs <- rsnak(n.samp, dfs)
qqnorm(rvs)
qqline(rvs, col = "red")
```

<img src="github_extra/figure/testit-1.png" title="plot of chunk testit" alt="plot of chunk testit" width="600px" height="500px" />

Using the additivity
property of cumulants, we can compute the cumulants of $Y$ easily if we have the cumulants of
the $X_i$. These in turn can be computed from the raw moments.  See
[wikipedia](https://en.wikipedia.org/wiki/Chi_distribution "chi distribution") for the raw moments
of the Chi distribution. The following function computes the cumulants:


```r
# for the moment2cumulant function:
library(PDQutils)
# compute the first ord.max raw cumulants of the
# sum of Nakagami variates
snak_cumulants <- function(dfs, ord.max = 10) {
    # first compute the raw moments
    moms <- lapply(dfs, function(k) {
        ords <- 1:ord.max
        moms <- 2^(ords/2) * exp(lgamma((k + ords)/2) - 
            lgamma(k/2))
        # we are dividing the chi by sqrt the d.f.
        moms <- moms/(k^(ords/2))
        moms
    })
    # turn moments into cumulants
    cumuls <- lapply(moms, moment2cumulant)
    # sum the cumulants
    tot.cumul <- Reduce("+", cumuls)
    return(tot.cumul)
}
```

We can now implement the 'dpq' functions trivially using the Gram-Charlier and Cornish-Fisher
approximations, as follows:


```r
library(PDQutils)

dsnak <- function(x, dfs, ord.max = 10, ...) {
    raw.moment <- cumulant2moment(snak_cumulants(dfs, 
        ord.max))
    retval <- dapx_gca(x, raw.moment, ...)
    return(retval)
}
psnak <- function(q, dfs, ord.max = 10, ...) {
    raw.moment <- cumulant2moment(snak_cumulants(dfs, 
        ord.max))
    retval <- papx_gca(q, raw.moment, ...)
    return(retval)
}
qsnak <- function(p, dfs, ord.max = 10, ...) {
    raw.cumul <- snak_cumulants(dfs, ord.max)
    retval <- qapx_cf(p, raw.cumul, ...)
    return(retval)
}
```

The qqplot should look better now:


```r
qqplot(qsnak(ppoints(n.samp), dfs = dfs), rvs, main = "Q-Q against qsnak (C-F appx.)")
qqline(rvs, distribution = function(p) qsnak(p, dfs), 
    col = "red")
```

<img src="github_extra/figure/improvedqq-1.png" title="plot of chunk improvedqq" alt="plot of chunk improvedqq" width="600px" height="500px" />


