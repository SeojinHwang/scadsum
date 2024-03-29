% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indeplasso.R
\name{indeplasso}
\alias{indeplasso}
\title{Independent LASSO using summary statistics (a.k.a. soft-thresholding)}
\usage{
indeplasso(coef, lambda = exp(seq(log(0.001), log(0.1), length.out = 20)))
}
\arguments{
\item{coef}{vector of regression coefficients (\eqn{r})}

\item{lambda}{a vector of \eqn{\lambda}s}
}
\value{
A list with the following
\item{lambda}{Same as \code{lambda} in input}
\item{beta}{A matrix of estimates of \eqn{\beta}}
}
\description{
Independent LASSO using summary statistics (a.k.a. soft-thresholding)
}
\details{
A function to find the minimum of \eqn{\beta} in  
\deqn{f(\beta)=\beta'\beta - 2\beta'r + 2\lambda||\beta||_1}
where \eqn{r} is the vector of regression coefficients. The analytical solution
is given by
\deqn{\hat{\beta}=sign(r)(max(|r| - \lambda))}
}
