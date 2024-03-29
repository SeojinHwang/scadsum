% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridgeR.R
\name{ridgeR}
\alias{ridgeR}
\title{Ridge using summary statistics}
\usage{
ridgeR(
  lambda1,
  lambda2 = 0,
  X,
  b,
  thr = 1e-04,
  trace = 0,
  maxiter = 10000,
  blocks = NULL,
  x = NULL
)
}
\description{
Coordinate descent algorithm to solve: 
0.5 x'X'Xx - x'b + 0.5 lambda1 ||x||_2^2 + 0.5 lambda2 ||x||_2^2
Function to get elastic net solutions given X, a reference panel, and
b, regression coefficients
}
\keyword{internal}
