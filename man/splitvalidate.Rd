% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitvalidate.R,
%   R/splitvalidate.lassosum.pipeline.R
\name{splitvalidate}
\alias{splitvalidate}
\alias{splitvalidate.lassosum.pipeline}
\title{Function to perform split-validation using output from lassosum.pipeline with external phenotype}
\usage{
splitvalidate(...)

\method{splitvalidate}{lassosum.pipeline}(
  ls.pipeline,
  test.bfile = NULL,
  keep = NULL,
  remove = NULL,
  pheno = NULL,
  covar = NULL,
  trace = 1,
  split = NULL,
  rematch = !is.null(test.bfile),
  ...
)
}
\arguments{
\item{...}{parameters to pass to \code{\link{validate.lassosum.pipeline}}}

\item{ls.pipeline}{A lassosum.pipeline object}

\item{test.bfile}{The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset}

\item{keep}{Participants to keep (see \code{\link{lassosum}} for more details)}

\item{remove}{Participants to remove}

\item{pheno}{A vector of phenotype OR a \code{data.frame} with 3 columns, the first 2 columns being headed "FID" and "IID", OR a filename for such a data.frame}

\item{covar}{A matrix of covariates OR a \code{data.frame} with 3 or more columns, the first 2 columns being headed "FID" and "IID", OR a filename for such a data.frame}

\item{trace}{Controls amount of output}

\item{rematch}{Forces a rematching of the ls.pipline beta's with the new .bim file}
}
\description{
Function to perform split-validation using output from lassosum.pipeline with external phenotype
}
\details{
Performs split-validation. Randomly split the test data into half for validation 
and half for prediction. Standardize the best cross-predicted pgs and stack together.
}
