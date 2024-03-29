# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Count number of lines in a text file
#' 
#' @param fileName Name of file
#' @keywords internal
#' 
countlines <- function(fileName) {
    .Call('_scadsum_countlines', PACKAGE = 'scadsum', fileName)
}

#' Multiply genotypeMatrix by a matrix
#' 
#' @param fileName location of bam file
#' @param N number of subjects 
#' @param P number of positions 
#' @param input the matrix
#' @param col_skip_pos which variants should we skip
#' @param col_skip which variants should we skip
#' @param keepbytes which bytes to keep
#' @param keepoffset what is the offset
#' @return an armadillo genotype matrix 
#' @keywords internal
#' 
multiBed3 <- function(fileName, N, P, input, col_skip_pos, col_skip, keepbytes, keepoffset, trace) {
    .Call('_scadsum_multiBed3', PACKAGE = 'scadsum', fileName, N, P, input, col_skip_pos, col_skip, keepbytes, keepoffset, trace)
}

#' Multiply genotypeMatrix by a matrix (sparse)
#' 
#' @param fileName location of bam file
#' @param N number of subjects 
#' @param P number of positions 
#' @param input the matrix
#' @param col_skip_pos which variants should we skip
#' @param col_skip which variants should we skip
#' @param keepbytes which bytes to keep
#' @param keepoffset what is the offset
#' @return an armadillo genotype matrix 
#' @keywords internal
#' 
multiBed3sp <- function(fileName, N, P, beta, nonzeros, colpos, ncol, col_skip_pos, col_skip, keepbytes, keepoffset, trace) {
    .Call('_scadsum_multiBed3sp', PACKAGE = 'scadsum', fileName, N, P, beta, nonzeros, colpos, ncol, col_skip_pos, col_skip, keepbytes, keepoffset, trace)
}

#' Performs lasso
#'
#' @param lambda1 lambda
#' @param lambda2 lambda
#' @param X genotype Matrix
#' @param r correlations
#' @param x beta coef
#' @param thr threshold 
#' @param yhat A vector
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @return conv
#' @keywords internal
#' 
elnet <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter) {
    .Call('_scadsum_elnet', PACKAGE = 'scadsum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter)
}

repelnet <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_repelnet', PACKAGE = 'scadsum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}

#' sj0715
#' Performs scad
#'
#' @param lambda1 lambda
#' @param lambda2 lambda
#' @parem gamma gamma
#' @param X genotype Matrix
#' @param r correlations
#' @param x beta coef
#' @param thr threshold 
#' @param yhat A vector
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @return conv
#' @keywords internal
#' 
scad <- function(lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter) {
    .Call('_scadsum_scad', PACKAGE = 'scadsum', lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter)
}

repscad <- function(lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_repscad', PACKAGE = 'scadsum', lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}

#' sj0715
#' Performs mcp
#'
#' @param lambda1 lambda
#' @param lambda2 lambda
#' @parem gamma gamma
#' @param X genotype Matrix
#' @param r correlations
#' @param x beta coef
#' @param thr threshold 
#' @param yhat A vector
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @return conv
#' @keywords internal
#' 
mcp <- function(lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter) {
    .Call('_scadsum_mcp', PACKAGE = 'scadsum', lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter)
}

repmcp <- function(lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_repmcp', PACKAGE = 'scadsum', lambda1, lambda2, gamma, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}

#' Performs ridge 
#'
#' @param lambda1 lambda
#' @param lambda2 lambda
#' @param X genotype Matrix
#' @param r correlations
#' @param x beta coef
#' @param thr threshold 
#' @param yhat A vector
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @return conv
#' @keywords internal
#' 
ridge <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter) {
    .Call('_scadsum_ridge', PACKAGE = 'scadsum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter)
}

repridge <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_repridge', PACKAGE = 'scadsum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}

#' imports genotypeMatrix
#' 
#' @param fileName location of bam file
#' @param N number of subjects 
#' @param P number of positions 
#' @param col_skip_pos which variants should we skip
#' @param col_skip which variants should we skip
#' @param keepbytes which bytes to keep
#' @param keepoffset what is the offset
#' @return an armadillo genotype matrix 
#' @keywords internal
#' 
genotypeMatrix <- function(fileName, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, fillmissing) {
    .Call('_scadsum_genotypeMatrix', PACKAGE = 'scadsum', fileName, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, fillmissing)
}

#' normalize genotype matrix
#' 
#' @param genotypes a armadillo genotype matrix
#' @return standard deviation
#' @keywords internal
#' 
normalize <- function(genotypes) {
    .Call('_scadsum_normalize', PACKAGE = 'scadsum', genotypes)
}

#' Runs lasso with various parameters
#' 
#' @param lambda1 a vector of lambdas (lambda2 is 0)
#' @param fileName the file name of the reference panel
#' @param r a vector of correlations
#' @param N number of subjects
#' @param P number of position in reference file
#' @param col_skip_posR which variants should we skip
#' @param col_skipR which variants should we skip
#' @param keepbytesR required to read the PLINK file
#' @param keepoffsetR required to read the PLINK file
#' @param thr threshold
#' @param x a numeric vector of beta coefficients
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @param Constant a constant to multiply the standardized genotype matrix
#' @return a list of results
#' @keywords internal
#'  
runElnet <- function(lambda, shrink, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_runElnet', PACKAGE = 'scadsum', lambda, shrink, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec)
}

#' sj0717
#' Runs scad with various parameters
#' 
#' @param lambda1 a vector of lambdas (lambda2 is 0)
#' @param gamma 
#' @param fileName the file name of the reference panel
#' @param r a vector of correlations
#' @param N number of subjects
#' @param P number of position in reference file
#' @param col_skip_posR which variants should we skip
#' @param col_skipR which variants should we skip
#' @param keepbytesR required to read the PLINK file
#' @param keepoffsetR required to read the PLINK file
#' @param thr threshold
#' @param x a numeric vector of beta coefficients
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @param Constant a constant to multiply the standardized genotype matrix
#' @return a list of results
#' @keywords internal
#'  
runScad <- function(lambda, shrink, gamma, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_runScad', PACKAGE = 'scadsum', lambda, shrink, gamma, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec)
}

#' sj0717
#' Runs mcp with various parameters
#' 
#' @param lambda1 a vector of lambdas (lambda2 is 0)
#' @param gamma 
#' @param fileName the file name of the reference panel
#' @param r a vector of correlations
#' @param N number of subjects
#' @param P number of position in reference file
#' @param col_skip_posR which variants should we skip
#' @param col_skipR which variants should we skip
#' @param keepbytesR required to read the PLINK file
#' @param keepoffsetR required to read the PLINK file
#' @param thr threshold
#' @param x a numeric vector of beta coefficients
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @param Constant a constant to multiply the standardized genotype matrix
#' @return a list of results
#' @keywords internal
#'  
runMcp <- function(lambda, shrink, gamma, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_runMcp', PACKAGE = 'scadsum', lambda, shrink, gamma, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec)
}

#' Runs ridge with various parameters
#' 
#' @param lambda1 a vector of lambdas (lambda2 is 0)
#' @param fileName the file name of the reference panel
#' @param r a vector of correlations
#' @param N number of subjects
#' @param P number of position in reference file
#' @param col_skip_posR which variants should we skip
#' @param col_skipR which variants should we skip
#' @param keepbytesR required to read the PLINK file
#' @param keepoffsetR required to read the PLINK file
#' @param thr threshold
#' @param x a numeric vector of beta coefficients
#' @param trace if >1 displays the current iteration
#' @param maxiter maximal number of iterations
#' @param Constant a constant to multiply the standardized genotype matrix
#' @return a list of results
#' @keywords internal
#'  
runRidge <- function(lambda, shrink, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec) {
    .Call('_scadsum_runRidge', PACKAGE = 'scadsum', lambda, shrink, fileName, r, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, thr, x, trace, maxiter, startvec, endvec)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_scadsum_RcppExport_registerCCallable', PACKAGE = 'scadsum')
})
