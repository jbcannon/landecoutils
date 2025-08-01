% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clip_lidR_ctgs.R
\name{check_for_lax}
\alias{check_for_lax}
\title{Check for and create a lax index from a directory of LAS files}
\usage{
check_for_lax(input, n_cores = 1, write_lax = TRUE)
}
\arguments{
\item{input}{can be a LASCatalog from \code{lidR} package or a path to a directory containing .LAS or .LAZ files to index}

\item{n_cores}{number of cores to create doSNOW cluster}

\item{write_lax}{indicates if .lax file should be written (\code{TRUE}), or only
checked for (\code{FALSE})}
}
\description{
Check for LAX index files and create them as necessary.  Indexing greatly
increases the processingAn earlier version of this function overwrote the
original LAS resulting in possible corruption of the original file if the
process was interrupted. This updated version uses temporary files to avoid
this problem.
}
\examples{
check_for_lax('E:/my/las/dir/', n_cores=4)
}
