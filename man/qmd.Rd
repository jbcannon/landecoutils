% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forestry_measures.R
\name{qmd}
\alias{qmd}
\title{Calculate Quadratic Mean Diameter from a list of DBH measurements}
\usage{
qmd(dbh, na.rm = TRUE)
}
\arguments{
\item{dbh}{numeric - vector representing measure of diameter at breast height}

\item{na.rm}{boolean - should NAs be ignored.}
}
\description{
This function calculates quadratic mean diameter (qmd), and returns
a value in the same units as the inputs
}
\examples{
library(landecoutils)
print(qmd(c(40.1, NA, 30.3, 29.4, 23.5), na.rm=TRUE))
}
