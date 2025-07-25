% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forestry_measures.R
\name{dbh_equiv}
\alias{dbh_equiv}
\title{Calculate the single-trunk DBH equivalent of a multi-trunked stem.}
\usage{
dbh_equiv(dbh, na.rm = TRUE)
}
\arguments{
\item{dbh}{numeric - vector representing measure of mulit-trunked stem diameter at breast height}

\item{na.rm}{boolean - should NAs be ignored.}
}
\description{
This function calculates the single-trunk equivalent DBH of a multi-trunked
stem by calculating the sum of the basal area of all stems, then calculuating the
single-tree DBH that would produce the same basal area. Following the method
of Steward and Salazar, 1992.
}
\details{
J.L. Stewart, R. Salazar. A review of measurement options for multipurpose trees
Agrofor. Syst., 19 (1992), pp. 173-183, doi: 10.1007/BF00138507
}
\examples{
library(landecoutils)
print(dbh_equiv(c(10.1, 6.4, NA, 2.3), na.rm=TRUE))
}
