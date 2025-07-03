# A few additional utilities to for forestry measurements

#' Calculate Quadratic Mean Diameter from a list of DBH measurements
#'
#' This function calculates quadratic mean diameter (qmd), and returns
#' a value in the same units as the inputs
#' @param dbh numeric - vector representing measure of diameter at breast height
#' @param na.rm boolean - should NAs be ignored.
#' @examples
#' library(landecoutils)
#' print(qmd(c(40.1, NA, 30.3, 29.4, 23.5), na.rm=TRUE))
#' @export
qmd = function(dbh, na.rm=TRUE) {
  if(na.rm) dbh = dbh[!is.na(dbh)]
  if(length(dbh)==0) return(NA)
  qmd = sqrt(sum(dbh^2)/length(dbh))
  return(qmd)
}


#' Calculate the single-trunk DBH equivalent of a multi-trunked stem.
#'
#' This function calculates the single-trunk equivalent DBH of a multi-trunked
#' stem by calculating the sum of the basal area of all stems, then calculuating the
#' single-tree DBH that would produce the same basal area. Following the method
#' of Steward and Salazar, 1992.
#'
#' J.L. Stewart, R. Salazar. A review of measurement options for multipurpose trees
#' Agrofor. Syst., 19 (1992), pp. 173-183, doi: 10.1007/BF00138507
#' @param dbh numeric - vector representing measure of mulit-trunked stem diameter at breast height
#' @param na.rm boolean - should NAs be ignored.
#' @examples
#' library(landecoutils)
#' print(dbh_equiv(c(10.1, 6.4, NA, 2.3), na.rm=TRUE))
#' @export
dbh_equiv = function(dbh, na.rm=TRUE) {
  if(na.rm) dbh = dbh[!is.na(dbh)]
  if(length(dbh)==0) return(NA)
  dbh_equiv = 2*sqrt(sum((dbh/2)^2))
  return(dbh_equiv)
}



#' Calculate basal area
#'
#' This function calculates basal area from DBH using simple circular assumption.
#' Out put is in squared units of inputs.
#' @param dbh numeric - vector representing measure of diameter at breast height
#' @param na.rm boolean - should NAs be ignored.
#' @examples
#' library(landecoutils)
#' print(ba(c(.101, .064, NA, 0.023), na.rm=TRUE))
#' @export
ba = function(dbh, na.rm = TRUE) {
  if(na.rm) dbh = dbh[!is.na(dbh)]
  ba = pi*(dbh/200)^2
  return(ba)
}
