#' @title path to included examples substances.xlsx file
#'
#' @return path to substances.xlsx file
#'
#' @export
substances.path <- function() system.file("extdata", "substances.xlsx", package = "PesticideLoadIndicator", mustWork=TRUE);

#' @title load included substances.xlsx file
#'
#' @return substances.xlsx file as data.frame
#'
#' @export
substances.load <- function() read.excel(substances.path());

#' @title path to  includedexamples products.xlsx file
#'
#' @return path to products.xlsx file
#'
#' @export
products.path <- function() system.file("extdata", "products.xlsx", package = "PesticideLoadIndicator", mustWork=TRUE);

#' @title load included products.xlsx file
#'
#' @return products.xlsx file as data.frame
#'
#' @export
products.load <- function() read.excel(products.path());
