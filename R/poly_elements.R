#' @title Number of \code{poly} elements.
#' @description This function aims calculate the number of terms of a polynomial interactions.
#' @param n The number of variables.
#' @param d Degreess of polynomial interaction.
#' @export
poly_elements <- function(n, d) {
  x <- sapply(1:d, combination_with_repetition, n = n) 
  return(sum(x))
}
