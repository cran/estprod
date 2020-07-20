#' @title Combination with repetition.
#' @description From combinatorial math, this function aims calculates combinations with repetitions.
#' @param n The number of elements (variables).
#' @param r The size of the groups (degreess of the polynomial interaction).
#' @export
combination_with_repetition <- function(n, r){
  factorial(r+n-1)/(factorial(n-1)*factorial(r))
} 
