#' From a probability form as a list create text
#'
#' @param form A prob.form object
#' @param values A data frame containing variables appearing in form and their values
#' @param letter An optional character string to put at the start of the probability, e.g., "p" or "Pr"
#' @returns A character vector
#' @export
#' @examples
#' make_prob_names(list(out = c("Y", "X"), cond = "Z"),
#' values = expand.grid(Y = 0:1, X = 0:1, Z = 0:1), letter = "Pr")
#'


make_prob_names <- function(form, values, letter = "P") {

  nameplace <- if(length(form$cond) > 0) {
    paste0(letter, "(", paste(paste0(form$out, " = %s"), collapse = ", ")," | ",
           paste(paste0(form$cond, " = %s"), collapse = ", "), ")"
    )
  } else {
    paste0(letter, "(", paste(paste0(form$out, " = %s"), collapse = ", "), ")"
    )
  }

  do.call(sprintf, c(nameplace, as.list(values[, unlist(form), drop = FALSE])))

}
