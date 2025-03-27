#' Convert counterfactual relations to the H-representation of a polytope
#'
#' @param prob A list of probability forms
#' @param relation A list of relations for the corresponding probability forms
#' @param thisr Multiplication factor
#' @returns An hrep object, see \link[rcdd]{scdd}
#' @export


cf_relations_to_hrep <- function(prob, relation, thisr) {
  ## make the matrix
  rvarprobs <- unique(unlist(relation))
  Rmat <- matrix("0", nrow = length(unlist(prob)), ncol = length(rvarprobs))

  rdex <- do.call(rbind, lapply(1:length(relation), \(i) {
    rr <- relation[[i]]
    cbind(i, match(rr, rvarprobs))
  }))
  if(missing(thisr)) thisr <- rep("1", nrow(rdex))
  Rmat[rdex] <- thisr[rdex[, 1]]

  vrep <- rcdd::makeV(points = t(Rmat))

  vrep <- rcdd::redundant(vrep)$output
  hrep2 <- rcdd::scdd(vrep, representation = "V")
  foo <- rcdd::q2d(hrep2$output)

  list(Vrep = vrep, Hrep = foo)

}


#' Parse Hrep and remove redundancies
#'
#' @param hrep H-representation
#' @param prob Character vector of probabilities
#' @param sets List of sets of probabilities that sum to 1
#' @export
#' @returns A character vector

hrep_to_text <- function(hrep, prob, sets, simplify = FALSE) {

  l <- hrep[ , 1]
  b <- hrep[ , 2]
  v <- hrep[ , - c(1, 2)]
  a <- -v
  # axb <- a %*% x - b
  # all(axb <= 0)
  # all(l * axb == 0)


  inequalities.text <- rep(NA_character_, nrow(a))
  for(i in 1:nrow(a)) {
    thisp <- prob[as.logical(a[i,] != 0)]
    # ## checking conditions for standard probabilistic constraints that we omit
    if(simplify) {
      if(l[i] == 0) {  ## inequalities
        if((sum(a[i,] == -1) == 1 & b[i] == 0) |
           (sum(a[i,] == 1) == 1 & b[i] == 1)) {
          next
        }
        if(any(sapply(sets, \(se) length(setdiff(thisp, se)) == 0))){
          if(all(a[i,][a[i,]!=0] == b[i]) & abs(b[i]) == 1) next
        }
      } else { ## equalities

        if(any(sapply(sets, \(se) length(setdiff(thisp, se)) == 0)) &
           (all(sign(a[i,][a[i,]!=0]) == -1) & b[i] == -1) |
           (all(sign(a[i,][a[i,]!=0]) == 1) & b[i] == 1)){
          next
        }
      }
    }

    #thisp <- prob[as.logical(abs(a[i,]))]

    inp <- paste(paste(a[i,a[i,]!=0], "*", thisp), collapse = " + ")

    inequalities.text[i] <- if(l[i] == 0) paste(inp, "<=", b[i]) else{
      paste(l[i], "*(", inp, " - ", b[i], ") = 0")
    }
  }

  inequalities.text[!is.na(inequalities.text)]

}

