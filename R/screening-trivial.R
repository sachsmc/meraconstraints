#' Identify constraints that are nontrivial at deterministic distributions
#'
#' Enumerates all deterministic (0/1) counterfactual distributions in each
#' district and checks which constraints are violated by at least one of them.
#' Constraints that are never violated are trivially implied by the
#' probability simplex. See [prob_simplex_screen()] for a faster equivalent.
#'
#' @param drels A list of relations as returned by [dhvm_linear_relations()]
#' @param hreps A list of H-representations, one for each district, as in the
#'   `Hrep` element returned by [dhvm_relations_to_polytope()]
#' @returns A list of logical vectors, one for each district, with one element
#'   per constraint that is `TRUE` if the constraint is possibly nontrivial and
#'   `FALSE` if it is implied by the probability simplex
#' @export
prob_simplex_extremes <- function(drels, hreps) {

  checks <- vector(mode = "list", length = length(hreps))
  for(i in 1:length(hreps)) {


    hrep <- hreps[[i]]

    if(length(drels$district_parms$poterms[[i]]$do) == 0) {
      checks[[i]] <- rep(FALSE, nrow(hrep))
      next
    }
    ctabin <- split(drels$problevels[[i]], drels$problevels[[i]][, drels$district_parms$poterms[[i]]$do])
    constgrps <- lapply(ctabin, \(cin) {
      # res <- c(glue::glue_data(cin,
      #                   paste(c("p", paste0("{", drels$poterms[[i]]$out, "}"),
      #                         "d", paste0("{", drels$poterms[[i]]$do, "}")),
      #                         collapse = "")))
      res <- do.call(paste, c(as.list(cin), sep = "."))
      res
    })

    rawprobs <- unlist(drels$probform[[i]])

    set1<- sapply(constgrps,"[", 1)
    combos <- combinations(rep(length(constgrps[[1]]), length(constgrps))) + 1
    #combos <- matrix(runif(length(set1) * 10000), ncol = length(set1))
    check <- matrix(NA, nrow = nrow(combos), ncol = nrow(hrep))
    for(j in 1:nrow(combos)) {
      theseprobs <- rep(0, length(rawprobs))
      names(theseprobs) <- names(rawprobs)
      for(k in 1:ncol(combos)) {
        theseprobs[constgrps[[k]][combos[j, k]]] <- 1
      }

      check[j, ] <- hrep_test(hrep, theseprobs, tol = 1e-8, which = TRUE)
      #allprobs[j, ] <- theseprobs

    }

    checks[[i]] <- colMeans(check) < 1
  }

  names(checks) <- names(hreps)
  checks
}

#' Screen constraints that are implied by the probability simplex
#'
#' For each constraint, computes its maximum and minimum over the product of
#' probability simplices (one per level of the intervened variables) and
#' checks whether it can be violated. Constraints that cannot be violated are
#' trivially implied by the probability simplex.
#'
#' @inheritParams prob_simplex_extremes
#' @param tol Numeric tolerance for the checks
#' @inherit prob_simplex_extremes return
#' @export
prob_simplex_screen <- function(drels, hreps, tol = 1e-8) {

  checks <- vector(mode = "list", length = length(hreps))
  for (i in seq_along(hreps)) {

    hrep <- hreps[[i]]

    if (length(drels$district_parms$poterms[[i]]$do) == 0) {
      checks[[i]] <- rep(FALSE, nrow(hrep))
      next
    }
    ctabin <- split(drels$problevels[[i]],
                    drels$problevels[[i]][, drels$district_parms$poterms[[i]]$do])
    blocks <- lapply(ctabin, \(cin) do.call(paste, c(as.list(cin), sep = ".")))
    pnames <- names(unlist(drels$probform[[i]]))
    blockidx <- lapply(blocks, match, pnames)

    A  <- -hrep[, -(1:2), drop = FALSE]
    b  <- hrep[, 2]
    eq <- hrep[, 1] == 1

    blockmax <- do.call(cbind, lapply(blockidx, \(ix) apply(A[, ix, drop = FALSE], 1, max)))
    blockmin <- do.call(cbind, lapply(blockidx, \(ix) apply(A[, ix, drop = FALSE], 1, min)))
    hi <- rowSums(blockmax)
    lo <- rowSums(blockmin)

    checks[[i]] <- ifelse(eq,
                          abs(hi - b) > tol | abs(lo - b) > tol,  # equality violated somewhere
                          hi > b + tol)                           # inequality violated somewhere
  }
  names(checks) <- names(hreps)
  checks
}
