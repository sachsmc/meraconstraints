#' Convert counterfactual relations to the H-representation of a polytope
#'
#' @param relations A list of relations for the corresponding probability forms as returned by [dag_get_linear_relations()]
#' @returns A list with the v-representation and the h-representation objects, one for each district, see \link[rcdd]{scdd}
#' @export


dhvm_relations_to_polytope <- function(relations) {
  ## make the matrix

  prob <- relations$probform
  relation <- relations$relation
  foos <- vreps <- vector(mode = "list", length = length(prob))
  for(j in 1:length(prob)) {
    rvarprobs <- unique(unlist(relation[[j]]))
    Rmat <- matrix("0", nrow = length(unlist(prob[[j]])), ncol = length(rvarprobs))

    rRm <- 1
    for(k in 1:length(relation[[j]])) {
      rr <- relation[[j]][[k]]
      rrm <- match(rr, rvarprobs)
      ri <- seq(rRm, rRm + length(prob[[j]][[k]]) - 1)
      rdexi <- cbind(rep(ri, each = length(rrm)), rep(rrm, length(ri)))
      Rmat[rdexi] <- "1"
      rRm <- rRm + length(prob[[j]][[k]])
    }



    vreps[[j]] <- rcdd::makeV(points = t(Rmat))

    #vrep <- rcdd::redundant(vrep)$output
    hrep2 <- rcdd::scdd(vreps[[j]], representation = "V")
    foos[[j]] <- rcdd::q2d(hrep2$output)
  }

  list(Vrep = vreps, Hrep = foos)

}


#' Parse Hrep and remove redundancies
#'
#' @param hrep H-representation
#' @param prob Character vector of probabilities
#' @param sets List of sets of probabilities that sum to 1
#' @export
#' @returns A character vector

hrep_to_text <- function(hrep, prob) {

  l <- hrep[ , 1]
  b <- hrep[ , 2]
  v <- hrep[ , - c(1, 2)]
  a <- -v
  # axb <- a %*% x - b
  # all(axb <= 0)
  # all(l * axb == 0)
  if(is.list(prob)) prob <- unlist(prob)

  inequalities.text <- rep(NA_character_, nrow(a))
  for(i in 1:nrow(a)) {
    thisp <- prob[as.logical(a[i,] != 0)]

    #thisp <- prob[as.logical(abs(a[i,]))]

    termfactl <- (a[i,a[i,]!=0])
    termfact <- ifelse(termfactl < 0, paste0(" ", termfactl, "* "),
                       paste0(" +", termfactl, "* "))
    termfact[abs(termfactl) == 1] <- gsub("1*", "", termfact[abs(termfactl) == 1], fixed = TRUE)

    termfact[1] <- gsub(" - ", "-", gsub("+", "", termfact[1], fixed = TRUE))
    inp <- trimws(paste(rbind(termfact, thisp), collapse = ""))

    inequalities.text[i] <- if(l[i] == 0) paste(inp, "<=", b[i]) else{
      if(l[i] == 1) {
        paste(inp, " = ", c("0" = "0", "-1" = "-1", "1" = "1")[as.character(b[i])])
      } else {
      paste(l[i], "*(", inp, " - ", b[i], ") = 0")
      }
    }
  }

  inequalities.text[!is.na(inequalities.text)]

}


#' @export

hrep_test <- function(Hrep, probform, tol = 1e-10, which = FALSE) {

  l <- Hrep[ , 1]
  b <- Hrep[ , 2]
  v <- Hrep[ , - c(1, 2)]
  a <- (- v)

  axb <- a %*% probform - b

  chk <- axb <= tol & abs(l * axb) <= tol
  if(which) {
    chk
  } else {
    all(chk)
  }

}

#' Check if two Hreps are the same
#'
#' Lexigraphical sort the matrices, and then compare them elementwise. Is symmetric.
#'
#' @param H1 first hrep
#' @param H2 second hrep
#' @returns Logical vector of length 1
#' @export
#' @examples
#' ivhrep <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
#' 0, 1, 1, 0, 1, 0, 1, 0, 0, -1, -1, 0, 1, 0, 1, 0, -1, 0, 0, 0,
#' -1, 0, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, 0, 0,
#' 1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
#' 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0,
#' 1, 0, -1, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
#' -1, 0, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#' 0, 0, 1), dim = c(14L, 10L), representation = "H")
#'
#' ivd <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
#' ivrels2 <- dhvm_linear_relations(ivd)
#' hrepiv2 <- cf_relations_to_hrep(ivrels2)
#' hreps_are_equal(hrepiv2$Hrep[[1]], ivhrep)
#'
#'
hreps_are_equal <- function(H1, H2) {

  if(!identical(dim(H1), dim(H2))) {
    return(FALSE)
  }

  ## lexi sort and compare
  all( H1[do.call(order, apply(H1, 2, \(x) x, simplify = FALSE)), ] ==
    H2[do.call(order, apply(H2, 2, \(x) x, simplify = FALSE)), ] )


}


prob_simplex_extremes0 <- function(drels, hreps) {

  checks <- vector(mode = "list", length = length(hreps))
  for(i in 1:length(hreps)) {
  sepprobs <- lapply(strsplit(unlist(drels$probform[[i]]), ")", fixed = TRUE), \(x) paste0(x, ")"))
  rawprobs <- sort(unique(unlist(sepprobs)))

  # probmap1 <- lapply(sepprobs, \(prc) {
  #   match(prc, rawprobs)
  # })

  constgrps <- list()
  hrep <- hreps[[i]]

  for(j in 1:length(drels$problist[[i]])) {
    ptabin <- unique(drels$problevels[[i]][, c(names(drels$problist[[i]])[j], drels$problist[[i]][[j]])])
    y <- names(drels$problist[[i]])[j]
    ctabin <- split(ptabin, ptabin[, drels$problist[[i]][[j]]])
    constgrps <- c(constgrps, lapply(ctabin, \(cin) {
      c(glue::glue_data(cin, sprintf("P(%s = {%s}|%s)", y, y,
                                     paste(paste0(drels$problist[[i]][[j]],
                                                  " = {", drels$problist[[i]][[j]], "}"), collapse = ","))))
    }))
  }


  set1<- sapply(constgrps,"[", 1)
  combos <- combinations(rep(2, length(set1)))
  #combos <- matrix(runif(length(set1) * 10000), ncol = length(set1))
  check <- matrix(NA, nrow = nrow(combos), ncol = nrow(hrep))
  for(j in 1:nrow(combos)) {
      theseprobs <- rep(NA_integer_, length(rawprobs))
      names(theseprobs) <- rawprobs
      theseprobs[set1] <- combos[j,]
      for(kk in constgrps) {
        for(k in kk[-1]) {
          theseprobs[k] <- 1 - theseprobs[kk[1]]
        }
      }

      theseoutprobs <- sapply(sepprobs,\(vv) {
        prod(theseprobs[vv])
      })
      names(theseoutprobs) <- unlist(drels$probform[[i]])

      check[j, ] <- hrep_test(hrep, theseoutprobs, tol = 1e-8, which = TRUE)


    }

  checks[[i]] <- colMeans(check) < 1
  }

  checks
}



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

  checks
}


#' @export
combinations <- function(p) {

  tot = prod(p)
  cp = cumprod(p)
  retval = rep.int(0, tot * length(p))
  dim(retval) = c(tot, length(p))
  for (i in seq_along(p)) {
    retval[, i] <- rep(seq_len(p[i]) - 1, each = cp[i]/p[i])
  }
  retval

}


