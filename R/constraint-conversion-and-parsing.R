#' Convert counterfactual relations to the H-representation of a polytope
#'
#' @param prob A list of probability forms
#' @param relation A list of relations for the corresponding probability forms
#' @returns An hrep object, see \link[rcdd]{scdd}
#' @export


cf_relations_to_hrep <- function(prob, relation) {
  ## make the matrix

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

    # rdex <- do.call(rbind, lapply(1:length(relation[[j]]), \(i) {
    #   rr <- relation[[j]][[i]]
    #   cbind(i, match(rr, rvarprobs))
    # }))

    #if(missing(thisr)){
    #  thisr <- rep("1", nrow(rdex))
    #}
    #Rmat[rdex] <- thisr[rdex[, 1]]

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

hrep_to_text <- function(hrep, prob, sets, simplify = FALSE) {

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

    termfactl <- (a[i,a[i,]!=0])
    termfact <- ifelse(termfactl < 0, paste0(" ", termfactl, "*"),
                       paste0(" + ", termfactl, "*"))
    termfact[abs(termfactl) == 1] <- gsub("1*", " ", termfact[abs(termfactl) == 1], fixed = TRUE)

    termfact[1] <- gsub("+", "", termfact[1], fixed = TRUE)
    inp <- paste(rbind(termfact, thisp), collapse = "")

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


prob_simplex_extremes <- function(drels, hreps) {

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


