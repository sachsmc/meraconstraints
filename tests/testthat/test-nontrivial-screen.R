## Exact screen of nontrivial constraints (Proposition 6 of the paper)

screen_model <- function(g) {
  rels <- dhvm_linear_relations(g)
  poly <- dhvm_relations_to_polytope(rels)
  list(rels = rels, poly = poly, nt = nontrivial_screen(g, rels, poly))
}

test_that("Verma graph has nontrivial equality and inequality constraints", {

  vm <- screen_model(dhvm(caugi(A %-->% B %-->% C %-->% D, U %-->% B + D)))
  H <- vm$poly$Hrep[["B_D"]]
  fl <- vm$nt[["B_D"]]
  expect_equal(sum(fl & H[, 1] == 0), 8)
  expect_equal(sum(fl & H[, 1] == 1), 2)
  expect_equal(sum(vm$nt[["C"]]), 0)

})

test_that("districts without external parents have no nontrivial constraints", {

  vm <- screen_model(dhvm(caugi(U %-->% X + Y)))
  expect_equal(sum(unlist(vm$nt)), 0)

})

test_that("inequality screen does not depend on how the rows are written", {

  sm <- screen_model(dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)))
  H <- sm$poly$Hrep[[1]]
  eq <- H[, 1] == 1

  ## add random combinations of the equality rows to the inequality rows, and rescale
  set.seed(1)
  H2 <- H
  lam <- matrix(sample(-2:2, sum(!eq) * sum(eq), replace = TRUE), sum(!eq))
  H2[!eq, -1] <- sample(1:3, sum(!eq), replace = TRUE) * (H[!eq, -1] + lam %*% H[eq, -1])
  poly2 <- sm$poly
  poly2$Hrep[[1]] <- H2
  poly2$HrepQ[[1]] <- rcdd::d2q(H2)

  expect_equal(nontrivial_screen(dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
                                 sm$rels, poly2)[[1]][!eq],
               sm$nt[[1]][!eq])
  ## the simplex screen is not invariant
  expect_false(isTRUE(all.equal(prob_simplex_screen(sm$rels, poly2$Hrep)[[1]][!eq],
                                prob_simplex_screen(sm$rels, sm$poly$Hrep)[[1]][!eq])))

})

test_that("inequality screen agrees with linear programming", {

  ## a row is nontrivial iff its maximum over the simplices intersected with the
  ## equality constraints exceeds its bound
  models <- list(
    dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
    dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y)),
    dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
    dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
    dhvm(caugi(A %-->% B %-->% C %-->% D, U %-->% B + D))
  )
  for (g in models) {
    sm <- screen_model(g)
    for (j in seq_along(sm$poly$Hrep)) {
      H <- sm$poly$Hrep[[j]]
      eq <- H[, 1] == 1
      Q <- rbind(H[eq, , drop = FALSE], cbind(0, 0, diag(ncol(H) - 2)))
      lp <- vapply(which(!eq), \(r) {
        out <- rcdd::lpcdd(rcdd::d2q(Q), rcdd::d2q(-H[r, -(1:2)]), minimize = FALSE)
        rcdd::q2d(out$optimal.value) > H[r, 2] + 1e-9
      }, logical(1))
      expect_equal(sm$nt[[j]][!eq], lp)
    }
  }

})
