## Examples from "Deriving Complete Constraints in Hidden Variable Models"

constraint_counts <- function(cn) {
  lapply(seq_along(cn$polytope$Hrep), \(j) {
    H <- cn$polytope$Hrep[[j]]
    fl <- cn$nontrivial[[j]]
    c(ineq = sum(H[, 1] == 0), eq = sum(H[, 1] == 1),
      fl_ineq = sum(fl & H[, 1] == 0), fl_eq = sum(fl & H[, 1] == 1))
  })
}

test_that("conditioning sets are minimal (Proposition 3)", {

  rg <- dhvm(caugi(Z %-->% X, X %-->% Y, W %-->% Y, H %-->% X + Y), latent.pattern = "^H")
  cs <- dhvm_district_params(rg)$plist[[1]]
  expect_setequal(cs$X, "Z")
  expect_setequal(cs$Y, c("X", "Z", "W"))

})

test_that("districts and c-degrees (Figures 1 and 7, and Section 5)", {

  g1 <- dhvm(caugi(A %-->% B + E, B %-->% C, C %-->% D, D %-->% F, E %-->% C, G %-->% F,
                   U1 %-->% B + D, U2 %-->% A + E + C, U3 %-->% F + D))
  expect_equal(unname(dhvm_c_degree(g1)), c(1, 2, 0))

  g7 <- dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                   U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))
  expect_equal(unname(dhvm_c_degree(g7)), c(1, 2))
  expect_equal(unname(dhvm_c_degree(dhvm_merge_latents(g7))), c(1, 1))

  tri <- dhvm(caugi(U1 %-->% V1 + V3, U2 %-->% V1 + V2, U3 %-->% V2 + V3))
  expect_equal(unname(dhvm_c_degree(tri)), 3)
  expect_equal(unname(dhvm_c_degree(dhvm_merge_latents(tri))), 1)

})

test_that("constraint counts match Table 1", {

  iv <- meraconstraints(dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)))
  expect_equal(constraint_counts(iv)[[1]], c(ineq = 12, eq = 2, fl_ineq = 4, fl_eq = 0))
  expect_equal(nrow(unique(rcdd::q2d(iv$polytope$Vrep[[1]]))), 12)

  iv3 <- meraconstraints(dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3)))
  expect_equal(constraint_counts(iv3)[[1]], c(ineq = 48, eq = 3, fl_ineq = 36, fl_eq = 0))
  expect_equal(nrow(unique(rcdd::q2d(iv3$polytope$Vrep[[1]]))), 28)

  right <- meraconstraints(dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)))
  expect_equal(constraint_counts(right),
               list(c(ineq = 24, eq = 6, fl_ineq = 8, fl_eq = 0),
                    c(ineq = 4, eq = 2, fl_ineq = 0, fl_eq = 0)))

  xab <- meraconstraints(dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)))
  expect_equal(constraint_counts(xab)[[1]], c(ineq = 24, eq = 8, fl_ineq = 8, fl_eq = 0))

  f7 <- dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                   U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))
  f7m <- meraconstraints(dhvm_merge_latents(f7))
  expect_equal(constraint_counts(f7m),
               list(c(ineq = 16, eq = 5, fl_ineq = 0, fl_eq = 0),
                    c(ineq = 1128, eq = 12, fl_ineq = 1096, fl_eq = 6)))

})

test_that("sequential IV gives two copies of the IV constraints (Figure 3)", {

  sq <- meraconstraints(dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4 %-->% V5,
                                   U1 %-->% V2 + V3, U2 %-->% V4 + V5)))
  expect_equal(names(sq$polytope$Hrep), c("V2_V3", "V4_V5"))
  expect_equal(constraint_counts(sq),
               rep(list(c(ineq = 12, eq = 2, fl_ineq = 4, fl_eq = 0)), 2))

  ## V4 _||_ (V1, V2) | V3 and V5 _||_ (V1, V2) | V3, V4, which give (V1, V2) _||_ (V4, V5) | V3
  sets <- attr(sq$independences, "sets")
  expect_length(sets, 2)
  expect_true(all(vapply(sets, \(s) setequal(s$indep[[2]], c("V1", "V2")), logical(1))))
  expect_true(setequal(sets[[1]]$given, "V3") && setequal(sets[[2]]$given, c("V3", "V4")))

})

test_that("front door has no nontrivial constraints (Section 3.3)", {

  fd <- meraconstraints(dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y)))
  expect_equal(names(fd$polytope$Hrep), c("A_Y", "M"))
  expect_equal(constraint_counts(fd)[[1]], c(ineq = 8, eq = 3, fl_ineq = 0, fl_eq = 0))
  expect_equal(sum(fd$nontrivial$M), 0)
  expect_identical(fd$trivial_if_false, fd$nontrivial)

  ## the screen that ignores the equality constraints flags 4 of them
  old <- prob_simplex_screen(fd$functional_relations, fd$polytope$Hrep)
  txt <- hrep_to_text(fd$polytope$Hrep[[1]], fd$functional_relations$probstarform[[1]])
  expect_setequal(txt[old[[1]]], c(
    "-P*(A = 0,Y = 0|M = 0) - P*(A = 0,Y = 1|M = 0) + P*(A = 0,Y = 0|M = 1) <= 0",
    "P*(A = 0,Y = 0|M = 0) + P*(A = 0,Y = 1|M = 0) + P*(A = 1,Y = 0|M = 1) <= 1",
    "P*(A = 0,Y = 0|M = 0) + P*(A = 0,Y = 1|M = 0) - P*(A = 0,Y = 0|M = 1) - P*(A = 0,Y = 1|M = 1)  =  0",
    "-P*(A = 0,Y = 0|M = 0) - P*(A = 0,Y = 1|M = 0) - P*(A = 1,Y = 0|M = 1) - P*(A = 1,Y = 1|M = 1)  =  -1"
  ))

})

test_that("closed-form simplex screening agrees with enumeration", {

  models <- list(
    dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
    dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3)),
    dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y)),
    dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V1 %-->% V4, U2 %-->% V2 + V4)),
    dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
    dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B))
  )
  for (m in models) {
    rels <- dhvm_linear_relations(m)
    hreps <- dhvm_relations_to_polytope(rels)$Hrep
    expect_equal(prob_simplex_screen(rels, hreps), prob_simplex_extremes(rels, hreps))
  }

})

test_that("nested Markov equivalent graphs are distinguished (Figure 5)", {

  left  <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V1 %-->% V4, U2 %-->% V2 + V4))
  right <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4))
  cn_left  <- meraconstraints(left)
  cn_right <- meraconstraints(right)
  expect_equal(constraint_counts(cn_left)[[1]], c(ineq = 16, eq = 6, fl_ineq = 0, fl_eq = 0))

  ## P*(V2, V4 | V1, V3) = P(V2 | V1) P(V4 | V1, V2, V3)
  pstar <- function(probs, cn) {
    probs$prob.test <- probs$prob.joint / (probs$prob.V3_V2 * probs$prob.V1)
    merge(cn$functional_relations$problevels[["V2_V4"]], probs, sort = FALSE)$prob.test
  }

  set.seed(2024)
  viol_left <- viol_right <- FALSE
  for (i in 1:50) {
    probs <- dhvm_sample_distribution(left, parmsamp = \(k) rnorm(k, sd = 5), nsamp_latent = 200)
    viol_left <- viol_left | !hrep_test(cn_left$polytope$Hrep[["V2_V4"]], pstar(probs, cn_left), which = TRUE)
    viol_right <- viol_right | !hrep_test(cn_right$polytope$Hrep[["V2_V4"]], pstar(probs, cn_right), which = TRUE)
  }
  expect_equal(sum(viol_left), 0)
  expect_gt(sum(viol_right), 0)
  ## only inequality constraints of the right graph can be violated
  expect_true(all(cn_right$polytope$Hrep[["V2_V4"]][viol_right, 1] == 0))

})

test_that("Pearl and Bonet inequalities are not sufficient with a ternary instrument", {

  cn <- meraconstraints(dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3)))
  lev <- cn$functional_relations$problevels[[1]][, c("X", "Y", "Z")]

  ## q_yx(z) = P(Y = y, X = x | Z = z)
  q <- rbind(c(19, 19, 1, 1), c(1, 19, 19, 1), c(19, 1, 1, 19)) / 40
  colnames(q) <- c("00", "01", "10", "11")
  p <- q[cbind(lev$Z + 1, match(paste0(lev$Y, lev$X), colnames(q)))]

  ## Pearl's instrumental inequality holds
  for (x in 0:1) {
    expect_lte(sum(sapply(0:1, \(y) max(q[, paste0(y, x)]))), 1)
  }

  ok <- hrep_test(cn$polytope$Hrep[[1]], p, which = TRUE)
  expect_false(all(ok))
  expect_equal(sum(!ok), 2)

  ## mixing with the uniform distribution moves it back inside the model
  expect_true(hrep_test(cn$polytope$Hrep[[1]], 0.5 * p + 0.5 * rep(1/4, length(p))))

})

test_that("triangle graph has no nontrivial constraints after merging (Section 5)", {

  tri <- dhvm(caugi(U1 %-->% V1 + V3, U2 %-->% V1 + V2, U3 %-->% V2 + V3))
  cn <- meraconstraints(dhvm_merge_latents(tri))
  expect_length(cn$independences, 0)
  expect_equal(sum(unlist(cn$trivial_if_false)), 0)

})
