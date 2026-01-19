test_that("Known constraints are recovered", {
  ivhrep <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                      0, 1, 1, 0, 1, 0, 1, 0, 0, -1, -1, 0, 1, 0, 1, 0, -1, 0, 0, 0,
                      -1, 0, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, 0, 0,
                      1, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0,
                      1, 0, -1, 1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
                      -1, 0, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 1), dim = c(14L, 10L), representation = "H")


  ivdag <- graph_from_literal(Z -+ X, X -+ Y, U -+ X:Y) |> initialize_graph()
  ivrels <- dag_get_linear_relations(ivdag)
  hrepiv <- cf_relations_to_hrep(ivrels)
  expect_equal(hrepiv$Hrep[[1]], ivhrep)

  ivd <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y))
  ivrels2 <- dhvm_linear_relations(ivd)
  hrepiv2 <- cf_relations_to_hrep(ivrels2)
  expect_equal(hrepiv2$Hrep[[1]], ivhrep)



  prop7 <- dhvm(caugi(U1 %-->% V1 + V3 + V4,
                      U2 %-->% V2 + V4 + V5,
                      V3 + V4 + V5 %-->% V1 + V2))

  prop7 <- dhvm(caugi(U1 %-->% V1 + V2,
                      U2 %-->% V4 + V5,
                      V4 %-->% V5,
                      V3 + V4 + V5 %-->% V1 + V2, V3 %-->% V4))

  p7rels <- dhvm_linear_relations(prop7)
  hrepp7 <- cf_relations_to_hrep(p7rels)

  hrep_to_text(hrepp7$Hrep[[1]], p7rels$probform[[1]])[prob_simplex_extremes(p7rels, hrepp7$Hrep)[[1]]]

})
