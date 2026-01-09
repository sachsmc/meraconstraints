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

})
