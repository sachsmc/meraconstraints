test_that("conditional independence", {

  dag <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, U %-->% V2 + V4))
  cd1 <- dhvm_list_independences(dag)
  expect_true(
    setequal(cd1[[1]]$indep, c("V1", "V3")) &
     setequal(cd1[[1]]$given, "V2")
  )


  dag2 <- dhvm(caugi(A %-->% C, B %-->% C, C %-->% D))
  cd2 <- dhvm_list_independences(dag2)

  ## A _||_ B, and D _||_ (A, B) | C in chained form
  expect_true(setequal(cd2[[1]]$indep, c("A", "B")) & is.null(cd2[[1]]$given))
  sets2 <- attr(cd2, "sets")
  expect_true(any(vapply(sets2, \(s) s$indep[[1]] == "D" && setequal(s$indep[[2]], c("A", "B")) &&
                           setequal(s$given, "C"), logical(1))))
  expect_true(all(vapply(cd2, \(ii) d_separated(dag2$caugi_graph, X = ii$indep[1], Y = ii$indep[2],
                                                Z = ii$given), logical(1))))


  dag3 <- dhvm(caugi(V1 %-->% V2, V2 %-->% V3, V3 %-->% V4, U1 %-->% V1+V3, U2 %-->% V2+V4))
  expect_true(length(dhvm_list_independences(dag3)) == 0)

})
