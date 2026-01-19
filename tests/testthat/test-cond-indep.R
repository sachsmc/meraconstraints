test_that("conditional independence", {

  dag <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, U %-->% V2 + V4))
  cd1 <- dag_list_independences(dag)
  expect_true(
    setequal(cd1[[1]]$indep, c("V1", "V3")) &
     setequal(cd1[[1]]$given, "V2")
  )


  dag2 <- dhvm(caugi(A %-->% C, B %-->% C, C %-->% D))
  cd2 <- dag_list_independences(dag2)

  expect_true(setequal(cd2[[1]]$indep, c("A", "B")) & is.null(cd2[[1]]$given) &
    setequal(cd2[[2]]$indep, c("A", "D")) & setequal(cd2[[2]]$given, "C") &
    setequal(cd2[[3]]$indep, c("B", "D")) & setequal(cd2[[3]]$given, "C") )


  dag3 <- dhvm(caugi(V1 %-->% V2, V2 %-->% V3, V3 %-->% V4, U1 %-->% V1+V3, U2 %-->% V2+V4))
  expect_true(length(dag_list_independences(dag3)) == 0)

})
