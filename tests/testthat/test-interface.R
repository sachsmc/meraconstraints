test_that("Construction of dhvm works", {

  expect_no_error(test <- dhvm(caugi(Z %-->% X %-->% Y, U %-->% X + Y), latent.variables = "U"))

})
