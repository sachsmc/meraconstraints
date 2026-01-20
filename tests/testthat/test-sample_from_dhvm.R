test_that("Sampling is done correctly", {

  set.seed(120)
  ivdag <- caugi(Z %-->% X %-->% Y, U %-->% X + Y) |> dhvm()
  ivsam <- dhvm_sample_distribution(ivdag)
  expect_equal(sum(ivsam$prob.joint), 1)

  ## no latents

  dag2 <- caugi(Z %-->% X %-->% Y) |> dhvm()
  d2sam <- dhvm_sample_distribution(dag2)
  expect_equal(sum(d2sam$prob.joint), 1)
  ## iv constraints satisfied

  ivconst <- dhvm_meraconstraints(ivdag)
  ivsam$prob.test <- ivsam$prob.Y_X * ivsam$prob.X_Z
  ivsam <- merge(ivconst$functional_relations$problevels[[1]], ivsam, sort = FALSE)

  expect_true(hrep_test(ivconst$polytope$Hrep[[1]], ivsam$prob.test))

  ## iv constraints violated

  set.seed(118)
  cdag <- caugi(Z %-->% X %-->% Y,
                #Z %-->% Y,
                U %-->% Z +  X + Y) |> dhvm()
  csam <- dhvm_sample_distribution(cdag)
  csam$prob.test <- csam$prob.joint / csam$prob.Z
  csam <- merge(ivconst$functional_relations$problevels[[1]], csam, sort = FALSE)

  expect_false(hrep_test(ivconst$polytope$Hrep[[1]], csam$prob.test))

  ## more than 2 levels

  ivdag3 <- caugi(Z %-->% X %-->% Y, U %-->% X + Y) |> dhvm(nvals = c(Z = 3))
  ivsam3 <- dhvm_sample_distribution(ivdag3)

  expect_equal(nrow(ivsam3), 2 * 2 * 3)
  expect_equal(sum(ivsam3$prob.joint), 1)

  ivdag33 <- caugi(Z %-->% X %-->% Y, U %-->% X + Y) |> dhvm(nvals = c(Z = 3, X = 3))
  ivsam33 <- dhvm_sample_distribution(ivdag33)

  expect_equal(nrow(ivsam33), 2 * 3 * 3)
  expect_equal(sum(ivsam33$prob.joint), 1)

  ivdag333 <- caugi(Z %-->% X %-->% Y, U %-->% X + Y) |> dhvm(nvals = c(Z = 3, X = 3, Y = 3))
  ivsam333 <- dhvm_sample_distribution(ivdag333)

  expect_equal(nrow(ivsam333), 3 * 3 * 3)
  expect_equal(sum(ivsam333$prob.joint), 1)

  nolatent <- dhvm_sample_distribution(caugi(Z %-->% X %-->% Y, C %-->% Z + Y) |> dhvm(nvals = c(Z = 2, X = 2, Y = 2, C = 4)))

  expect_equal(nrow(nolatent), 2^3 * 4)
  expect_equal(sum(nolatent$prob.joint), 1)

  ## multiple latents

  ivd2 <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4 %-->% V5, U1 %-->% V2 + V3, U2 %-->% V4 + V5))
  ivsam2 <- dhvm_sample_distribution(ivd2, nsamp_latent = 100)
  expect_equal(sum(ivsam2$prob.joint), 1)

})
