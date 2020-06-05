test_that("check_prediction_robust", {
  test = sample((1:nrow(X)), floor(nrow(X)/4))
  fit=roben(X[-test,], Y[-test,], E[-test,], clin[-test,], iterations = 5000)
  out = predict(fit, X[test,], E[test,], clin[test,], Y[test,])
  expect_equal(names(out$error), "PMAD")
  expect_length(out$y.pred, length(test))
})

test_that("check_prediction_robust", {
  test = sample((1:nrow(X)), floor(nrow(X)/4))
  fit=roben(X[-test,], Y[-test,], E[-test,], clin[-test,], iterations = 5000)
  expect_error(predict(fit, X[test,-1], E[test,], clin[test,], Y[test,]), "X.new")
  expect_error(predict(fit, X[test,], E[test,], clin[test,-1], Y[test,]), "clinical covariates")
})


test_that("check_prediction_nonrobust", {
  test = sample((1:nrow(X)), floor(nrow(X)/4))
  fit=roben(X[-test,], Y[-test,], E[-test,], clin[-test,], iterations = 5000, robust=FALSE)
  out = predict(fit, X[test,], E[test,], clin[test,], Y[test,])
  expect_equal(names(out$error), "PMSE")
  expect_length(out$y.pred, length(test))
})
