test_that("check dimension of input", {
  # skip_on_cran()
  n = 25; X = matrix(floor(runif(n*3,0,100)), n, 3);
  Y= matrix(floor(runif(n-1,0,10)), n-1, 1)
  expect_error(robin(X, Y, E=NULL), "Length of Y does not match", ignore.case = TRUE)

  Y= matrix(floor(runif(n,0,10)), n, 1)
  expect_error(robin(X, Y, E=NULL, clin=1), "clin has a different number of rows.", ignore.case = TRUE)
  expect_error(robin(X, Y, E=NULL), "E factors must be provided.", ignore.case = TRUE)
  expect_error(robin(X, Y, E=1), "E has a different number of rows", ignore.case = TRUE)

})


test_that("check design matrix", {
  # skip_on_cran()
  n = 25; env = 2; size = env+1
  X = matrix(floor(runif(n*3,0,100)), n, 3); Y= matrix(floor(runif(n,0,10)), n, 1)
  E = matrix(floor(runif(n*env,0,20)), n, env)
  fit=robin(X, Y, E)

  expect_equal(dim(fit$design), c(n, 3*size))
  expect_equal(fit$design[,1], X[,1])
  expect_equal(fit$design[,2], X[,1]*E[,1])
  expect_equal(fit$design[,3], X[,1]*E[,2])
  expect_equal(fit$design[,size+1], X[,2])
  expect_equal(fit$design[,size+2], X[,2]*E[,1])
})


test_that("check parameters", {
  # skip_on_cran()
  n = 25; env = 2; size = env+1
  X = matrix(floor(runif(n*3,0,100)), n, 3); Y= matrix(floor(runif(n,0,10)), n, 1)
  E = matrix(floor(runif(n*env,0,20)), n, env)
  fit=robin(X, Y, E)

  expect_equal(fit$iterations, 10000)
  expect_equal(fit$burn.in, 5000)

})
