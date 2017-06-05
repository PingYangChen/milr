# rfda v0.3.0

* Remove dependency: assertthat, purrr
* Performance tuning with rewriting some functions with `RcppParallel`.
* Fix the NULL return in `predict.milr` and `predict.softmax`.
* Add some tests for milr.
* Add vignette based on the published R Journal.

# rfda v0.2.0

* New options for milr: numLambda, tolerance.
* Add option for output the labels of bags/instances into `fitted` and `predict` functions.
* Use the coefficients of ridge regression acquired by `glmnet` as initial values for p > n problems.

# rfda v0.1.0

* First release on CRAN.
