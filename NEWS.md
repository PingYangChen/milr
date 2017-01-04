
# rfda v0.2.0

* Remove dependency: purrr 
* Performance tuning with rewriting some functions with C++.
* New options for milr: numLambda, tolerance.
* Add option for output the labels of bags/instances into `fitted` and `predict` functions.
* Use the coefficients of ridge regression acquired by `glmnet` as initial values for p > n problems.

# rfda v0.1.0

* First release on CRAN.
