Status
------

     License       |  Linux/osx Build  |   Windows Build   |

-------------------|-------------------|-------------------| [![GitHub license](https://img.shields.io/badge/lincense-MIT-blue.svg)](http://badges.mit-license.org) | [![Build status](https://travis-ci.org/ChingChuan-Chen/milr.svg?branch=master)](https://travis-ci.org/ChingChuan-Chen/milr/branches) | [![Build status](https://ci.appveyor.com/api/projects/status/2yms6ao3mf69fdht/branch/master?svg=true)](https://ci.appveyor.com/project/ChingChuan-Chen/milr/branch/master) |

milr
====

The multiple-instance logistic regression with lasso penalty.

Installation
------------

You can install:

-   install from CRAN

    ``` r
    install.packages("milr")
    devtools::install_github("PingYangChen/milr")
    ```

-   the latest development version from github with

    ``` r
    install.packages("devtools")
    devtools::install_github("PingYangChen/milr")
    ```

If you encounter a bug, please file a reproducible example on [github](https://github.com/PingYangChen/milr/issues).

examples
--------

    ```R
    set.seed(100)
    beta <- runif(5, -5, 5)
    trainData <- DGP(70, 3, beta)
    testData <- DGP(30, 3, beta)
    # default (not use LASSO)
    milr_result <- milr(trainData$Z, trainData$X, trainData$ID)
    coef(milr_result)      # coefficients
    fitted(milr_result)                    # fitted bag labels
    fitted(milr_result, type = "instance") # fitted instance labels
    summary(milr_result)   # summary milr
    predict(milr_result, testData$X, testData$ID)                    # predicted bag labels
    predict(milr_result, testData$X, testData$ID, type = "instance") # predicted instance labels
    ```
