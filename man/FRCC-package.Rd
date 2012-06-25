\name{FRCC-package}
\alias{FRCC-package}
\alias{FRCC}
\docType{package}
\title{
 Fast Regularized Canonical Correlation Analysis
}
\description{
This package implements the Fast Regularized Canonical Correlation Analysis described in [Cruz-Cano et al., 2012]. 
The main idea of the algorithm is using the minimum risk estimators of the correlation matrices described in [Schafer and Strimmer, 2008] during the calculation of the Canonical correlation Structure. 
It can be considered an extesion of the work for two set of variables (blocks) mentioned in [Tenenhaus and Tenenhaus, 2011]
}
\details{
\tabular{ll}{
Package: \tab FRCC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-13\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The function frcc provides the canonical structure for two sets of variables X and Y. The rest of the functions help to visualize and interpret the values of the canonical structure.
}
\author{
Raul Cruz-Cano

Maintainer: Raul Cruz-Cano <raulcruz@umd.edu>
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.

Schafer, J; Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation
and implications for functional genomics. Statistical Applications in Genetics and Molecular
Biology 4:14, Article 32.

Tenenhaus, A.; Tenenhaus, M. (2011). Regularized Generalized Canonical Correlation Analysis.
Psychometrika 76:2, DOI: 10.1007/S11336-011-9206-8.
}
\keyword{ Fast Regularized Canonical Correlation Analysis }
\keyword{ Regularized Canonical Correlation Analysis }
\examples{
# Examples of the functions included in this package are listed
# in the help file of each individual function.
}
