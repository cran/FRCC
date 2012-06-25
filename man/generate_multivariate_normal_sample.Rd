\name{generate_multivariate_normal_sample}
\alias{generate_multivariate_normal_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
It generates a sample from a multinormal distribution function 
}
\description{
It generates a sample from a multinormal distribution function with the cross-covariance matrix described in [Cruz-Cano et al. 2012].
}
\usage{
generate_multivariate_normal_sample(p, q, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Number of desired variables in the dataset X.
}
  \item{q}{
Number of desired variables in the dataset Y.
}
  \item{n}{
 sample size desired.
}
}

\value{
A list of n sample units with the values for the variables of the datsets X and Y.
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
}
\author{
Raul Cruz-Cano
}

\examples{
p<-10
q<-10
n<-50
res<-generate_multivariate_normal_sample(p,q,n)
X<-res$X
Y<-res$Y
rownames(X)<-c(1:n)
colnames(X)<-c(1:p)
colnames(Y)<- c(1:q)
my_res<-frcc(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
