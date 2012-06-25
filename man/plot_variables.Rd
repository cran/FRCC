\name{plot_variables}
\alias{plot_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variables in the Canonical Factor Loadings Space
}
\description{
 This function plots the variables used in the FRCCA as points in a two-dimensional 
 plane in which the axis are the canonical factor loadings selected by the user.
}
\usage{
plot_variables(res.mrcc, i, j, inner_circle_radius = 0.5, text_size = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.mrcc}{
  List containing a canonical structure provided by the function frcc.
}
  \item{i}{
 Canonical Factor Loadings which will be used as the horizontal axis.
}
  \item{j}{
  Canonical Factor Loadings  which will be used as the vertical axis.
}
  \item{inner_circle_radius}{
  Radius of the circle which is used to determine which variables are significant. Only the significant variables will be labled. 
}
  \item{text_size}{
 	Character expansion factor for the labels of the variables.
}
}

\value{
 This function just creates the variables plot. It does not return a value.
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
}
\author{
Raul Cruz-Cano
}

\examples{
# Example: Multivariate Normal Data
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
dev.new()
plot_variables(my_res,1,2,text_size=1.0)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
