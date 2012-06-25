\name{rearrange.frcc}
\alias{rearrange.frcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Rearranges the canonical structure according to the canonical correlations
}
\description{
 By using the minimum risk estimators of the correlation matrices instead of the 
 sample correlation matrices the FRCC algoeithm might disrupt the order of the
 canonical correlations and hence of the canonical structure. This is unacceptable for the algorithm used  to calculate the p-values 
 which requires the canonical correltions to be ordered in a descending order. This function rearranges the canonical structure according 
 to the canonical correlations from largest to smallest. 
}
\usage{
rearrange.frcc(res.frcc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.frcc}{
  
  List containing a canonical structure produced by the function frcc.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{res.frcc }{List containing the sorted canonical structure.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review.
}
\author{
Raul Cruz-Cano
}

\examples{
##   This is an internal function. No need for examples.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
