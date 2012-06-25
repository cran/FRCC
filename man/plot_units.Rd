\name{plot_units}
\alias{plot_units}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the experimental units in the Canonical Variates Space
}
\description{
 This function plots the experimental units used in the FRCCA as points in a two-dimensional 
 plane in which the axis are the canonical variates selected by the user
}
\usage{
plot_units(X, Y, res.mrcc, i, text_size = 0.8, point_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
numeric matrix (n by p) which contains the observations on the X variables.
}
  \item{Y}{
numeric matrix (n by p) which contains the observations on the Y variables.
}
  \item{res.mrcc}{
List containing a canonical structure provided by the function frcc for the dataset X and Y.
}
  \item{i}{
Canonical Variate which will be used for the axes (X for horizontal and Y for vertical).
}
  \item{text_size}{
	Character expansion factor for the labels of the experimental units.
}
  \item{point_size}{
	Character expansion factor for the point representing the experimental units.
}
}

\value{

This function just creates the units plot. It does not return a value.
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
}
\author{
Raul Cruz-Cano
}

\examples{
#Example: NCI-60 micrRNA Data
data("Topoisomerase_II_Inhibitors")
data("microRNA")
my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
for( i in 1:dim(microRNA)[2])
{
 colnames(microRNA)[i]<-substr(colnames(microRNA)[i], 1, 2)
}#end for i
dev.new()
plot_units(t(microRNA),-1*t(Topoisomerase_II_Inhibitors),my_res,1,1,text_size=0.01)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
