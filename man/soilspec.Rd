\name{soilspec}
\alias{soilspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Soil Specification Data
}
\description{
Contains the Soil Specification Data.
}
\usage{
data(soilspec)
}
\details{
The original purpose of the experiment was to determine the relationships
between several soil characteristics of the limestone grassland in Wales and the abundance
of certain plant species. These variables were measured in a random sample of 10 x 10 square
meters in the community of Anglesey, North Wales. The dataset comprises data from 45 samples
on 8 species of plants (H. pubescens, P. bertolonii, T. pretense, P. sanguisorba, R. squarrosus, H.
pilosella, B. media and T. drucei) and 3 soil characteristics (d=depth, P= extractable phosphate
and K=exchangeable potassium) and their interactions (d x P, d x K and P x K). Previous work
shows that these soil characteristics are influential in determining how much each of the existing
plant species can flourish. This set of plants was selected because they have a diverse response to
variation of the soil variables.
}
\value{
A matrix with the information corresponding to the 8 types of soils (columns 2-9) and the soil characteristics and their interactions (columns 10-16) for the 45 soil samples as described in [Cruz-Cano et al., 2012]. The first columns keeps track of the number of the site of origin. The orginal source of the data is [Gittins, 2005]. 
It was first used as an R dataset in [De'ath and Walsh, 2001]. 

}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
Gittins, R.;  Canonical Analysis - A Review with Applications in Ecology. Biomathematics 12, 1985.
De'ath, G; Walsh, CJ (2001). The pcurve Package Principal Curve Analysis. Documentation for  R: A language and environment for statistical computing. Foundation for Statistical Computing,Vienna, Austria. URL http://www.R-project.org, 2001.

}
\author{
Raul Cruz-Cano
}

\examples{
#Example #2 Soil Specification Data
data(soilspec)
list_of_units_to_be_used<-sample(1:nrow(soilspec),14)#We will only 14 soil samples
X<- soilspec[list_of_units_to_be_used,2:9]
Y<- soilspec[list_of_units_to_be_used,10:15]
colnames(X)<-c("H. pubescens", "P. bertolonii", "T. pretense", 
"P. sanguisorba", "R. squarrosus", "H. pilosella", "B. media","T. drucei")
colnames(Y)<- c("d","P","K","d x P", "d x K","P x K")
my_res<-frcc(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ soil }
