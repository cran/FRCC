\name{Topoisomerase_II_Inhibitors}
\alias{Topoisomerase_II_Inhibitors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NCI-60 Topoisomerase II Inhibitor Data.
}
\description{
Load a matrix with Topoisomerase II Inhibitor Drugs datatset. 
}
\usage{
data(Topoisomerase_II_Inhibitors)
}
%- maybe also 'usage' for other objects documented here.
\details{
The NCI-60 is a set of cell cultures grown under controlled conditions by the National Cancer Institute. The NCI-60 cell lines include experimental units from the breast (8), central nervous system (6),
colorectal (7), lung (9), prostate (2), ovarian (6) and renal (8) cancers. It also includes leukemia (6) and melanoma (8) cell lines. 
}
\value{
A matrix with the growth inhibitory responses of the 15 Topoisomerase II Inhibitor Drugs for the 60 cell lines in the NCI-60 dataset as described in [Cruz-Cano et al., 2012].  This is a subset of the A118 drugs dataset orginal described  in [DTP, 2009].
}
\references{
Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis, under review, 2012.
DTP (2009). DTP Human Tumor Cell Line Screen.. Standard mechanism. Available from: http://dtp.nci.nih.gov/
}
\author{
Raul Cruz-Cano
}

\examples{
#Example #3 NCI-60 micrRNA Data
data("Topoisomerase_II_Inhibitors")
data("microRNA")
my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Topoisomerase_II_Inhibitor}
\keyword{ NCI-60 }% __ONLY ONE__ keyword per line
