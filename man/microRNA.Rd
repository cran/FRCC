\name{microRNAs}
\alias{microRNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NCI-60 microRNA data
}
\description{
Contains the expression level of 365 microRNA data in the NCI-60 cell lines. 
}
\usage{
data(microRNA)
}
\details{
The NCI-60 is a set of cell cultures grown under controlled conditions by the National Cancer Institute. The NCI-60 cell lines include experimental units from the breast (8), central nervous system (6),
colorectal (7), lung (9), prostate (2), ovarian (6) and renal (8) cancers. It also includes leukemia (6) and melanoma (8) cell lines. MicroRNAs are a type of RNA molecules found in eukaryotic cells. Each microRNA is a short
RNA sequence (around 22 nucleotides) which is involved in the regulation of multiple target genes. A large number of published papers deal with the problem of finding the microRNA expression signature of different cancers with the goal of designing early detection methods and providing
therapeutic targets.
}
\value{
A matrix with the expression level of 365 microRNAs for the 60 cell lines in the NCI-60 dataset as described in [Cruz-Cano et al., 2012]. The orginal source of the dataset is [DTP, 2009]  
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
\keyword{microRNA}
\keyword{ NCI-60 }% __ONLY ONE__ keyword per line
