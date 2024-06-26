% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mitch.R
\name{mitch_calc}
\alias{mitch_calc}
\title{mitch_calc}
\usage{
mitch_calc(
  x,
  genesets,
  minsetsize = 10,
  cores = 1,
  resrows = 50,
  priority = NULL
)
}
\arguments{
\item{x}{a multicolumn numerical table with each column containing 
differential expression scores for a contrast.
Rownames must match genesets.}

\item{genesets}{lists of genes imported by the gmt_imprt function or
similar.}

\item{minsetsize}{the minimum number of genes required in a set for it to be
included in the statistical analysis.
Default is 10.}

\item{cores}{the number of parallel threads for computation. Defaults to 1.}

\item{resrows}{an integer value representing the number of top genesets for
which a detailed report is to be 
generated. Default is 50.}

\item{priority}{the prioritisation metric used to selecting top gene sets. 
Valid options are 'significance', 
'effect' and 'SD'.}
}
\value{
mitch res object with the following parts:
$input_profile: the supplied input differential profile
$input_genesets: the supplied input gene sets
$ranked_profile: the differential profile after ranking
$enrichment_result: the table of MANOVA/ANOVA enrichment results for each
gene set
$analysis_metrics:  several metrics that are important to the interpretation
of the results
$detailed_sets: a list of dataframes containing ranks of members of
prioritised gene sets.
}
\description{
This function performs multivariate gene set enrichment analysis.
}
\examples{
# Example using mitch to calculate multivariate enrichments and
# prioritise based on effect size 
data(myImportedData,genesetsExample)
resExample<-mitch_calc(myImportedData,genesetsExample,priority='effect',
minsetsize=5,cores=2)
}
\keyword{calc}
\keyword{calculate}
\keyword{manova}
\keyword{mitch}
