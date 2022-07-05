#' Spiked-in data set of peptides
#'
#' A dataset containing quantification of peptides using Progenesis.
#' True positives peptides spiked-in from the Universal Proteomics Standard Set
#' 1 (UPS1) at three different concentrations
#' and true negatives from \emph{Chlamydomonas reinhardtii} with the same
#' concentration in all samples. For details see
#' \insertCite{berg2019evaluation;textual}{mavis} and if you use this dataset
#' please cite the same paper.
#' @format A data frame with 10599 rows and 13 variables:
#' \describe{
#'   \item{identifier}{id column for features, true positives contains UPS and
#'   true negatives contains Cre}
#'   \item{fmol25_* }{Technical replicates with true positives spiked-in from 25
#'   fmol UPS1 peptides}
#'   \item{fmol50_* }{Technical replicates with true positives spiked-in from 50
#'   fmol UPS1 peptides}
#'   \item{fmol100_* }{Technical replicates with true positives spiked-in from
#'   100 fmol UPS1 peptides}
#' }
#' @source \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2619-6}
#' @references
#' \insertAllCited{}
"ups_prog"
