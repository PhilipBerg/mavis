#' Spiked-in data set of peptides from total proteomics
#'
#' A dataset containing quantification of reversibly oxidized cysteines using MaxQuant \insertCite{tyanova2016maxquant}{mavis}.
#' True positives peptides spiked-in from the Universal Protein Standard at nine different concentrations
#' and true negatives from \emph{Saccharomyces cerevisiae} with the same
#' concentration in all samples. For details see
#' \insertCite{ramus2016benchmarking}{mavis} and if you use this dataset
#' please cite the same paper.
#' @format A data frame with 2235 rows and 7 variables:
#' \describe{
#'   \item{identifier}{id column for features, true positives contains UPS and
#'   true negatives contains YEAST}
#'   \item{\[a-i\]_* }{Biological replicates with true positives spiked-in at 0.05,
#'   0.125, 0.250, 0.5, 2.5, 5.0, 12.5, 25.0 50.0 fmol/\eqn{\mu}g yeast lysate, with increasing spike-in by alphabetic order}
#' }
#' @source \url{https://www.sciencedirect.com/science/article/pii/S187439191530186X}
#' @references
#' \insertAllCited{}
"ramus_max"
