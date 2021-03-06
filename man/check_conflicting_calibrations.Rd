% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_calibrations_bladj.R
\name{check_conflicting_calibrations}
\alias{check_conflicting_calibrations}
\title{function to check for conflicting calibrations
if calibrations are younger or older relative to descendants and ancestors, respectively}
\usage{
check_conflicting_calibrations(phy, calibration_distribution)
}
\arguments{
\item{phy}{A phylo object}

\item{calibration_distribution}{is a list of node ages distributions, named with the node number from phy}
}
\description{
function to check for conflicting calibrations
if calibrations are younger or older relative to descendants and ancestors, respectively
}
\details{
it removes them if needed, but bladj works as long as it has an age for the root
}
