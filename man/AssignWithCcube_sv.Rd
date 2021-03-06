% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccube_sv.R
\name{AssignWithCcube_sv}
\alias{AssignWithCcube_sv}
\title{Computing responsibilities (assignment probabilities) with modified Ccube_sv VBEM-step}
\usage{
AssignWithCcube_sv(res, ssm, tol = 1e-08, maxiter = 100, epi = 0.001,
  verbose = F)
}
\arguments{
\item{res}{A reference Ccube results list}

\item{ssm}{A data frame of SVs to be assigned. Ideally, the data has been processed by CcubeSV model. So it should have ccube_mult1 and ccube_mult2 columns.}

\item{tol}{convergence threshold}

\item{maxiter}{maximum iterations, default number is 100.}

\item{epi}{sequencing error, default is 1e-3}

\item{verbose}{show VBEM progress}
}
\value{
a standard CcubeSV model with recomputed responsibilities and logResponsibilities
}
\description{
Computing responsibilities (assignment probabilities) with modified Ccube_sv VBEM-step
}
