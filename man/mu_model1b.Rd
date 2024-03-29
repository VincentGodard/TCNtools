% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balco2017production.R
\name{mu_model1b}
\alias{mu_model1b}
\title{Muons production profile according to model 1B of Balco (2017)}
\usage{
mu_model1b(z, h, Rc, Sphi, cst)
}
\arguments{
\item{z}{depth below the surface z (g/cm2)}

\item{h}{site atmospheric pressure h (hPa)}

\item{Rc}{cutoff rigidity (GV)}

\item{Sphi}{solar modulation parameter}

\item{cst}{list of nuclide specific parameters (see below)}
}
\description{
Compute muons production profile as a function of depth below the surface z (g/cm2)
and site atmospheric pressure h (hPa). According to model 1B from Balco (2017)
}
\details{
Adapted from G. Balco matlab code

Balco, G. (\strong{2017}). Production rate calculations for cosmic-ray-muon-produced 10Be and 26Al benchmarked against geological calibration data.
\emph{Quaternary Geochronology}, 39, 150–173.
https://doi.org/10.1016/j.quageo.2017.02.001

cst is a list containing nuclide-specific constants, as follows:
cst$Natoms -  atom number density of the target atom (atoms/g)
cst$k_neg - summary cross-section for negative muon capture (atoms/muon)
cst$sigma190 - 190 GeV x-section for fast muon production (cm2)
}
