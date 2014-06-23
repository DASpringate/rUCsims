\name{stable_run}
\alias{stable_run}
\title{This function generates stable parameters across a range of ORs and then produces data for replicates at these parameter levels.}
\usage{
stable_run(Ubs, pY_target, pE_target, obs_stable, obs_run, sims_fn = lin_sim,
  iter_fn = iter_fn_1, cores = 1, ...)
}
\arguments{
  \item{Ubs}{vector of odds ratios}

  \item{pY_target}{numeric target for the total population
  prevalence of Y}

  \item{pE_target}{numeric target for the total population
  prevalence of E}

  \item{obs_stable}{number of observations for the
  stabilising run}

  \item{obs_run}{number of observations for the main run}

  \item{sims_fn}{function defining the simulation model}

  \item{iter_fn}{function defining how new stabilising
  parameters are generated}

  \item{cores}{number of processor cores to use (passed to
  mclapply)}

  \item{\dots}{arguments to be passed to sims_fn}
}
\value{
list containing data and aggregated data
}
\description{
This function generates stable parameters across a range of
ORs and then produces data for replicates at these
parameter levels.
}

