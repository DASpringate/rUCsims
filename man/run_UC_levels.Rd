\name{run_UC_levels}
\alias{run_UC_levels}
\title{wrapper function for lin_sim and lin_models
Runs both functions over a range of orUB levels
returns a dataframe of parameter values from the simulation and the models
calculates the adjustments for binary and continuous UCs}
\usage{
run_UC_levels(orUB_levels, sims_fn = lin_sim, model_fn = bin_UC_models,
  UC = "U1", ...)
}
\description{
wrapper function for lin_sim and lin_models Runs both
functions over a range of orUB levels returns a dataframe
of parameter values from the simulation and the models
calculates the adjustments for binary and continuous UCs
}

