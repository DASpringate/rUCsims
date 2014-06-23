\name{rep_sims}
\alias{rep_sims}
\title{This is the replication function.
Runs run_UC_levels over reps replicates and binds to a single dataframe
Runs simulations in parallel if cores > 1 (needs to have multicore package installed)
Set cores to 1 if running on a windows machine}
\usage{
rep_sims(orUB_levels, reps = 10, cores = 1, ...)
}
\description{
This is the replication function. Runs run_UC_levels over
reps replicates and binds to a single dataframe Runs
simulations in parallel if cores > 1 (needs to have
multicore package installed) Set cores to 1 if running on a
windows machine
}

