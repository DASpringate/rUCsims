rUCsims 
========

A package for simulating confounded datasets and accounting for unmeasured confounders
-----------------------------------------------------

This package is still in early development

Install via github and devtools:

```R
library(devtools)
install_github("DASpringate/rUCsims")
library(rUCsims)
```

This package implements the following methods for dealing with unmeasured confounders:

* Lin et al [reference]()
* Extension to the Lin method assuming continuous confounders

Not yet implemented:

* Rosenbaum and rubin [reference]()
* Bayesian methods [reference]()



examples (See lin_analysis.R):

### Run the Lin method across a set of levels of unmeasured confounders

```R
UBs <- c(0.1, 0.333, 0.5, 2, 3, 5, 10, 20)
run1 <- run_sims(UBs, reps = 20, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 100000, probY_0 = 0.1, probE_0 = 0.1, probU = 0.1)
run2 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 100000, probY_0 = 0.5, probE_0 = 0.5, probU = 0.5)
run3 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 100000, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9)

run4 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 100000, rho = 0.5, probY_0 = 0.1, probE_0 = 0.1, probU = 0.1)
run5 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 100000, rho = 0.5, probY_0 = 0.5, probE_0 = 0.5, probU = 0.5)
run6 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 10000, rho = 0.5, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9)

run7 <- run_sims(UBs, reps = 100, cores = 20, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1",  
                 obs = 10000, rho = 0.5, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9)


## plot the results:
bin_dat <- do.call(`rbind`, lapply(c("run1", "run2", "run3", "run4", "run5", "run6"), function(x){
    d <- eval(parse(text = x))$aggregates
    d$run <- x
    d
}))

molten <- melt(bin_dat[, c("orUB", "mOReyFull", "mOReyPart", "bOReyPartAdj", "cOReyPartAdj", "run")], id.vars = c("orUB", "run"))
p <- ggplot(molten, aes(x = orUB, y = value, colour = variable))
p + geom_line() + 
    facet_wrap(~ run) +
    scale_x_continuous(trans=log_trans(), breaks = c(0.1, 0.5, 2, 5, 10, 20)) + 
    scale_y_continuous("Odds ratio")
```


### Run the lin method over datasets with stabilised prevalences:

```R
run8 <- stable_run(Ubs = c(0.1, 0.5, 1, 2, 10), pY_target = 0.6, pE_target = 0.4, probE_0 = 0.5, probY_0 = 0.5, 
                          obs_stable = 1000000, obs_run = 10000, orEY_0 = 10,
                          probU = 0.6, modUB = 4, orME = 7, orMY_0 = 2, rho = 0, tol_Y = 0.0005, tol_E = 0.0005, cores = 1)

# Plot the results
molten = melt(run8$aggregates[, c("orUB", "mOReyFull", "mOReyPart", "bOReyPartAdj", "cOReyPartAdj")], id.vars = c("orUB"))
p <- ggplot(molten, aes(x = orUB, y = value, colour = variable))
p + geom_line() + 
    scale_x_continuous(trans=log_trans(), breaks = c(0.1, 0.5, 2, 5, 10, 20)) + 
    scale_y_continuous("Odds ratio")
```


