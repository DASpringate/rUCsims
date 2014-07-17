#' This function generates the simulated exposure and response data, based on a binary and an 
#' optional continuous confounder.
#' 
#' The exposure is generated as a function of the confounders and the response is then generated 
#' as a function of the exposure and the confounders 
#' 
#' The binary confounder is labeled "U1" and the continuous "U2". 
#' Confounders can be correlated, using the `rho` argument.  
#' Linear moderation can also be introduced via the `modUB` argument between the continuous 
#' confounder and the exposure/outcome.
#' 
#' The method argument specifies how the exposure (E) and outcome (Y) variables are generated.  
#' The "continuous_confounder" option means that E and Y are based on both the binary and 
#' continuous confounder. The "no_continuous_confounder" option means that E and Y are a function
#' of just the binary confounder.  If this is selected, the unmeasured confounder U2 is still generated
#' and may be correlated with U1, although it will be uncorrelated with both E and Y.
#' @export
#'
#' @param obs the number of observations in the simulated data (number of rows)
#' @param orUE the odds ratio OR(ue)
#' @param probE_0 probability that E=1 when (all confounders)= 0
#' @param probY_0  p(y=1/exposure=0 & confounders=0)
#' @param orEY_0 odds ratio OR(ey)
#' @param probU  p(u=1) prevalence of the binary unmeasured confounder
#' @param modUB moderating effect of the UC on OR(ey)
#' @param orME odds ratio OR(me)
#' @param orMY_0 odds ratio OR(my/e=0)
#' @param rho correlation coefficient between U1 and U2
#' @param orUY_0
#' @param method 
#' @return  object of class 'UCsims' which is a list containing a dataframe of the simulated data, a vector of parameter values and the method.
#' 
simulate_UC <- function(obs, orUE, probE_0, probY_0, orEY_0, probU, 
                    modUB, orME, orMY_0, rho, orUY_0 = NULL, 
                    method = c("continuous_confounder", "no_continuous_confounder")){
    oddsE_0 <- probE_0 / (1 - probE_0)
    oddsY_0 <- probY_0 / (1 - probY_0)
    probY_1 <- orEY_0 * probY_0 / (1 - probY_0 + orEY_0 * probY_0)
    if(is.null(orUY_0)) orUY_0 <- orUE
    UP <- runif(obs)
    Ulatent <- qnorm(1 - UP)
    # binary unmeasured confounder:
    U_binary <- as.integer(UP < probU)
    # Continuous measured confounder correlated with Ulatent at value rho
    M_continuous <- 0 + rho * Ulatent + sqrt(1 - rho * rho) * rnorm(obs)
    
    method <- match.arg(method)
    parameters = sapply(names(formals()), function(x) eval(parse(text=x)))
    parm_names <- names(parameters)
    parameters = as.numeric(parameters[names(parameters) != "method"])
    names(parameters) <- parm_names[parm_names != "method"]
    
    if(method == "continuous_confounder"){
        # Generate data for the exposure as a function of the set of confounders
        exposure <-  as.integer(runif(obs) < plogis(log(oddsE_0) + log(orUE) * U_binary))
        # Now generate data for the outcome as a function of the exposure and all confounders
        y <- as.integer(runif(obs) < plogis(log(oddsY_0) + log(orEY_0) * exposure + log(orUY_0) * U_binary + log(modUB)*(U_binary)*(exposure)))
                
    } else if(method == "no_continuous_confounder"){
        # Generate data for the exposure as a function of the set of confounders
        # Model: binary UC with measured confounder and linear moderation
        exposure <-  as.integer(runif(obs) < plogis(log(oddsE_0) + log(orUE) * U_binary + log(orME) * M_continuous))
        # Now generate data for the outcome as a function of the exposure and all confounders
        y <- as.integer(runif(obs) < plogis(log(oddsY_0) + log(orEY_0) * exposure + log(orUY_0) * U_binary + 
                                                log(modUB) * (U_binary) * exposure + log(orMY_0) * M_continuous))
    } else stop("You need to specify a method for determining exposure and outcome!")
    structure(list(data = data.frame(E = exposure, Y = y, U1 = U_binary, U2 = M_continuous),
                          parameters = c(parameters, 
                                         oddsE_0 = oddsE_0, oddsY_0 = oddsY_0, probY_1 = probY_1),
                          method = method), class = "UCsims")
    
}


#' wrapper function for simulate_UC to stop any perfectly collinear datasets being generated
catch_collinear <- function(orUE, BC, MC, subgroups,  ...){
    dat <- simulate_UC(orUE = orUE, ...)
    if(mean(dat$data$Y[dat$data$E == 0 & dat$data[[BC]] == 1]) %in% c(0, 1) ||
           mean(dat$data$Y[dat$data$E == 1 & dat$data[[BC]] == 1]) %in% c(0, 1)) {
        message("Regenerating data due to collinearity")
        catch_collinear(orUE, BC, MC, subgroups, ...)
    } else dat
}

#' wrapper function for adjust functions
#' Runs both functions over a range of orUE levels
#' returns a dataframe of parameter values from the simulation and the models
#' calculates the adjustments for binary and continuous UCs
#' @export
run_UC_levels <- function(orUE_levels, method = "no_moderator", adjust_fn = lin_adjust, 
                          BC = "U1", MC = "U2", subgroups = 5, ...){
    as.data.frame(t(sapply(orUE_levels, function(lev){
        dat <- catch_collinear(orUE = lev, BC = BC, method = method, ...) #sims_fn(orUE = lev, ...)
        adjust_fn(dat, subgroups = subgroups, BC = BC, MC = MC) 
    })))
}


#' This is the replication function.
#' Runs run_UC_levels over reps replicates and binds to a single dataframe
#' Runs simulations in parallel if cores > 1 )
#' Set cores to 1 if running on a windows machine
#' @export
rep_sims <- function(orUE_levels, reps = 10, cores = 1, ...){
    do.call(`rbind`, mclapply(1:reps, function(x){
        d <- run_UC_levels(orUE_levels = orUE_levels, ...)
        d$reps <- x
        cat(".")
        d
    }, mc.cores = cores))
}


#' Helper function to automate running different conditions
#' @export
run_sims <- function(Ubs, aggregator = median, rawdata = TRUE, agg_subset = NULL, ...){
    dat <- rep_sims(orUE_levels = Ubs, ...)
    if(!is.null(agg_subset)){
        dat_reduced <- dat[, agg_subset]         
    }else dat_reduced <- dat
    
    if(rawdata){
        list(data = dat, aggregates = aggregate(dat_reduced, by = list(orUE = dat_reduced$orUE), FUN = aggregator))
    } else aggregate(dat_reduced, by = list(orUE = dat_reduced$orUE), FUN = aggregator)
}




#' Helper function for overall_prev
#' This determines how the new values are generated
#' @export
iter_fn_1 <- function(prob_in, target, par_mean){
    prob_out <- prob_in * 2 * target / (target + par_mean)
    if(prob_out > 1) prob_out <- 1
    if(prob_out < 0) prob_out <- 0
    prob_out
}

#' This function uses an iterative algorithm to obtain overall prevalences of Y and E.
#' @export
overall_prev <- function(method, iter_fn = iter_fn_1, pY_target, pE_target, probE_0 = 0.2, 
                         probY_0 = 0.7, tol_Y = 0.0005, tol_E = 0.0005, verbose = TRUE,  ...){
    sim_dat <- simulate_UC(..., probE_0 = probE_0, probY_0 = probY_0, method = method)
    diff_Y <- 1
    diff_E <- 1
    iters <- 0
    while(diff_Y > tol_Y || diff_E > tol_E){
        old_pE0 <- sim_dat$parameters[["probE_0"]]
        old_pY0 <- sim_dat$parameters[["probY_0"]]
        new_pE0 <- iter_fn(old_pE0, pE_target, mean(sim_dat$data$E))
        new_pY0 <- iter_fn(old_pY0, pY_target, mean(sim_dat$data$Y))
        diff_E <- abs(mean(sim_dat$data$E) - pE_target)
        diff_Y <- abs(mean(sim_dat$data$Y) - pY_target)
        sim_dat <- simulate_UC(..., probE_0 = new_pE0, probY_0 = new_pY0, method = method)
        iters <- iters + 1
        if(verbose){
            message("iteration #", iters, 
                    "\n\tprobE_0 = ", old_pE0,  ", newPE = ", new_pE0, ", targetE = ", pE_target, ", meanE = ", mean(sim_dat$data$E), 
                    "\n\tprobY_0 = ", old_pY0,  ", newPY = ", new_pY0, ", targetY = ", pY_target, ", meanY = ", mean(sim_dat$data$Y))
        } else {
            cat(".")
        }
        
    }
    message("tolerance met in ", iters, " iterations.")
    sim_dat
}


#' This function generates stable parameters across a range of ORs and then produces data for replicates at these parameter levels.
#' 
#' @export
#' 
#' @param Ubs vector of odds ratios
#' @param pY_target numeric target for the total population prevalence of Y
#' @param pE_target numeric target for the total population prevalence of E
#' @param obs_stable number of observations for the stabilising run
#' @param obs_run number of observations for the main run
#' @param method the simulation method to be used.  See `simulate_UC`
#' @param iter_fn function defining how new stabilising parameters are generated
#' @param adjust_fn the method to be used for adjustment for unmeasured confounders
#' @param cores number of processor cores to use (passed to mclapply)
#' @param \dots arguments to be passed to sims_fn
#' @return list containing data and aggregated data
stable_run <- function(Ubs, pY_target, pE_target, obs_stable, obs_run, method,
                       iter_fn = iter_fn_1, adjust_fn, cores = 1, reps = 10, agg_subset = NULL, subgroups = NULL, ...){
    out <- do.call(`rbind`, lapply(Ubs, function(U){
        message("orU = ", U)
        message("Iterating to find stable parameters...")
        stable_params <- lapply(overall_prev(orUE = U, iter_fn = iter_fn, method = method,
                                             pY_target = pY_target, pE_target = pE_target, obs = obs_stable, verbose = FALSE, ...)$parameters,
                                function(x) x)
        stable_params$obs <- obs_run
        stable_params <- stable_params[names(stable_params) %in% names(formals(simulate_UC))]
        stable_params[["orUE"]] <- NULL
        message("Running simulations with stablised parameters (", reps, " reps)...")
        do.call(`run_sims`, c(stable_params, cores = cores, Ubs = U, method = method, 
                              adjust_fn = adjust_fn, reps = reps, subgroups = subgroups))$data
    }))
    if(!is.null(agg_subset)){
        dat_reduced <- out[, agg_subset]
    } else dat_reduced <- out
    list(data = out, aggregates = aggregate(dat_reduced, by = list(orUE = dat_reduced$orUE), FUN = median))
}


