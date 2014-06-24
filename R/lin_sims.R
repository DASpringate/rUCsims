#' This function generates the simulated data given the parameters
#' Return is an object of class 'UCsims' which is a list containing a dataframe of the simulated data and a vector of parameter values.
#' For generality, the confounders are labelled U1 and U2.
#' In this case, U1 is binary and U2 is continuous
#' @export
lin_sim <- function(obs = 1000, orUB = 0.1, probE_0 = 0.2, probY_0 = 0.7, orEY_0 = 1, probU = 0.4, 
                    modUB = 1, orME = 2, orMY_0 = 2, rho = 0.0, orUBe = NULL, orUBy_0 = NULL){
    oddsE_0 <- probE_0 / (1 - probE_0)
    oddsY_0 <- probY_0 / (1 - probY_0)
    probY_1 <- orEY_0 * probY_0 / (1 - probY_0 + orEY_0 * probY_0)
    if(is.null(orUBe)) orUBe <- orUB
    if(is.null(orUBy_0)) orUBy_0 <- orUB
    UP <- runif(obs)
    Ulatent <- qnorm(1 - UP)
    # binary unmeasured confounder:
    U_binary <- as.integer(UP < probU)
    # Continuous measured confounder correlated with Ulatent at value rho
    M_continuous <- 0 + rho * Ulatent + sqrt(1 - rho * rho) * rnorm(obs)
    # Generate data for the exposure as a function of the set of confounders
    exposure <-  as.integer(runif(obs) < plogis(log(oddsE_0) + log(orUBe) * U_binary))
    # Now generate data for the outcome as a function of the exposure and all confounders
    y <- as.integer(runif(obs) < plogis(log(oddsY_0) + log(orEY_0) * exposure + log(orUBy_0) * U_binary))
    structure(list(data = data.frame(E = exposure, Y = y, U1 = U_binary, U2 = M_continuous),
                   parameters = c(sapply(names(formals()), function(x) eval(parse(text=x))), 
                                  oddsE_0 = oddsE_0, oddsY_0 = oddsY_0, probY_1 = probY_1)), class = "UCsims")
}

#' This function take a UCsims object and runs a set of logistic models on the data
#' if coeffs == TRUE, it returns a list of model coefficients 
#' if coef == FALSE, it returns a list of glm model objects
#' @export
bin_UC_models <- function(UCsim_obj, coefs = TRUE, ...){
    full_mod <- glm(Y ~ E + U1 + U2, family = binomial("logit"), data = UCsim_obj$data, ...)
    reduced_mod <- glm(Y ~ E + U2, family = binomial("logit"), data = UCsim_obj$data, ...)
    expose0_mod <- glm(Y ~ U1 + U2, family = binomial("logit"), data = UCsim_obj$data[ UCsim_obj$data$E == 0,], ...)
    expose1_mod <- glm(Y ~ U1 + U2, family = binomial("logit"), data = UCsim_obj$data[ UCsim_obj$data$E == 1,], ...)
    if(coefs) {
        lapply(list(full_mod = full_mod, reduced_mod = reduced_mod, expose0_mod = expose0_mod, expose1_mod = expose1_mod), coef)
    } else list(full_mod = full_mod, reduced_mod = reduced_mod, expose0_mod = expose0_mod, expose1_mod = expose1_mod)
}


#' wrapper function for the sims_fn to stop any perfectly collinear datasets being generated
#' @export
catch_collinear <- function(orUB, sims_fn, UC, ...){
    dat <- sims_fn(orUB = orUB, ...)
    if(mean(dat$data$Y[dat$data$E == 0 & dat$data[[UC]] == 1]) %in% c(0, 1) ||
           mean(dat$data$Y[dat$data$E == 1 & dat$data[[UC]] == 1]) %in% c(0, 1)) {
        message("Regenerating data due to collinearity")
        catch_collinear(orUB, sims_fn, UC, ...)
    } else dat
}



#' For testing - all parameters set to 1
null_test_model <- function(UCsim_obj){
    list(full_mod = c(E = 1, U1 = 1, U2 = 1),
         reduced_mod = c(E = 1, U1 = 1, U2 = 1),
         expose0_mod = c(E = 1, U1 = 1, U2 = 1),
         expose1_mod = c(E = 1, U1 = 1, U2 = 1))
}


#' wrapper function for lin_sim and lin_models
#' Runs both functions over a range of orUB levels
#' returns a dataframe of parameter values from the simulation and the models
#' calculates the adjustments for binary and continuous UCs
#' @export
run_UC_levels <- function(orUB_levels, sims_fn = lin_sim, model_fn = bin_UC_models, UC = "U1", ...){
    as.data.frame(t(sapply(orUB_levels, function(lev){
        dat <- catch_collinear(orUB = lev, sims_fn = sims_fn, UC = UC, ...) #sims_fn(orUB = lev, ...)
        mods <- model_fn(dat)
        meanU_e0 <- mean(dat$data[[UC]][dat$data$E == 0])
        meanU_e1 <- mean(dat$data[[UC]][dat$data$E == 1])
        lambda0 <- mods$expose0_mod[[UC]]
        lambda1 <- mods$expose1_mod[[UC]]
        orTau0 <- exp(lambda0)
        orTau1 <- exp(lambda1)
        mOReyPart <- exp(mods$reduced_mod[["E"]])
        mOReyFull <- exp(mods$full_mod[["E"]])
        bA <- (orTau1 * meanU_e1 + (1 - meanU_e1)) / (orTau0 * meanU_e0 + (1 - meanU_e0))
        cA <- (meanU_e1 - meanU_e0) * (lambda0 + lambda1) / 2
        bOReyPartAdj = mOReyPart / bA      # OR adjusted assuming binary UC
        cOReyPartAdj <- mOReyPart / exp(cA) # OR adjusted for continuous UC)
        c(dat$parameters, meanE = mean(dat$data$E), meanY = mean(dat$data$Y), meanU = mean(dat$data[[UC]]), 
          meanU_e0 = meanU_e0, meanU_e1 = meanU_e1,
          mOReyFull = mOReyFull, 
          mOReyPart = mOReyPart,
          orTau0 = orTau0, orTau1 = orTau1, 
          bA = bA, bOReyPartAdj = bOReyPartAdj,     
          cA = cA, cOReyPartAdj = cOReyPartAdj,
          bOReydiff = bOReyPartAdj - mOReyFull,
          cOReydiff = cOReyPartAdj - mOReyFull) 
    })))
}

#' This is the replication function.
#' Runs run_UC_levels over reps replicates and binds to a single dataframe
#' Runs simulations in parallel if cores > 1 (needs to have multicore package installed)
#' Set cores to 1 if running on a windows machine
#' @export
rep_sims <- function(orUB_levels, reps = 10, cores = 1, ...){
    if(cores > 1){ # set to 1 for a windows machine.
        library(multicore)
        do.call(`rbind`, mclapply(1:reps, function(x){
            d <- run_UC_levels(orUB_levels = orUB_levels, ...)
            d$reps <- x
            cat(".")
            d
        }, mc.cores = cores))
    } else {
        do.call(`rbind`, lapply(1:reps, function(x){
            d <- run_UC_levels(orUB_levels = orUB_levels, ...)
            d$reps <- x
            cat(".")
            d
        }))
    }
}

#' Helper function to automate running different conditions
#' @export
run_sims <- function(Ubs, aggregator = median, rawdata = TRUE, ...){
    dat <- rep_sims(orUB_levels = Ubs, ...)
    dat_reduced <- dat[, c("orUB", "probE_0", "probY_0", "orEY_0", "probU", "mOReyFull", "mOReyPart", "bOReyPartAdj", "cOReyPartAdj")]
    if(rawdata){
        list(data = dat, aggregates = aggregate(dat_reduced, by = list(orUB = dat_reduced$orUB), FUN = aggregator))
    } else aggregate(dat_reduced, by = list(orUB = dat_reduced$orUB), FUN = aggregator)
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
overall_prev <- function(sims_fn = lin_sim, iter_fn = iter_fn_1, pY_target, pE_target, probE_0 = 0.2, 
                         probY_0 = 0.7, tol_Y = 0.0005, tol_E = 0.0005, verbose = TRUE,  ...){
    sim_dat <- sims_fn(..., probE_0 = probE_0, probY_0 = probY_0)
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
        sim_dat <- sims_fn(..., probE_0 = new_pE0, probY_0 = new_pY0)
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
#' @param sims_fn function defining the simulation model
#' @param iter_fn function defining how new stabilising parameters are generated
#' @param cores number of processor cores to use (passed to mclapply)
#' @param \dots arguments to be passed to sims_fn
#' @return list containing data and aggregated data
stable_run <- function(Ubs, pY_target, pE_target, obs_stable, obs_run, 
                       sims_fn = lin_sim, iter_fn = iter_fn_1, model_fn = bin_UC_models, cores = 1, reps = 10, ...){
    out <- do.call(`rbind`, lapply(Ubs, function(U){
        message("orU = ", U)
        message("Iterating to find stable parameters...")
        stable_params <- lapply(overall_prev(orUB = U, sims_fn = sims_fn, iter_fn = iter_fn, 
                                      pY_target = pY_target, pE_target = pE_target, obs = obs_stable, verbose = FALSE, ...)$parameters,
                                function(x) x)
        stable_params$obs <- obs_run
        stable_params <- stable_params[names(stable_params) %in% names(formals(sims_fn))]
        stable_params[["orUB"]] <- NULL
        message("Running simulations with stablised parameters (", reps, " reps)...")
        do.call(`run_sims`, c(stable_params, cores = cores, Ubs = U, sims_fn = sims_fn, model_fn = model_fn, reps = reps))$data
    }))
    dat_reduced <- out[, c("orUB", "probE_0", "probY_0", "orEY_0", "probU", "mOReyFull", "mOReyPart", "bOReyPartAdj", "cOReyPartAdj")]
    list(data = out, aggregates = aggregate(dat_reduced, by = list(orUB = dat_reduced$orUB), FUN = median))
}



