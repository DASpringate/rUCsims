###########################
### RUN THIS CODE FIRST ###
###########################

source("R/UC_sims.R")
source("R/R_and_R_sims.R")
source("R/lin_sims.R")



library(ggplot2)
library(scales)
library(reshape2)
library(parallel)
if(Sys.info()['sysname'] == "Windows"){
    my_cores <- 1 
} else {
    my_cores <- 12    
} 


#############################
### Simulate stata tables ###
#############################

# generate the simulations:
# This code replicates the original stata tables and produces the split line graphs (lin_chart.pdf)

UBs <- c(0.1, 0.333, 0.5, 2, 3, 5, 10, 20)
run1 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 100000, probY_0 = 0.1, probE_0 = 0.1, probU = 0.1)
run2 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 100000, probY_0 = 0.5, probE_0 = 0.5, probU = 0.5)
run3 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 100000, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9)

run4 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 100000, rho = 0.5, probY_0 = 0.1, probE_0 = 0.1, probU = 0.1)
run5 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 100000, rho = 0.5, probY_0 = 0.5, probE_0 = 0.5, probU = 0.5)
run6 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
                 obs = 10000, rho = 0.5, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9)

run7 <- run_sims(UBs, reps = 100, cores = my_cores, method = "no_moderator", adjust_fn = lin_adjust, BC = "U1",  
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

ggsave("figure/lin_chart.pdf")



##############################################################
### Test for effect of changes in UB on overall prevalence ###
##############################################################


#' This code produces the graphs that test for the effect of changes in UB on overall prevelance
#' Also produces the prev_test.pdf charts


test_prev <- function(UBs, ...){
    d <- as.data.frame(t(sapply(UBs, function(U){
        dat <- simulate_UC(orUB = U, ...)
        colMeans(dat$data)
    })))
    d$UB <- UBs
    d
}

UBs <- c(0.01, 0.05, 0.1, 0.333, 0.5, 1, 2, 3, 5, 10, 20)

prev_list <- list("Py0=0.1, rho=0" = test_prev(UBs, obs = 100000, rho = 0.0, probY_0 = 0.1, probE_0 = 0.1, probU = 0.1),
                  "Py0=0.5, rho=0" = test_prev(UBs, obs = 100000, rho = 0.0, probY_0 = 0.5, probE_0 = 0.5, probU = 0.5),
                  "Py0=0.9, rho=0" = test_prev(UBs, obs = 100000, rho = 0.0, probY_0 = 0.9, probE_0 = 0.9, probU = 0.9))

prev_dat <- do.call(`rbind`, lapply(names(prev_list), function(L){
    d <- prev_list[[L]]
    d$run <- L
    d
}))

molten <- melt(prev_dat[,c(1,2,3,5,6)], id.vars=c("UB", "run"))
p <- ggplot(molten, aes(x = UB, y = value, colour = variable))
p + geom_line() + 
    facet_wrap(~ run) +
    scale_x_continuous(trans=log_trans(), breaks = c(0.01, 0.1, 0.5, 2, 5, 10, 20)) + 
    scale_y_continuous("Prevalence")
ggsave("figure/prev_test.pdf")


#####################################################
### Runs with stabilised initial parameter values ###
#####################################################

# Set your parameter values here:

# Lin method
run8 <- stable_run(Ubs = c(0.1, 0.5, 1, 2, 10), 
                   obs_stable = 100000,  cores = 12, obs_run = 10000, reps = 10,
                   tol_Y = 0.005, tol_E = 0.005, pY_target = 0.5, pE_target = 0.5, 
                   probE_0 = 0.5, probY_0 = 0.5, orEY_0 = 1, probU = 0.3, 
                   modUB = 1, orME = 1, orMY_0 = 1, rho = 0.5, 
                   method = "continuous_confounder", adjust_fn = lin_adjust,
                   agg_subset = c("orUE", "probE_0", "probY_0", "orEY_0", 
                                  "probU", "mOReyFull", "mOReyPart", 
                                  "bOReyPartAdj", "cOReyPartAdj"))
# R and R method
run9 <- stable_run(Ubs = c(0.1, 0.5, 1, 2, 10), 
                   obs_stable = 100000,  cores = 12, obs_run = 10000, reps = 10,
                   tol_Y = 0.005, tol_E = 0.005, pY_target = 0.5, pE_target = 0.5, 
                   probE_0 = 0.5, probY_0 = 0.5, orEY_0 = 1, probU = 0.3, 
                   modUB = 1, orME = 1, orMY_0 = 1, rho = 0.5, subgroups = 5,
                   method = "continuous_confounder", adjust_fn = r_and_r_adjust)

# combine data from the two:
test_dat <- run8$aggregates
test_dat$r_and_r <- run9$aggregates$OReyPartAdj

#plot
molten = melt(test_dat[, c("orUE", "mOReyFull", "mOReyPart", "bOReyPartAdj", "cOReyPartAdj", "r_and_r")], id.vars = c("orUE"))
p <- ggplot(molten, aes(x = orUE, y = value, colour = variable))
p + geom_line() + geom_point() + 
    scale_x_continuous(trans=log_trans(), breaks = c(0.1, 0.5, 2, 5, 10, 20)) + 
    scale_y_continuous("Odds ratio")
ggsave("figure/lin_vs_rr_performance_with_rho0.5.pdf")









