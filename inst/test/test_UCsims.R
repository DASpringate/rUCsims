

simulate_UC(obs = 10000)
simulate_UC(obs = 10000, method = "moderator")
catch_collinear(0.5, BC = "U1", MC = "U2", subgroups = NULL, obs = 10000)
catch_collinear(0.5, BC = "U1", MC = "U2", subgroups = NULL, method = "mod", obs = 10000)

run_UC_levels(c(0.1, 1, 2), obs = 10000)
run_UC_levels(c(0.1, 1, 2), obs = 100000, method = "mod", adjust_fn = r_and_r_adjust)

run_UC_levels(c(0.1, 1, 2), obs = 100000, method = "no_mod", adjust_fn = lin_adjust, 
              probE_0 = 0.1, probY_0 = 0.2, orEY_0 = 1, probU = 0.3, 
              modUB = 1, orME = 1, orMY_0 = 1, rho = 0.0, orUBe = NULL, orUBy_0 = 5)


rep_sims(c(0.1, 1, 2), obs = 10000)
run_sims(c(0.1, 1, 2), agg_subset = c("orUB", "probE_0", "probY_0", "orEY_0", 
                                      "probU", "mOReyFull", "mOReyPart", 
                                      "bOReyPartAdj", "cOReyPartAdj"), obs = 10000)
run_sims(c(0.1, 1, 2, 10), obs = 100000, probE_0 = 0.1, probY_0 = 0.2, orEY_0 = 1, probU = 0.3, 
         modUB = 1, orME = 1, orMY_0 = 1, rho = 0.0, orUBe = NULL, orUBy_0 = 5, 
         method = "mod", adjust_fn = r_and_r_adjust, cores = 12)

overall_prev(method = "moderator", iter_fn = iter_fn_1, pY_target = 0.5, pE_target = 0.5, probE_0 = 0.2, 
                         probY_0 = 0.7, tol_Y = 0.0005, tol_E = 0.0005, verbose = TRUE)

run8 <- stable_run(Ubs = c(0.1, 0.5, 1, 2, 10), 
                   obs_stable = 100000,  cores = 12, obs_run = 10000, reps = 10,
                   tol_Y = 0.005, tol_E = 0.005, pY_target = 0.5, pE_target = 0.5, 
                   probE_0 = 0.5, probY_0 = 0.5, orEY_0 = 1, probU = 0.3, method = "no_moderator", 
                   modUB = 10, orME = 1, orMY_0 = 1, rho = 0, agg_subset = c("orUB", "probE_0", "probY_0", "orEY_0", 
                                                                             "probU", "mOReyFull", "mOReyPart", 
                                                                             "bOReyPartAdj", "cOReyPartAdj"))

run_sims(c(0.1, 0.5, 1, 2, 10), obs = 10000, probE_0 = 0.5, probY_0 = 0.5, orEY_0 = 1, probU = 0.3, 
         modUB = 10, orME = 1, orMY_0 = 1, rho = 0.0, orUBe = NULL, orUBy_0 = 5, 
         method = "mod", adjust_fn = lin_adjust, cores = 12)





# tests
a <- simulate_UC(obs = 100000, orUE = 5, probE_0 = 0.1, probY_0 = 0.2, orEY_0 = 1, probU = 0.3, 
                 modUB = 1, orME = 1, orMY_0 = 1, rho = 0.0, orUY_0 = 5, method = "no")

aa <- rep_sims(orUE_levels= c(0.1, 0.5, 1, 2, 10), 
               cores = 1, obs = 10000, reps = 10,
               probE_0 = 0.5, probY_0 = 0.5, orEY_0 = 1, probU = 0.3, 
               modUB = 1, orME = 1, orMY_0 = 1, rho = 0, 
               method = "no_continuous_confounder", adjust_fn = r_and_r_adjust)


exp(glm(E ~ U1 + U2, family = binomial("logit"), data = a$data)$coefficients)
a$parameters
exp(glm(Y ~ E + U1 + U2 + E:U1, family = binomial("logit"), data = a$data)$coefficients)
exp(glm(Y ~ E + U2, family = binomial("logit"), data = a$data[ a$data$U1 == 0,])$coefficients)
exp(glm(Y ~ E + U2, family = binomial("logit"), data = a$data[ a$data$U1 == 1,])$coefficients)
sims_obj <- rr_sim(obs = 100000)

sg <- subgroup_analysis(a)


