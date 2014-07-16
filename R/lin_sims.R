#' This function take a UCsims object and runs a set of logistic models on the data
#' if coeffs == TRUE, it returns a list of model coefficients 
#' if coef == FALSE, it returns a list of glm model objects
bin_UC_models <- function(UCsim_obj, coefs = TRUE, ...){
    full_mod <- glm(Y ~ E + U1 + U2, family = binomial("logit"), data = UCsim_obj$data, ...)
    reduced_mod <- glm(Y ~ E + U2, family = binomial("logit"), data = UCsim_obj$data, ...)
    expose0_mod <- glm(Y ~ U1 + U2, family = binomial("logit"), data = UCsim_obj$data[ UCsim_obj$data$E == 0,], ...)
    expose1_mod <- glm(Y ~ U1 + U2, family = binomial("logit"), data = UCsim_obj$data[ UCsim_obj$data$E == 1,], ...)
    if(coefs) {
        lapply(list(full_mod = full_mod, reduced_mod = reduced_mod, expose0_mod = expose0_mod, expose1_mod = expose1_mod), coef)
    } else list(full_mod = full_mod, reduced_mod = reduced_mod, expose0_mod = expose0_mod, expose1_mod = expose1_mod)
}


#' For testing - all parameters set to 1
null_test_model <- function(UCsim_obj){
    list(full_mod = c(E = 1, U1 = 1, U2 = 1),
         reduced_mod = c(E = 1, U1 = 1, U2 = 1),
         expose0_mod = c(E = 1, U1 = 1, U2 = 1),
         expose1_mod = c(E = 1, U1 = 1, U2 = 1))
}

#' Performs the adjustment for unmeasured confounding according to Lin
#' 
#' @export
#' 
lin_adjust <- function(sims_obj, BC = "U1", MC = NULL, subgroups = NULL){
    mods <- bin_UC_models(sims_obj)
    meanU_e0 <- mean(sims_obj$data[[BC]][sims_obj$data$E == 0])
    meanU_e1 <- mean(sims_obj$data[[BC]][sims_obj$data$E == 1])
    lambda0 <- mods$expose0_mod[[BC]]
    lambda1 <- mods$expose1_mod[[BC]]
    orTau0 <- exp(lambda0)
    orTau1 <- exp(lambda1)
    mOReyPart <- exp(mods$reduced_mod[["E"]])
    mOReyFull <- exp(mods$full_mod[["E"]])
    bA <- (orTau1 * meanU_e1 + (1 - meanU_e1)) / (orTau0 * meanU_e0 + (1 - meanU_e0))
    cA <- (meanU_e1 - meanU_e0) * (lambda0 + lambda1) / 2
    bOReyPartAdj = mOReyPart / bA      # OR adjusted assuming binary UC
    cOReyPartAdj <- mOReyPart / exp(cA) # OR adjusted for continuous UC)
    c(sims_obj$parameters, meanE = mean(sims_obj$data$E), meanY = mean(sims_obj$data$Y), meanU = mean(sims_obj$data[[BC]]), 
      meanU_e0 = meanU_e0, meanU_e1 = meanU_e1,
      mOReyFull = mOReyFull, 
      mOReyPart = mOReyPart,
      orTau0 = orTau0, orTau1 = orTau1, 
      bA = bA, bOReyPartAdj = bOReyPartAdj,     
      cA = cA, cOReyPartAdj = cOReyPartAdj,
      bOReydiff = bOReyPartAdj - mOReyFull,
      cOReydiff = cOReyPartAdj - mOReyFull) 
}


