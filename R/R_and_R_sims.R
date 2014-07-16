

#' performs the R and R sensitivity on the subgroups of the stratified measured confounder
#' 
#' Note that R&R are wrong in saying that alpha and gamma are for E=0; they are actually for E=1
#' To get the log odds and log odds-ratio for E=0 alpha and gamma need to be inverted
#' (because they are in logs, this is done by changing the signs)
#' @export
#' 
r_and_r_adjust <- function(sims_obj, subgroups = 5, BC = "U1", MC = "U2"){
    if(subgroups > 1){
        sims_obj$data$MC_strata <- cut(sims_obj$data[[MC]], breaks = subgroups, labels = FALSE)
    } else sims_obj$data$MC_strata <- 1
    # sensitivity parameters
    pi_ <- 1 - sims_obj$parameters[["probU"]]
    alpha <- -log(sims_obj$parameters[["orUE"]]) # negative to correct for an error in the paper
    delta0 <- log(1 / sims_obj$parameters[["orUY_0"]])
    delta1 <- log(1 / (x=sims_obj$parameters[["orUY_0"]] * sims_obj$parameters[["modUB"]]))
    
    parms <- data.frame(t(sapply(1:subgroups, function(subgroup){
        dat <- sims_obj$data[ sims_obj$data$MC_strata == subgroup,]
        
        c(nsub = nrow(dat), meanY <- mean(dat$Y), meanE = mean(dat$E), 
          meanE_u0 = mean(dat$E[dat[[BC]] == 0]),  meanE_u1 = mean(dat$E[dat[[BC]] == 1]),
          meanU = mean(dat[[BC]]), meanU_e0 = mean(dat[[BC]][dat$E == 0]), meanU_e1 = mean(dat[[BC]][dat$E == 1]),
          meanY = mean(dat$Y), meanY_e0 = mean(dat$Y[dat$E == 0]), meanY_e1 = mean(dat$Y[dat$E == 1]),
          meanY_u0 = mean(dat$Y[dat[[BC]] == 0]), meanY_u1 = mean(dat$Y[dat[[BC]] == 1]),
          meanY_e0u0 = mean(dat$Y[dat$E == 0 & dat[[BC]] == 0]), meanY_e1u0 = mean(dat$Y[dat$E == 1 & dat[[BC]] == 0]),
          meanY_e0u1 = mean(dat$Y[dat$E == 0 & dat[[BC]] == 1]), meanY_e1u1 = mean(dat$Y[dat$E == 1 & dat[[BC]] == 1]))
    })))
    parms$mean_not_E <- 1 - parms$meanE_u0
    a_ <- exp(alpha) * parms$mean_not_E
    b_ <- (parms$mean_not_E - pi_) * exp(alpha) + parms$mean_not_E + pi_ - 1
    c_ <- parms$mean_not_E - 1
    parms$gamma <- -log((-b_ + sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_)) # log odds of E=1 when U=0
    parms$odds_not_E_0 <- parms$mean_not_E / (1 - parms$mean_not_E)
    wt0 <- pi_ / (pi_ + ((1 - pi_) * (1 + exp(parms$gamma))/(1 + exp(parms$gamma + alpha)))) # estimated adjusted prob that U=0 when E=0?
    wt1 <- pi_ / (pi_ + ((1 - pi_) * exp(alpha) * (1 + exp(parms$gamma)) / (1 + exp(parms$gamma + alpha)))) # estimated adjusted prob that u=0 when E=1?
    parms$mean_not_Y_e0 <- 1 - parms$meanY_e0
    aa0 <- exp(delta0) * parms$mean_not_Y_e0
    bb0 <- (parms$mean_not_Y_e0 - wt0) * exp(delta0) + parms$mean_not_Y_e0 + wt0 - 1
    cc0 <- parms$mean_not_Y_e0 - 1
    parms$x0 <- (-bb0 + sqrt(bb0 * bb0 - 4 * aa0 * cc0)) / (2 * aa0) # odds that Y=1 when U=0
    parms$mean_not_Y_e1 <- 1 - parms$meanY_e1
    aa1 <- exp(delta1) * parms$mean_not_Y_e1
    bb1 <- (parms$mean_not_Y_e1 - wt1) * exp(delta1) + parms$mean_not_Y_e1 + wt1 - 1
    cc1 <- parms$mean_not_Y_e1 - 1
    parms$x1 <- (-bb1 + sqrt(bb1 * bb1 - 4 * aa1 * cc1)) / (2 * aa1) # odds that Y=1 when U=1
    parms$beta0 <- log(parms$x0) # log odds that Y=0 when U=0
    parms$beta1 <- log(parms$x1) # log odds that Y=0 when U=1
    parms$tau0 <- (1 - pi_) * (exp(parms$beta0 + delta0) / (1 + exp(parms$beta0 + delta0))) + pi_ * exp(parms$beta0) / 
        (1 + exp(parms$beta0)) #estimated adjusted %Y=1 at e=0
    parms$tau1 <- (1 - pi_) * (exp(parms$beta1 + delta1) / (1 + exp(parms$beta1 + delta1))) + pi_ * exp(parms$beta1) / 
        (1 + exp(parms$beta1))  #estimated adjusted %y=1 at e=1
    parms$adjOReyRR <- parms$tau1 * (1 - parms$tau0) / (parms$tau0 * (1 - parms$tau1))
    parms$psub <- parms$nsub / sum(parms$nsub)
    
    parms$totalProbY_e0 <- parms$psub * parms$tau0 # compute contribution to overall adjusted P(Y=1/E=0)
    parms$totalProbY_e1 <- parms$psub * parms$tau1
    
    parm_totals <- colMeans(parms[, c("totalProbY_e0", "totalProbY_e1")])
    
    totalAdjOREY <- parm_totals[["totalProbY_e1"]] * (1 - parm_totals[["totalProbY_e0"]]) / 
        (parm_totals[["totalProbY_e0"]] * (1 - parm_totals[["totalProbY_e1"]])) #compute overall adjusted OR(EY)
    full_mod <- glm(as.formula(paste("Y ~ E +", BC, "+", MC)), 
                    family = binomial("logit"), data = sims_obj$data)$coefficients
    reduced_mod <- glm(as.formula(paste("Y ~ E +", MC)), 
                       family = binomial("logit"), data = sims_obj$data)$coefficients
    
    c(sims_obj$parameters, meanE = mean(sims_obj$data$E), meanY = mean(sims_obj$data$Y), meanU = mean(sims_obj$data[[BC]]),
      meanU_e0 = mean(sims_obj$data[[BC]][sims_obj$data$E == 0]), meanU_e1 = mean(sims_obj$data[[BC]][sims_obj$data$E == 1]),
      OReyFull = exp(full_mod[["E"]]), OReyPart = exp(reduced_mod[["E"]]), parm_totals, OReyPartAdj = totalAdjOREY, 
      bOReydiff = totalAdjOREY - exp(full_mod[["E"]]))
}












