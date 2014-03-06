#' Generates 2-way probability tables for an exposure and an outcome given a confounder
#' 
#' This function returns a table of probabilities when supplied with prevalences, odds ratios and interaction strength
#' 
#' E = exposure
#' Y = Outcome
#' U = the (unmeasured) confounder
#' 
#' @export
#' @param p_u1 Prevalence of the UC (u)
#' @param p_e1_u0 <- 0.5 # prevalence of e when u=0
#' @param p_y1_u0 <- 0.2 # prevalence of y when u=0
#' @param S <- 3 # Strength of the interaction term
#' @param or_ey_u0 <- 4 # Odds-ratio for ey when u=0
#' @param or_uy_e0 <- 3 # Odds-ratio for uy when e=0
#' @param or_ue_y0 <- 9 # Odds-ratio for ue when y=0
#' @param N <- 100 # Population size
two_way_table <- function(p_u1, p_e1_u0, p_y1_u0, S, or_ey_u0, or_uy_e0, or_ue_y0, N){
    n_u0 <- N * (1 - p_u1)
    n_u1 <- N * p_u1
    
    x <- 1 - or_ey_u0
    y <- n_u0 * (1 - p_e1_u0 - p_y1_u0 + or_ey_u0 * p_e1_u0 + or_ey_u0 * p_y1_u0)
    z <- -n_u0 * n_u0 * or_ey_u0 * p_e1_u0 * p_y1_u0
    
    d0 <- (-y + sqrt(y * y - 4 * x * z)) / (2 * x)
    b0 <- p_y1_u0 * n_u0 - d0
    c0 <- p_e1_u0 * n_u0 - d0
    a0 <- n_u0 - b0 - c0 - d0
    
    or_ey_u1 <- S * or_ey_u0
    or_uy_e1 <- S * or_uy_e0
    or_ue_y1 <- S * or_ue_y0
    
    X <- or_uy_e0 * b0 + or_ue_y0 * c0 + or_ue_y0 * S * or_uy_e0 * d0
    a1 <- n_u1 / (1 + X / a0)
    b1 <- a1 * or_uy_e0 * b0 / a0
    c1 <- a1 * or_ue_y0 * c0 / a0
    d1 <- a1 * or_ue_y0 * S * or_uy_e0 * d0 / a0
    
    as.table(array(c(a0,c0,b0,d0,a1,c1,b1,d1,
            a0+a1, c0+c1, b0+b1, d0+d1), 
          dim=c(2,2,3),
          dimnames = list(c("E=0", "E=1"), 
                          c("Y=0", "Y=1"), 
                          c("U=0", "U=1", "Total"))))
}

