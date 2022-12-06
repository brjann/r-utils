library(nlme)
library(boot)
getCI <- function(boots){
  res <- NULL
  for(i in 1:length(boots$t0)){
    ci.res <- boot.ci(boots, index = i, type = c("norm", "basic", "perc", "bca"))
    for(CI in c("normal", "basic", "percent", "bca")){
      x <- tail(ci.res[[CI]][1, ], 2)
      low = x[1]
      high = x[2]
      sig = ifelse(sign(high) == sign(low), '*', '')
      res <- rbind(res, data.frame(coeff = names(boots$t0)[i], CI = CI, est = boots$t0[i], low = low, high = high, sig = sig))
    }
  }
  row.names(res) <- NULL
  res
}

# ---- MULTIPLE MEDIATOR FUNCTIONS ----

# mediation_multi_boot <- function(X, Y, Ms, X_formula, model, ...){
#   function(data, index){
#     mediation_multi(data[index, ], X, Y, Ms, X_formula, model_safe(model), ...)["est"]
#   }
# }
# 
mediation_multi <- function(data, X, Y, Ms, X_formula, model, ...){
  f_b <- as.formula(paste(Y, "~", X_formula, "+", paste(Ms, collapse = "+")))
  m_b <- model(f_b, data, ...)
  f_c <- as.formula(paste(Y, "~", X_formula))
  m_c <- model(f_c, data, ...)
  
  CPRIM <- coef.stats(m_b, X)
  C <- coef.stats(m_c, X)
  c_ <- C["est"]
  
  ab_multi <- 0
  res <- data.frame()
  for(M in Ms){
    f_a <- as.formula(paste(M, "~", X_formula))
    m_a <- model(f_a, data, ...)
    abc <- get_ab(m_a, m_b, c_, X, M)
    # l_m_b[[M]] <- m_a
    ab_multi <- ab_multi + abc[3, "est"]
    res <- rbind(res, abc)
  }
  AB_MULTI = c(est = ab_multi, se = NA, CI_low = NA, CI_up = NA, pval = NA)
  AB_MULTI_PERC = c(est = ab_multi / c_, se = NA, CI_low = NA, CI_up = NA, pval = NA)
  res <- rbind(res, data.frame(rbind(C, CPRIM, AB_MULTI, AB_MULTI_PERC)))
  
  res
}

sobel.test <- function(a, a_se, b, b_se){
  sqrt(b^2 * a_se^2 + a^2 * b_se^2 + a_se^2 * b_se^2)
}

get_ab <- function(m_a, m_b, c_, X, M){
  A <- coef.stats(m_a, X)
  B <- coef.stats(m_b, M)
  
  a <- A["est"]
  b <- B["est"]
  a_se <- A["se"]
  b_se <- B["se"]

  ab <- a * b
  ab_se <- sobel.test(a, a_se, b, b_se)
  ab_pval <- (1 - pnorm(abs(ab/ab_se))) * 2
  perc <- ab / c_
  AB = c(est = ab, se = ab_se, CI_low = NA, CI_up = NA, pval = ab_pval)
  PERC <- c(est = perc, se = NA, CI_low = NA, CI_up = NA, pval = NA)
  res <- data.frame(rbind(A, B, AB, PERC))
  row.names(res) <- paste0(M, "_", row.names(res))

  res
}

coef.stats <- function(fit, coef){
  est <- fixef(fit)[coef]
  se <- sqrt(diag(fit$varFix))[coef]
  CI_low <- intervals(fit, which = "fixed")$fixed[coef, "lower"]
  CI_up <- intervals(fit, which = "fixed")$fixed[coef, "upper"]
  pval <- summary(fit)$tTable[coef,"p-value"]
  named_c(est, se, CI_low, CI_up, pval)
}

named_c <- function(...) {
  dots <- substitute(list(...))[-1]
  res <- c(...)
  names(res) <- c(sapply(dots, deparse))
  res
}

mediation_add_boot <- function(abc, boot_out, type = "bca"){
  CIs <- NULL
  for(i in 1:nrow(abc)){
    ci <- boot.ci(boot_out, index = i, type = type)
    CIs <- rbind(CIs, data.frame("boot_CI_low" = ci[[4]][4], "boot_CI_up" = ci[[4]][5]))
  }
  cbind(abc, CIs)
}