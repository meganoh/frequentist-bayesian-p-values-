#3 group simulation
fit_3groups <- function(i, sample_size, mod_iter, mod_warmup) {
  data <- data.frame(group = rep(c("control", "treatmentA", "treatmentB"), each = sample_size), 
                     value = c(rnorm(n = sample_size), 
                               rnorm(n = sample_size), 
                               rnorm(n = sample_size)))
  data$group <- factor(data$group)
  
  #freq model
  freq <- lm(value ~ group, data)
  freq_test <- joint_tests(freq)
  freq_pval <- freq_test$p.value
  
  #bayes model - null 
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - flat prior 
  bayes_flatprior <- update(bayes_flatprior_prefit, 
                            newdata = data,  recompile = FALSE)
  bayes_flat_test <- joint_tests(bayes_flatprior)
  bayes_flat_pval <- bayes_flat_test$p.value
  
  #bayes model - wide prior 
  bayes_wideprior <- update(bayes_wideprior_prefit, 
                            newdata = data, recompile = FALSE)
  
  bayes_wide_test <- joint_tests(bayes_wideprior)
  
  bayes_wide_pval <- bayes_wide_test$p.value
  
  
  #bayes model - tight prior 
  bayes_tightprior <- update(bayes_tightprior_prefit, 
                             newdata = data, recompile = FALSE)
  
  bayes_tight_test <- joint_tests(bayes_tightprior)
  
  bayes_tight_pval <- bayes_tight_test$p.value
  
  #bayes model - oosterwijk prior 
  bayes_oosterwijkprior <- update(bayes_oosterwijkprior_prefit, 
                                  newdata = data, recompile = FALSE)
  bayes_oosterwijk_test <- joint_tests(bayes_oosterwijkprior)
  bayes_oosterwijk_pval <- bayes_oosterwijk_test$p.value
  
  #bayes factors 
  #flat 
  bf_flat_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                               bayes_flatprior, log = TRUE)
  bf_flat <- bf_flat_test$bf
  
  #wider prior 
  bf_wide_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                               bayes_wideprior, log = TRUE)
  bf_wide <- bf_wide_test$bf
  
  #tighter prior 
  bf_tight_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                bayes_tightprior, log = TRUE)
  bf_tight <- bf_tight_test$bf
  
  #oosterwijk
  bf_oosterwijk_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                     bayes_oosterwijkprior, log = TRUE)
  bf_oosterwijk <- bf_oosterwijk_test$bf
  
  
  #likelihood ratio test 
  null <- as.data.frame(bayes_intercept$fit)
  lp_null <- mean(null$lp__)
  lp_null_max <- max(null$lp__)
  
  #flat 
  flat <- as.data.frame(bayes_flatprior$fit)
  lp_flat <- mean(flat$lp__)
  lp_flat_max <- max(flat$lp__)
  
  diff_flat <- lp_flat - lp_null
  diff_flat_max <- lp_flat_max - lp_null_max
  p_flat <- pchisq(diff_flat, df = 2, lower.tail = FALSE)
  p_flat_max <- pchisq(diff_flat_max, df = 2, lower.tail = FALSE)
  
  #wide
  wide <- as.data.frame(bayes_wideprior$fit)
  wide <- wide %>% 
    mutate(
      lp_wide_new = wide$lp__ - 
        dlst(x = wide$b_group1, df = 3, sigma = 0.5, mu = 0, log = TRUE) 
    )
  lp_wide <- mean(wide$lp_wide_new)
  lp_wide_max <- max(wide$lp_wide_new)
  
  diff_wide <- lp_wide - lp_null
  diff_wide_max <- lp_wide_max - lp_null_max
  
  p_wide <- pchisq(diff_wide, df = 1, lower.tail = FALSE)
  p_wide_max <-pchisq(diff_wide_max, df = 1, lower.tail = FALSE)
  
  #tight
  tight <- as.data.frame(bayes_tightprior$fit)
  tight <- tight %>% 
    mutate(
      lp_tight_new = tight$lp__ - 
        dlst(x = tight$b_group1, df = 3, sigma = 0.5, mu = 0, log = TRUE) 
    )
  lp_tight <- mean(tight$lp_tight_new)
  lp_tight_max <- max(tight$lp_tight_new)
  
  diff_tight <- lp_tight - lp_null
  diff_tight_max <- lp_tight_max - lp_null_max
  
  p_tight <- pchisq(diff_tight, df = 1, lower.tail = FALSE)
  p_tight_max <-pchisq(diff_tight_max, df = 1, lower.tail = FALSE)
  
  #oosterwijk
  oosterwijk <- as.data.frame(bayes_oosterwijkprior$fit)
  oosterwijk <- oosterwijk %>% 
    mutate(
      lp_oosterwijk_new = oosterwijk$lp__ -
        dlst(x = oosterwijk$b_group1, df = 3, sigma = 0.102, mu = 0.350, log = TRUE) -
        dlst(x = oosterwijk$b_group2, df = 3, sigma = 0.102, mu = 0.350, log = TRUE) 
    )
  
  lp_oosterwijk <- mean(oosterwijk$lp_oosterwijk_new)
  lp_oosterwijk_max <- max(oosterwijk$lp_oosterwijk_new)
  
  diff_oosterwijk <- lp_oosterwijk - lp_null
  diff_oosterwijk_max <- lp_oosterwijk_max - lp_null_max
  
  p_oosterwijk <- pchisq(diff_oosterwijk, df = 2, lower.tail = FALSE)
  p_oosterwijk_max <-pchisq(diff_oosterwijk_max, df = 2, lower.tail = FALSE)
  
  out_freq <- data.frame(n = sample_size,
                         test = "freq", pval = freq_pval)
  out_flatbayes <- data.frame(n = sample_size,
                              test = "bayes_flatprior", pval = bayes_flat_pval, 
                              bf = bf_flat,
                              lp_diff = diff_flat,
                              lp_diff_max = diff_flat_max,
                              p_likelihood = p_flat,
                              p_likelihood_max = p_flat_max)
  out_widebayes <- data.frame(n = sample_size, 
                              test = "bayes_wideprior", pval = bayes_wide_pval, 
                              bf = bf_wide,
                              lp_diff = diff_wide,
                              lp_diff_max = diff_wide_max,
                              p_likelihood = p_wide,
                              p_likelihood_max = p_wide_max)
  out_tightbayes <- data.frame(n = sample_size, 
                               test = "bayes_tightprior", pval = bayes_tight_pval, 
                               bf = bf_tight,
                               lp_diff = diff_tight,
                               lp_diff_max = diff_tight_max,
                               p_likelihood = p_tight,
                               p_likelihood_max = p_tight_max)
  out_oosterwijkbayes <- data.frame(n = sample_size, 
                                    test = "bayes_oosterwijkprior", pval = bayes_oosterwijk_pval, 
                                    bf = bf_oosterwijk,  
                                    lp_diff = diff_oosterwijk,
                                    lp_diff_max = diff_oosterwijk_max,
                                    p_likelihood = p_oosterwijk,
                                    p_likelihood_max = p_oosterwijk_max) 
  out <- bind_rows(out_freq,
                   out_flatbayes, 
                   out_widebayes, 
                   out_tightbayes,
                   out_oosterwijkbayes)
  Sys.sleep(1)
  return(out)
}

run_3groups <-  function(iter, sample_size, mod_iter, mod_warmup, cores) {
  results_3groups <- mclapply(X = 1:iter, 
                              FUN = fit_3groups, sample_size, mod_iter, mod_warmup, 
                              mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(results_3groups, file = paste0("simulation_3groups", "_n", sample_size, "_iter", iter, ".rda"))
}



