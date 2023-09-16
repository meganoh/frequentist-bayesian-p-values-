#2 group simulation
fit_2groups <- function(i, sample_size, mod_iter, mod_warmup) {
  
  data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                     value = c(rnorm(n = sample_size), 
                               rnorm(n = sample_size)))
  
  data$group <- factor(data$group)
  
  means <- data %>% group_by(group) %>% 
    summarise(mean = mean(value)) %>% 
    pivot_wider(names_from = group, values_from = mean) %>% 
    mutate(diff = treatmentA - control)
  diff <- means$diff
  
  var <- data %>% group_by(group) %>% 
    summarise(sd = sd(value)) %>% 
    pivot_wider(names_from = group, values_from = sd) %>% 
    mutate(sd = sqrt((control^2 + treatmentA^2)/2)) %>% 
    mutate(se = sd * sqrt(2/sample_size))
  se <- var$se
  
  #freq model
  
  freq <- lm(value ~ group, data)
  
  freq_pval <- summary(freq)$coef[2, 4]

  #bayes model - null 
  
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - flat prior 
  bayes_flatprior <- update(bayes_flatprior_prefit, 
                            newdata = data,  recompile = FALSE)
  
  bayes_flat_test <- joint_tests(bayes_flatprior)
  
  
  bayes_flat_pval <- bayes_flat_test$p.value
  
  #bayes model - student prior 
  bayes_studentprior <- update(bayes_studentprior_prefit, 
                               newdata = data, recompile = FALSE)
  
  
  bayes_student_test <- joint_tests(bayes_studentprior)
  
  
  bayes_student_pval <- bayes_student_test$p.value
  
  #bayes model - oosterwijk prior 
  bayes_oosterwijkprior <- update(bayes_oosterwijkprior_prefit, 
                                  newdata = data, recompile = FALSE)
  bayes_oosterwijk_test <- joint_tests(bayes_oosterwijkprior)
  bayes_oosterwijk_pval <- bayes_oosterwijk_test$p.value
  
  #LOO & WAIC comparison
  #flat
  loo_flat <- LOO(bayes_flatprior, bayes_intercept)
  loo_elpd_diff_flat <- loo_flat$diffs[2]
  loo_se_diff_flat <- loo_flat$diffs[4]
  
  waic_flat <- WAIC(bayes_flatprior, bayes_intercept)
  waic_elpd_diff_flat <- waic_flat$diffs[2]
  waic_se_diff_flat <- waic_flat$diffs[4]
  
  #student
  loo_student <- LOO(bayes_studentprior, bayes_intercept)
  loo_elpd_diff_student <- loo_student$diffs[2]
  loo_se_diff_student <- loo_student$diffs[4]
  
  waic_student <- WAIC(bayes_studentprior, bayes_intercept)
  waic_elpd_diff_student <- waic_student$diffs[2]
  waic_se_diff_student <- waic_student$diffs[4]
  
  #oosterwijk
  loo_oosterwijk <- LOO(bayes_oosterwijkprior, bayes_intercept)
  loo_elpd_diff_oosterwijk <- loo_oosterwijk$diffs[2]
  loo_se_diff_oosterwijk <- loo_oosterwijk$diffs[4]
  
  waic_oosterwijk <- WAIC(bayes_oosterwijkprior, bayes_intercept)
  waic_elpd_diff_oosterwijk <- waic_oosterwijk$diffs[2]
  waic_se_diff_oosterwijk <- waic_oosterwijk$diffs[4]
  
  
  #bayes factors 
  #flat 
  bf_flat_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                               bayes_flatprior, log = TRUE)
  
  
  bf_flat <- bf_flat_test$bf
  
  #student
  
  
  bf_student_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                  bayes_studentprior, log = TRUE)
  
  
  bf_student <- bf_student_test$bf
  
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
  p_flat <- pchisq(diff_flat, df = 1, lower.tail = FALSE)
  p_flat_max <- pchisq(diff_flat_max, df = 1, lower.tail = FALSE)
  
  #student
  student <- as.data.frame(bayes_studentprior$fit)
  student <- student %>% 
    mutate(
      lp_student_new = student$lp__ - 
        dlst(x = student$b_group1, df = 3, sigma = 0.5, mu = 0, log = TRUE) 
    )
  lp_student <- mean(student$lp_student_new)
  lp_student_max <- max(student$lp_student_new)
  
  diff_student <- lp_student - lp_null
  diff_student_max <- lp_student_max - lp_null_max
  
  p_student <- pchisq(diff_student, df = 1, lower.tail = FALSE)
  p_student_max <-pchisq(diff_student_max, df = 1, lower.tail = FALSE)
  
  #oosterwijk
  oosterwijk <- as.data.frame(bayes_oosterwijkprior$fit)
  oosterwijk <- oosterwijk %>% 
    mutate(
      lp_oosterwijk_new = oosterwijk$lp__ -
        dlst(x = oosterwijk$b_group1, df = 3, sigma = 0.102, mu = 0.350, log = TRUE)
    )
  lp_oosterwijk <- mean(oosterwijk$lp_oosterwijk_new)
  lp_oosterwijk_max <- max(oosterwijk$lp_oosterwijk_new)
  
  diff_oosterwijk <- lp_oosterwijk - lp_null
  diff_oosterwijk_max <- lp_oosterwijk_max - lp_null_max
  
  p_oosterwijk <- pchisq(diff_oosterwijk, df = 1, lower.tail = FALSE)
  p_oosterwijk_max <-pchisq(diff_oosterwijk_max, df = 1, lower.tail = FALSE)
  
  out_freq <- data.frame(n = sample_size, diff, se, 
                         test = "freq", pval = freq_pval)
  
  out_flatbayes <- data.frame(n = sample_size, diff, se, 
                              test = "bayes_flatprior", pval = bayes_flat_pval, 
                              loo_diff = loo_elpd_diff_flat, 
                              loo_se = loo_se_diff_flat,
                              waic_diff = waic_elpd_diff_flat,
                              waic_se = waic_se_diff_flat,
                              bf = bf_flat,
                              lp_diff = diff_flat,
                              lp_diff_max = diff_flat_max,
                              p_likelihood = p_flat,
                              p_likelihood_max = p_flat_max)
  out_studentbayes <- data.frame(n = sample_size, diff, se, 
                                 test = "bayes_studentprior", pval = bayes_student_pval, 
                                 loo_diff = loo_elpd_diff_student, 
                                 loo_se = loo_se_diff_student,
                                 waic_diff = waic_elpd_diff_student,
                                 waic_se = waic_se_diff_student,
                                 bf = bf_student,
                                 lp_diff = diff_student,
                                 lp_diff_max = diff_student_max,
                                 p_likelihood = p_student,
                                 p_likelihood_max = p_student_max)
  out_oosterwijkbayes <- data.frame(n = sample_size, diff, se, 
                                    test = "bayes_oosterwijkprior", pval = bayes_oosterwijk_pval,
                                    loo_diff = loo_elpd_diff_oosterwijk, 
                                    loo_se = loo_se_diff_oosterwijk,
                                    waic_diff = waic_elpd_diff_oosterwijk,
                                    waic_se = waic_se_diff_oosterwijk,
                                    bf = bf_oosterwijk,
                                    lp_diff = diff_oosterwijk,
                                    lp_diff_max = diff_oosterwijk_max,
                                    p_likelihood = p_oosterwijk, 
                                    p_likelihood_max = p_oosterwijk_max) 
  out <- bind_rows(out_freq,
                   out_flatbayes, 
                   out_studentbayes, 
                   out_oosterwijkbayes)
  return(out)
}

run_2groups <-  function(iter, sample_size, mod_iter, mod_warmup) {
  results_2groups <- mclapply(X = 1:iter, 
                              FUN = fit_2groups, sample_size, mod_iter, mod_warmup, 
                              mc.cores = 16, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(results_2groups, file = paste0("simulation_2groups", "_n", sample_size, "_iter", iter, ".rda"))
}