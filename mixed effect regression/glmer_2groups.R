source("convergence.R")

glmer_2groups <- function(i, sample_size, re_sd, mod_iter, mod_warmup){
  prob <- qlogis(rbeta(1, 2, 2))

  data <- data.frame(group = rep(c("control", "treatmentA"), each = sample_size), 
                     prob = prob,
                     re_ppt = rnorm(n = sample_size*2, 0, re_sd), 
                     n = 25) #number of trials 
  
  data <- data %>% mutate(x = plogis(prob + re_ppt)) #exp log the probability range from 0 to 1 
  data <- data %>% mutate(y = rbinom(n = nrow(data), size = n, prob = x))
  
  data <- data %>% mutate(id = row_number())
  data$group <- factor(data$group)
  
  freq <- glmer(cbind(y, n-y) ~ group + (1|id), data, family=binomial())
  freq_test <- joint_tests(freq)
  freq_pval <- freq_test$p.value
  freq_pval_df <- pf(freq_test$F.ratio, freq_test$df1, df2 = sample_size*2, lower.tail = FALSE)
  
  conv <- mod_check(freq)
  freq_df <- freq %>% tidy()
  
  #null model 
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - flat prior 
  bayes_flatprior <- update(bayes_flatprior_prefit, 
                            newdata = data, recompile = FALSE)
  bayes_flat_test <- joint_tests(bayes_flatprior)
  bayes_flat_pval <- bayes_flat_test$p.value
  bayes_flat_pval_df <- pf(bayes_flat_test$F.ratio, bayes_flat_test$df1, df2 = sample_size*2, lower.tail = FALSE)
  
  #bayes model - tighter .2 prior 
  bayes_tighterprior <- update(bayes_tighterprior_prefit, 
                               newdata = data, recompile = FALSE)
  bayes_tighter_test <- joint_tests(bayes_tighterprior)
  bayes_tighter_pval <- bayes_tighter_test$p.value
  bayes_tighter_pval_df <- pf(bayes_tighter_test$F.ratio, bayes_tighter_test$df1, df2 = sample_size*2, lower.tail = FALSE)
  
  #bayes model - wider .5 prior 
  bayes_widerprior <- update(bayes_widerprior_prefit, 
                             newdata = data, recompile = FALSE)
  bayes_wider_test <- joint_tests(bayes_widerprior)
  bayes_wider_pval <- bayes_wider_test$p.value
  bayes_wider_pval_df <- pf(bayes_wider_test$F.ratio, bayes_wider_test$df1, df2 = sample_size*2, lower.tail = FALSE)
  
  #LOO & WAIC comparison
  #flat
  loo_flat <- LOO(bayes_intercept, bayes_flatprior)
  loo_elpd_diff_flat <- loo_flat$diffs[2]
  loo_se_diff_flat <- loo_flat$diffs[4]
  
  waic_flat <- WAIC(bayes_intercept, bayes_flatprior)
  waic_elpd_diff_flat <- waic_flat$diffs[2]
  waic_se_diff_flat <- waic_flat$diffs[4]
  
  #tighter
  loo_tighter <- LOO(bayes_intercept, bayes_tighterprior)
  loo_elpd_diff_tighter <- loo_tighter$diffs[2]
  loo_se_diff_tighter <- loo_tighter$diffs[4]
  
  waic_tighter <- WAIC(bayes_intercept, bayes_tighterprior)
  waic_elpd_diff_tighter <- waic_tighter$diffs[2]
  waic_se_diff_tighter <- waic_tighter$diffs[4]
  
  #wider
  loo_wider <- LOO(bayes_intercept, bayes_widerprior)
  loo_elpd_diff_wider <- loo_wider$diffs[2]
  loo_se_diff_wider <- loo_wider$diffs[4]
  
  waic_wider <- WAIC(bayes_intercept, bayes_widerprior)
  waic_elpd_diff_wider <- waic_wider$diffs[2]
  waic_se_diff_wider <- waic_wider$diffs[4]
  
  #bayes factors 
  #flat 
  bf_flat_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                               bayes_flatprior, log = TRUE)
  bf_flat <- bf_flat_test$bf
  
  #tighter
  bf_tighter_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                  bayes_tighterprior, log = TRUE)
  bf_tighter <- bf_tighter_test$bf
  
  #wider
  bf_wider_test <- bridgesampling::bayes_factor(bayes_intercept, 
                                                bayes_widerprior, log = TRUE)
  bf_wider <- bf_wider_test$bf
  
  out_freq <- data.frame(n = sample_size,
                         test = "freq", pval = freq_pval, pval_df = freq_pval_df, 
                         joint_tests = list(freq_test))
  out_flatbayes <- data.frame(n = sample_size, 
                              test = "bayes_flatprior", pval = bayes_flat_pval,  pval_df = bayes_flat_pval_df,
                              loo_diff = loo_elpd_diff_flat, 
                              loo_se = loo_se_diff_flat,
                              waic_diff = waic_elpd_diff_flat,
                              waic_se = waic_se_diff_flat,
                              bf = bf_flat, 
                              joint_tests = list(bayes_flat_test))
  out_tighterbayes <- data.frame(n = sample_size, 
                                 test = "bayes_tighterprior", pval = bayes_tighter_pval, pval_df = bayes_tighter_pval_df, 
                                 loo_diff = loo_elpd_diff_tighter, 
                                 loo_se = loo_se_diff_tighter,
                                 waic_diff = waic_elpd_diff_tighter,
                                 waic_se = waic_se_diff_tighter,
                                 bf = bf_tighter, 
                                 joint_tests = list(bayes_tighter_test))
  out_widerbayes <- data.frame(n = sample_size, 
                               test = "bayes_widerprior", pval = bayes_wider_pval, pval_df = bayes_wider_pval_df, 
                               loo_diff = loo_elpd_diff_wider, 
                               loo_se = loo_se_diff_wider,
                               waic_diff = waic_elpd_diff_wider,
                               waic_se = waic_se_diff_wider,
                               bf = bf_wider, 
                               joint_tests = list(bayes_wider_test)) 
  out <- bind_rows(out_freq,
                   out_flatbayes, 
                   out_tighterbayes, 
                   out_widerbayes)
  
  out$freqConv = conv
  out = out %>% mutate(frequentist = list(freq_df))
  return(out)
  Sys.sleep(2)
  
}

glmerrun_2groups <-  function(iter, sample_size, re_sd, mod_iter, mod_warmup, cores) {
  glmerresults_2groups <- mclapply(X = 1:iter, 
                              FUN = glmer_2groups, sample_size,  re_sd, mod_iter, mod_warmup, 
                              mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(glmerresults_2groups, file = paste0("glmerresults_2groups", "_n", sample_size, "_iter", iter, "sd", re_sd, ".rda"))
}
