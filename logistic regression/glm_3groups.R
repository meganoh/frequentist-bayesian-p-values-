glm_3groups <- function(i, sample_size, mod_iter, mod_warmup){
  prob <- rbeta(1, 2, 2)
  
  data <- data.frame(group = rep(c("control", "treatmentA", "treatmentB"), each = sample_size), 
                     value = c(rbinom(n = sample_size, size = 1, prob = prob), 
                               rbinom(n = sample_size, size = 1, prob = prob), 
                               rbinom(n = sample_size, size = 1, prob = prob)))
  data <- data %>% mutate(id = row_number())
  data$group <- factor(data$group)
  
  #freq model
  freq <- glm(value ~ group, data, family=binomial(link='logit'))
  freq_test <- joint_tests(freq)
  freq_pval <- freq_test$p.value
  freq_pval_df <- pf(freq_test$F.ratio, freq_test$df1, df2 = sample_size*2, lower.tail = FALSE)

  freq_df <- freq %>% tidy()
  #null model 
  bayes_intercept <- update(bayes_intercept_prefit, 
                            newdata = data, recompile = FALSE)
  
  #bayes model - flat prior 
  bayes_flatprior <- update(bayes_flatprior_prefit, 
                            newdata = data,  recompile = FALSE)
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
                              bf = bf_flat, 
                              joint_tests = list(bayes_flat_test))
  out_tighterbayes <- data.frame(n = sample_size, 
                                 test = "bayes_tighterprior", pval = bayes_tighter_pval, pval_df = bayes_tighter_pval_df, 
                                 bf = bf_tighter, 
                                 joint_tests = list(bayes_tighter_test))
  out_widerbayes <- data.frame(n = sample_size, 
                               test = "bayes_widerprior", pval = bayes_wider_pval, pval_df = bayes_wider_pval_df, 
                               bf = bf_wider, 
                               joint_tests = list(bayes_wider_test))
  out <- bind_rows(out_freq,
                   out_flatbayes, 
                   out_tighterbayes, 
                   out_widerbayes)
  out = out %>% mutate(frequentist = list(freq_df))
  
  return(out)
  Sys.sleep(2)
}

glmrun_3groups <-  function(iter, sample_size, mod_iter, mod_warmup, cores) {
  glmresults_3groups <- mclapply(X = 1:iter, 
                              FUN = glm_3groups, sample_size, mod_iter, mod_warmup, 
                              mc.cores = cores, mc.preschedule = FALSE, mc.cleanup = TRUE)
  save(glmresults_3groups, file = paste0("glm_3groups", "_n", sample_size, "_iter", iter, ".rda"))
}