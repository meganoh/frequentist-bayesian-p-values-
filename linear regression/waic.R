waic <- print(waic_oosterwijk)
waic[0]
waic = pm.waic(trace, model)
data.frame(waic_oosterwijk$loos[1])
df <- waic_oosterwijk$loo[1] %>% tibble()
waic_oosterwijk$bayes_oosterwijkprior
