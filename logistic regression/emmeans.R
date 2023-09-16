

### CODE BEGINS ###



library("emmeans")

library("brms")

library("bayestestR")

options(contrasts = c("contr.equalprior_deviations", "contr.poly"))



data("Machines", package = "MEMSS")



fit1 <- brm(score ~ Machine + (Machine|Worker), data=Machines)



# generate a summary of the results

summary(fit1)



em1 <- emmeans(fit1, c("Machine"))

em1



em1@V



summary(em1, frequentist = TRUE)



test(contrast(em1, "consec"), joint = TRUE)

joint_tests(fit1)


### CODE ENDS ###


