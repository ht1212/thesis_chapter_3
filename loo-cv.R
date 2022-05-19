
#-------------------------MODEL COMPARISONS-------------------------------------

#extract pointwise log likelihood 
ll1 <- extract_log_lik(model.4, parameter_name = "log_lik", merge_chains = FALSE)
ll2 <- extract_log_lik(model.5, parameter_name = "log_lik", merge_chains = FALSE)
ll3 <- extract_log_lik(model.6, parameter_name = "log_lik", merge_chains = FALSE)
#ll4 <- extract_log_lik(model.4, parameter_name = "log_lik", merge_chains = FALSE)

# sample size and monte carlo error 
r_eff_ll1 <- relative_eff(exp(ll1), cores = 4)
r_eff_ll2 <- relative_eff(exp(ll2), cores = 4)
r_eff_ll3 <- relative_eff(exp(ll3), cores = 4)
#r_eff_ll4 <- relative_eff(exp(ll4), cores = 4)

# loo 
loo_1 <- loo(ll1, r_eff = r_eff_ll1, cores = 4)
loo_2 <- loo(ll2, r_eff = r_eff_ll2, cores = 4)
loo_3 <- loo(ll3, r_eff = r_eff_ll3, cores = 4)
#loo_4 <- loo(ll4, r_eff = r_eff_ll4, cores = 4)

# Compare
comp <- loo_compare(loo_1, loo_2, loo_3)#, loo_4) #  
print(comp, simplify=FALSE)


# WAIC 
# waic_1 <- waic(ll1)
# waic_2 <- waic(ll2)
# waic_3 <- waic(ll3)
# waic_4 <- waic(ll4)
# 
# 
# waic_diff <- loo_compare(waic_1, waic_2, waic_3, waic_4)
# 
# print(waic_diff, simplify = FALSE)
