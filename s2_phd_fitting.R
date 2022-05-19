
##************************************************************************************
# script to fit bayesian survival models to Zongo et al survival data through rStan * 
# and produce plots to compare model estimated survival and efficacy to trial data  *
#************************************************************************************
# H.Thompson 
# Jan 2022  
# Imperial College London 
# PhD Chapter 3 code 
# ht1212@imperial.ac.uk 

#---------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------set-up: load data and packages----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
source("packages.R")

set.seed(27)
options(mc.cores = parallel::detectCores())
color_scheme_set("darkgray")
#color_scheme_view()
bayesplot_theme_update(text = element_text(size=10, family = "sans"))

source("s1_load_data.R")
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------single exponential model with modifiers--------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("model.1.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, stan_data_llin_age_site)

#---------------------------FITTING---------------------------------------------
model.1 <- stan("model.1_new.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.1, pars=c("alpha","lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)

saveRDS(model.1, "model.1.RDS")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------monthly exponential model with modifiers-------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("s1_load_data.R")

source("model.2.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, split_data, split_data_2, split_data_3)

#---------------------------FITTING---------------------------------------------
model.2 <- stan("model.2_new.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.2, pars=c("alpha","lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)

saveRDS(model.2, "model.2.RDS")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------biweekly exponential model with modifiers-------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("s1_load_data.R")

source("model.3.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, split_data, split_data_2, split_data_3)

#---------------------------FITTING---------------------------------------------
model.3 <- stan("model.3_new.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.3, pars=c("alpha","lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)

saveRDS(model.3, "model.3.RDS")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------weekly exponential model with modifiers--------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("s1_load_data.R")

source("model.4.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, split_data, split_data_2, split_data_3)

#---------------------------FITTING---------------------------------------------
model.4 <- stan("model.4_new.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.4, pars=c("alpha","lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)

saveRDS(model.4, "model.4.RDS")

#-------------------------------------------------------------------------------------------------------------------------------------
#----------------------------plot the comparisons-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
source("model_plot_functions.R")

# model 1 predicted curves
source("model.1_plot.R")
fig_1a <- survival_plot(model.1_predited) + labs(title="A")
fig_1a 

# model 2 predicted curves
source("model.2_plot_new.R")
fig_1b <- survival_plot(model.2_predicted) + labs(title="B")
fig_1b 

# model 3 predicted curves
source("model.3_plot.R")
fig_1c <- survival_plot(model.3_predicted) + labs(title="C")
fig_1c 

# model 4 predicted curves
source("model.4_plot.R")
fig_1d <- survival_plot(model.4_predicted) + labs(title="D")
fig_1d 

# model 1 efficacy curve 
fig_1e <- duration_protection(model.1) + labs(title="E")
fig_1e  

# model 2 efficacy curve
fig_1f <- duration_protection(model.2) + labs(title="F")
fig_1f

# model 3 efficacy curve
fig_1g <- duration_protection(model.3) + labs(title="G")
fig_1g

# model 4 efficacy curve
fig_1h <- duration_protection(model.4) + labs(title="H")
fig_1h

#-combine into single plot---------
col_1 <-
  fig_1a + fig_1b + fig_1c + fig_1d  +
 # plot_annotation(tag_levels = list(c("A", "B", "C", "D")))+
  plot_layout(ncol=1, guides = "collect") &
  theme(legend.position = 'bottom', plot.tag = element_text(size = 14))

col_1


col_2 <-
  fig_1e + fig_1f + fig_1g + fig_1h  +
#  plot_annotation(tag_levels = list(c("E", "F", "G", "H")))+
  plot_layout(ncol=1) &
  theme(plot.tag = element_text(size = 14))

col_2

col_1 | col_2 #+
#  plot_annotation(tag_levels = "A")&
 # theme(plot.tag = element_text(size = 14))

ggsave("fig_1.png", width=9.5, height=14, dpi=600)

#-------------------------------------------------------------------------------------------------------------------------------------
#----------------------------select best fitting model?-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------

# through splitting the data into multiple observations per patient we can't use stan built in leave one out cross validation 
# due to the data sets now not being the same size. So decisions for the first sets of models are based on observed fitting to the 
# trial data. 

# model comparisons then examine model.4 and different functional forms of SP+AQ efficacy profiles 

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------weekly exponential model with modifiers---------------------------------------------------------------------
#---------------------------- Gamma functional form SP+AQ--------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("s1_load_data.R")

source("model.4.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, split_data, split_data_2, split_data_3)

#---------------------------FITTING---------------------------------------------
model.5 <- stan("model.5.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.5, pars=c("alpha","theta", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)


saveRDS(model.5, "model.5.RDS")

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------weekly exponential model with modifiers---------------------------------------------------------------------
#---------------------------- Hill function for SP+AQ------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
source("s1_load_data.R")

source("model.4.R")

rm(control_data, agg_data, data, spaq_data, stan_data, stan_data_llin, stan_data_plain, split_data, split_data_2, split_data_3)

#---------------------------FITTING---------------------------------------------
model.6 <- stan("model.6.stan",
                data = stan_data_llin_age,
                chains = 4,
                iter = 5000)

print(model.6, pars=c("alpha","l", "r", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),  digits_summary = 3)


saveRDS(model.6, "model.6.RDS")

#------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------LOO-CV WITH LOO PACKAGE----------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
source("loo-cv.R")

# weibull model - there is very small differences suggesting all models capture the decay in protection well but will select the weibull 
# models.  


#-----------------------------------------------------------------------------------------------------------------------------------------
#------------------------------Trace plots and histograms for chapter---------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
post_model.4 <- as.array(model.4)

# trace plot
mcmc_trace(post_model.4, pars=c("alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]", "alpha[5]", "alpha[6]","alpha[7]","alpha[8]","alpha[9]", 
                                "lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"), 
            facet_args = list(strip.position = "top",nrow=5))

ggsave("model.4_trace.png", dpi=600, width=8, height=7)


# parameter histogram
mcmc_hist(post_model.4, pars=c("alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]", "alpha[5]", "alpha[6]","alpha[7]","alpha[8]","alpha[9]",
                               "lambda", "k", "llin", "age_0", "age_1", "age_2", "age_3", "age_4"),                
          facet_args = list(strip.position = "top",nrow=5))

ggsave("model.4_posterior.png", dpi=600, width=8, height=7)

#------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------duration of protection--------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
#For model.1
dur_weib_protect(model.1)

#For mode.2
dur_weib_protect(model.2)

#For model.3
dur_weib_protect(model.3)

#For model.4
dur_weib_protect(model.4)
over_90(model.4)

#------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------posterior draws of p_protect files-----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
post_model.4 <- rstan::extract(model.4)

# time seq from old file 
t <- seq(0,60,0.2)

# extract weibull lambda and k  
output_mcmc = cbind(post_model.4[["lambda"]], post_model.4[["k"]])
par_med     = c(median(post_model.4[["lambda"]]), median(post_model.4[["k"]]))
par_med     = data.frame(lambda = par_med[1], k = par_med[2], draw=0)

params = as.data.frame(output_mcmc)

# keep every 200th row to get 50 different values to sample over 
params.new = params[seq(1, nrow(params), 200), ] %>% 
  dplyr::select(lambda = V1,
                k = V2) %>%
  dplyr::mutate(draw=seq(1,50,1))

params.new <- bind_rows(par_med, params.new)

# make p_protect txt file for all outputs for each parameter value 

# make a output matrix 
p_protect_out <- matrix(nrow=length(t), ncol=4)

# fill and save
for(i in 0:50){
  
  df <- params.new %>% filter(draw == i) 
  
  p_protect <- exp(-(t / df$lambda)^ df$k)
  
  # fill the matrix 
  p_protect_out[,1] <- t
  p_protect_out[,2] <- 0
  p_protect_out[,3] <- 200
  p_protect_out[,4] <- p_protect
  
  # write the matrix to a tab seperated txt file 
  write.table(p_protect_out, file=here::here("p_protect_files",paste0("p_protect_3_",df$draw)), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
}

#-to make the discussion plot-----------------------------------------------------------------------------------------------------
# model 4 efficacy curve - but 
source("model_plot_functions.R") #but comment back in the trail line
fig_1h <- duration_protection(model.4)
fig_1h
ggsave("discussion_plot.png", width=6, height=5, dpi=600)
