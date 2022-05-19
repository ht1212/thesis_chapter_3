# data subsets
control_data   <- agg_data %>% filter(treatment == 0) 
treatment_data <- agg_data %>% filter(treatment == 1)  
agg_data       <- bind_rows(control_data, treatment_data) %>% arrange(desc(malaria_event))

# making age 
#agg_data$age_at_enrolment = ifelse(agg_data$age_at_enrolment == 0, 0.5, agg_data$age_at_enrolment)

# using survSplit to create interval based data set for piecewise exponential 
breaks <- seq(30,60,30) # weekly base hazard 

split_data <- survSplit(Surv(event_time, malaria_event) ~ ., data = agg_data, cut = breaks, episode = "interval", start = "start")

split_data_2 <- mutate(split_data, exposure = event_time - start)#, interval_2 = factor(interval,  labels = paste("(", c(0,breaks), ",", c(breaks,60), "]", sep=""))) 

# tidy 
split_data_3 <- 
  split_data_2 %>%  
  dplyr::select(id, 
                interval, 
                time_in_interval = exposure, 
                malaria_event,
                event_time, 
                treatment, 
                use_of_itn,
                age_at_enrolment,
                site)  %>% 
  na.omit() %>% arrange(desc(malaria_event))

# augment the co-variates to correct format for stan model 
control_data        = split_data_3 %>% filter(treatment == 0) 
treatment_data      = split_data_3 %>% filter(treatment == 1) %>% mutate(time_in_interval = event_time)
split_data_3        = bind_rows(control_data, treatment_data) %>% arrange(desc(malaria_event))
# age groups 
# split_data_3$age_0  = ifelse(split_data_3$age_at_enrolment == 0, 1, 0)
# split_data_3$age_1  = ifelse(split_data_3$age_at_enrolment == 1, 1, 0)
# split_data_3$age_2  = ifelse(split_data_3$age_at_enrolment == 2, 1, 0)
# split_data_3$age_3  = ifelse(split_data_3$age_at_enrolment == 3, 1, 0)
# split_data_3$age_4  = ifelse(split_data_3$age_at_enrolment == 4, 1, 0)
# # site location 
# split_data_3$site_1 = ifelse(split_data_3$site == 1, 1, 0)
# split_data_3$site_2 = ifelse(split_data_3$site == 2, 1, 0)
# split_data_3$site_3 = ifelse(split_data_3$site == 3, 1, 0)

#stan data list 
# stan_data <- list( N              = nrow(split_data_3),
#                    RT             = split_data_3$time_in_interval,
#                    indicator      = split_data_3$malaria_event, 
#                    A              = length(breaks)+1,   
#                    treatment_id   = split_data_3$treatment, 
#                    interval       = split_data_3$interval,
#                    bed_net        = split_data_3$use_of_itn,
#                    age_1          = split_data_3$age_1,
#                    age_2          = split_data_3$age_2,
#                    age_3          = split_data_3$age_3,
#                    age_4          = split_data_3$age_4,
#                    site_2         = split_data_3$site_2, 
#                    site_3         = split_data_3$site_3, 
#                    Ti             = 60,                                     
#                    time_seq       = seq(1,60), 
#                    time_id        = split_data_3$time_in_interval,
#                    break_ids      = breaks) 
# 
# stan_data_plain <- list( N              = nrow(split_data_3),
#                          RT             = split_data_3$time_in_interval,
#                          indicator      = split_data_3$malaria_event, 
#                          A              = length(breaks)+1,   
#                          treatment_id   = split_data_3$treatment, 
#                          interval       = split_data_3$interval,
#                          # bed_net        = split_data_3$use_of_itn,
#                          # age_1          = split_data_3$age_1,
#                          # age_2          = split_data_3$age_2,
#                          # age_3          = split_data_3$age_3,
#                          # age_4          = split_data_3$age_4,
#                          # site_2         = split_data_3$site_2, 
#                          # site_3         = split_data_3$site_3, 
#                          Ti             = 60,                                     
#                          time_seq       = seq(1,60), 
#                          time_id        = split_data_3$time_in_interval,
#                          break_ids      = breaks) 
# 
# 
# stan_data_llin <- list( N               = nrow(split_data_3),
#                          RT             = split_data_3$time_in_interval,
#                          indicator      = split_data_3$malaria_event, 
#                          A              = length(breaks)+1,   
#                          treatment_id   = split_data_3$treatment, 
#                          interval       = split_data_3$interval,
#                          bed_net        = split_data_3$use_of_itn,
#                          # age_1          = split_data_3$age_1,
#                          # age_2          = split_data_3$age_2,
#                          # age_3          = split_data_3$age_3,
#                          # age_4          = split_data_3$age_4,
#                          # site_2         = split_data_3$site_2, 
#                          # site_3         = split_data_3$site_3, 
#                          Ti             = 60,                                     
#                          time_seq       = seq(1,60), 
#                          time_id        = split_data_3$time_in_interval,
#                          break_ids      = breaks) 
# 
#age and llin
stan_data_llin_age <- list( N               = nrow(split_data_3),
                            RT             = split_data_3$time_in_interval,
                            indicator      = split_data_3$malaria_event,
                            A              = length(breaks),
                            treatment_id   = split_data_3$treatment,
                            interval       = split_data_3$interval,
                            bed_net        = split_data_3$use_of_itn,
                            age            = split_data_3$age_at_enrolment,
                            Ti             = 60,
                            time_seq       = seq(1,60),
                            time_id        = split_data_3$time_in_interval,
                            break_ids      = breaks)

# age, llin and site
# stan_data_llin_age_site <- list( N         = nrow(split_data_3),
#                             RT             = split_data_3$time_in_interval,
#                             indicator      = split_data_3$malaria_event, 
#                             A              = length(breaks),   
#                             treatment_id   = split_data_3$treatment, 
#                             interval       = split_data_3$interval,
#                             bed_net        = split_data_3$use_of_itn,
#                             age            = split_data_3$age_at_enrolment,
#                             site           = split_data_3$site, 
#                             Ti             = 60,                                     
#                             time_seq       = seq(1,60), 
#                             time_id        = split_data_3$time_in_interval,
#                             break_ids      = breaks) 


