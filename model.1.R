
control_data   <- agg_data %>% filter(treatment == 0) 
treatment_data <- agg_data %>% filter(treatment == 1)  
agg_data       <- bind_rows(control_data, treatment_data) %>% arrange(desc(malaria_event))

# age groups 
agg_data$age_0  = ifelse(agg_data$age_at_enrolment == 0, 1, 0)
agg_data$age_1  = ifelse(agg_data$age_at_enrolment == 1, 1, 0)
agg_data$age_2  = ifelse(agg_data$age_at_enrolment == 2, 1, 0)
agg_data$age_3  = ifelse(agg_data$age_at_enrolment == 3, 1, 0)
agg_data$age_4  = ifelse(agg_data$age_at_enrolment == 4, 1, 0)
# site location 
agg_data$site_1 = ifelse(agg_data$site == 2, 1, 0)
agg_data$site_2 = ifelse(agg_data$site == 1, 1, 0)
agg_data$site_3 = ifelse(agg_data$site == 3, 1, 0)

# making age 
#agg_data$age_at_enrolment = ifelse(agg_data$age_at_enrolment == 0, 0.5, agg_data$age_at_enrolment)

agg_data <- agg_data %>% na.omit() %>% arrange(desc(malaria_event))

#stan data list 
stan_data <- list( N               = nrow(agg_data),
                   RT              = agg_data$event_time,
                   event_indicator = agg_data$malaria_event, 
                   treatment_id    = agg_data$treatment, 
                   Ti              = 60,                                     
                   time_seq        = seq(1,60), 
                   time_id         = agg_data$event_time,
                   bed_net         = agg_data$use_of_itn,
                   age_1           = agg_data$age_1,
                   age_2           = agg_data$age_2,
                   age_3           = agg_data$age_3,
                   age_4           = agg_data$age_4,
                   site_2          = agg_data$site_2, 
                   site_3          = agg_data$site_3) 

#stan data list plain
stan_data_plain <- list( N               = nrow(agg_data),
                   RT              = agg_data$event_time,
                   event_indicator = agg_data$malaria_event, 
                   treatment_id    = agg_data$treatment, 
                   Ti              = 60,                                     
                   time_seq        = seq(1,60), 
                   time_id         = agg_data$event_time) 

#stan data list bed nets
stan_data_llin <- list( N               = nrow(agg_data),
                        RT              = agg_data$event_time,
                        event_indicator = agg_data$malaria_event, 
                        treatment_id    = agg_data$treatment, 
                        Ti              = 60,                                     
                        time_seq        = seq(1,60), 
                        time_id         = agg_data$event_time,
                        bed_net         = agg_data$use_of_itn)

#age and llin 
stan_data_llin_age <- list( N               = nrow(agg_data),
                            RT              = agg_data$event_time,
                            event_indicator = agg_data$malaria_event, 
                            treatment_id    = agg_data$treatment, 
                            Ti              = 60,                                     
                            time_seq        = seq(1,60), 
                            time_id         = agg_data$event_time,
                            bed_net         = agg_data$use_of_itn,
                            age             = agg_data$age_at_enrolment) 

#age and llin and site
stan_data_llin_age_site <- list( N          = nrow(agg_data),
                            RT              = agg_data$event_time,
                            event_indicator = agg_data$malaria_event, 
                            treatment_id    = agg_data$treatment, 
                            Ti              = 60,                                     
                            time_seq        = seq(1,60), 
                            time_id         = agg_data$event_time,
                            bed_net         = agg_data$use_of_itn,
                            age             = agg_data$age_at_enrolment, 
                            site            = agg_data$site) 
