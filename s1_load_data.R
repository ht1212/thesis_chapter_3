
#load survival data 

agg_data <- readRDS("survival_data_processed.RDS")

# observed survival data
raw_survival <- survfit(Surv(event_time, malaria_event) ~ treatment , data=agg_data)

trial_plot <- ggsurvplot(raw_survival,
                         data = agg_data,
                         break.time.by = 7,
                         xlab="Time since SP+AQ dose (days)",
                         fun = "event",
                         censor = FALSE,
                         palette = c( "black", "black"),
                         ylab = "Proportion diagnosed",
                         legend = "none",
                         ggtheme = theme_bw(11),
                         conf.int = FALSE)

#trial_plot # check - yes as expected

# data frame to pass to ggplot for model comparison plots
t.d.p <- data.frame(time= trial_plot[["data.survplot"]][["time"]],
                    surv= 1-trial_plot[["data.survplot"]][["surv"]],
                    treatment = trial_plot[["data.survplot"]][["treatment"]])

source("trial_efficacy_data.R")
