
# ggplot(agg_data %>% filter(site==2)) + 
#   geom_histogram(aes(x=age_at_enrolment, y=..density..)) + 
#   facet_wrap(treatment ~ .)
library(rstan)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(shinystan)
library(bayesplot)
library(msm)
library(loo)
library(RColorBrewer)
library(compiler)
library(patchwork)

# # remove.packages("BH")
# # require(devtools)
# # install_version("BH", version = "1.69.0-1", repos = "http://cran.us.r-project.org")
# 
# remove.packages(c("StanHeaders", "rstan"))
# install.packages(c("rstan","StanHeaders"))
# #Then install the older version of BH:
# devtools::install_version("BH",version="1.72.0-3")