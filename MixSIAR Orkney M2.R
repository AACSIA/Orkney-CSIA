
title: "Orkney_mixingmodel protein"
author: "AB"
date: "2025-03-13"
output: html_document

#1. Load basic packages 


library(MixSIAR)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)
library(Hmisc)
library(ggpubr)
library(googlesheets4)
library (devtools)
library(R2jags)
library(coda)
library(rjags)
library(tibble)
library(janitor)



#2. Data sheet with orkney info on


Data <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1JVz8aStoMTwlUspwgTdmeGZE9bME_AsPrMahzm2QW_w/edit?gid=810027140#gid=810027140", sheet="Mixsiar")|> as.data.frame()

#Load mixture and source data
mix <- load_mix_data(filename="/Users/annie/OneDrive/Brazil/Mixsiar/HumansM2.csv", iso_names=c("ile_d13C", "leu_d13C", "val_d13C", "phe_d13C", "phe_d15N", "lys_d15N", "glu_d15N","glu_d13C","asp_d13C","ala_d13C"),
                     factors="Sample_ID",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

source <- load_source_data(filename="/Users/annie/OneDrive/Brazil/Mixsiar/SourcesM2.csv",
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)

#Load discrimination data including the offests for the EAAs if requried. Note d13C lysine does not have a value. 
discr <- load_discr_data(filename="/Users/annie/OneDrive/Brazil/Mixsiar/Discr_offsetsM2.csv", mix)

#4b Create the JAGS file and run model

# Write the JAGS model file
model_filename <- "MixSIAR_modelM2.txt"# Name of the JAGS model file
resid_err <-FALSE
process_err <-TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#4c Check model by running a test

jags.1 <- run_model(run="test", mix, source, discr, model_filename)

# Run full model in long mode

jags.1 <- run_model(run="long", mix, source, discr, model_filename)

#5. Obtain model outputs and check for convergence

#5a. Obtain the model outputs 

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_final",
                       sup_post = TRUE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "posterior_density",
                       sup_pairs = TRUE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = TRUE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE,
                       return_obj=TRUE)

# Get posterior chains into one tidy data frame

summary(jags.1)
attach.jags(jags.1)

#Run diagnostics and create a data frame with some summary stats.Check the stats in the consule to make sure the model has converged.  

output_JAGS(jags.1, mix, source, output_options)

output_JAGS(
  jags.1,
  mix,
  source,
  output_options = list(summary_save = TRUE, summary_name = "summary_statistics",
                        sup_post = FALSE, plot_post_save_pdf = FALSE, plot_post_name = "Whole diet",
                        sup_pairs = FALSE, plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot", sup_xy
                        = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot", gelman = TRUE, heidel =
                          FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics", indiv_effect =
                          FALSE, plot_post_save_png = TRUE, plot_pairs_save_png = TRUE, plot_xy_save_png =
                          TRUE, diag_save_ggmcmc = TRUE))
