############################################################
###### estimate all parameters ########
############################################################
args = commandArgs(trailingOnly=TRUE)
dir = args[1]

i <- as.numeric(Sys.getenv("PBS_ARRAYID")) # total number of simulations


#----------------------
#check results
#----------------------
source('mfpca_sim_summary.R')

# R_v3
results_R_v3_cov = coverage_cov_95(dir_home = dir, 
                               file_data = "/sim_data/sim_dat_3B_N_100_R_corr_3_Complete_0.8.Rdata", 
                               file_results ='/summary_est/summary_3D_est_all_R_v3_nsims', 
                               nsims=i)

results_R_v3_cor = coverage_cor_95(dir_home = dir, 
                               file_data = "/sim_data/sim_dat_3B_N_100_R_corr_3_Complete_0.8.Rdata", 
                               file_results ='/summary_est/summary_3D_est_all_R_v3_nsims', 
                               nsims=i)

results_R_v3_MI = coverage_MI_95(dir_home = dir, 
                               file_results ='/summary_est/summary_3D_est_all_R_v3_nsims', 
                               nsims=i)

save(results_R_v3_cov, results_R_v3_cor, results_R_v3_MI,
     file=paste0(dir, '/cov_prob_est/cov_prob_est_all_R_v3_totalSims', i, '.Rdata'))