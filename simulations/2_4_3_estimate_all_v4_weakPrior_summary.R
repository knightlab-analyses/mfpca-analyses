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

# R_v4
results_R_v4_cov = coverage_cov_95(dir_home = dir, 
                               file_data = "/sim_data/sim_dat_3B_N_100_R_corr_4_Complete_0.8.Rdata", 
                               file_results ='/summary_est/summary_3D_est_all_R_v4_nsims', 
                               nsims=i)

results_R_v4_cor = coverage_cor_95(dir_home = dir, 
                               file_data = "/sim_data/sim_dat_3B_N_100_R_corr_4_Complete_0.8.Rdata", 
                               file_results ='/summary_est/summary_3D_est_all_R_v4_nsims', 
                               nsims=i)

results_R_v4_MI = coverage_MI_95(dir_home = dir, 
                               file_results ='/summary_est/summary_3D_est_all_R_v4_nsims', 
                               nsims=i)

save(results_R_v4_cov, results_R_v4_cor, results_R_v4_MI,
     file=paste0(dir, '/cov_prob_est/cov_prob_est_all_R_v4_totalSims', i, '.Rdata'))