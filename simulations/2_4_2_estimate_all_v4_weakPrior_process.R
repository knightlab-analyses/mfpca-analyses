############################################################
###### estimate all parameters ########
############################################################
args = commandArgs(trailingOnly=TRUE)
dir = args[1]

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))


#----------------------
#check results
#----------------------
source('mfpca_rotation_v2.R')
source('sim_settings_v2.R') # for plotting figures
source('basis_setup_sparse.R')

# results_nsim = output_sim(mfpca_output=paste0(dir, '/results_est_all/3D_est_all_R_v4_nsims', i, '.Rdata'),
#                            dat_sim =  paste0(dir, "/sim_data/sim_dat_3B_N_100_R_corr_4_Complete_0.8.Rdata"),
#                            idx_sim=i, mean_lim=c(-15,15), FPC_lim=c(-5,5),
#                            fig_name=paste0(dir, '/figures_est_all/3D_est_all_R_v4_nsims', i, '.png'),
#                            fig_on=FALSE)

# plot figure using pdf()
results_nsim = output_sim(mfpca_output = paste0(dir, '/results_est_all/3D_est_all_R_v4_nsims', i, '.Rdata'),
                           dat_sim = paste0(dir, "/sim_data/sim_dat_3B_N_100_R_corr_4_Complete_0.8.Rdata"),
                           idx_sim = i, mean_lim = c(-15,15), FPC_lim = c(-5,5),
                           fig_name = paste0('figures/3D_est_all_R_v4_nsims', i, '.pdf'),
                           fig_on = TRUE)

save(results_nsim, file=paste0(dir, '/summary_est/summary_3D_est_all_R_v4_nsims', i, '.Rdata'))


