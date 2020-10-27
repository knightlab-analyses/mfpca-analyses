############################################################
###### estimate all parameters ########
############################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

library(parallel)
library(rstan)
options(mc.cores = parallel::detectCores())
Nsamples = 1000
Nchains = 3
Ncores=Nchains

source("sim_settings.R") # for reformat simulated data
load(paste0(dir, "/sim_data/sim_dat_3B_N_100_R_corr_2_Complete_0.8.Rdata"), dat <- new.env())
nsims = dat$nsims
N = dat$sim_param$N
theta_mu_true = dat$sim_param$params[[7]]
sigma_eps_true = dat$sim_param$params[[5]]
Theta_true = dat$sim_param$params[[6]]
cov_true = dat$sim_param$cov_alpha_true

model_file = "../code/mfpca_weakPrior.stan"
smod = stan_model(model_file)

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(i)

dat_new = sim_data_reformat(dat = dat, nsim=i)
pca_data <- list(N = dat_new$num_subjects, P = dat_new$num_blocks, K = dat_new$num_PCs, Q = dat_new$num_basis,
				 M = dat_new$num_visits_all, subject_starts = dat_new$subject_starts, subject_stops = dat_new$subject_stops,
				 subject_starts_each = dat_new$subject_starts_each, subject_stops_each = dat_new$subject_stops_each,
				 V = dat_new$visits_matrix, B = dat_new$spline_basis, Y = dat_new$response_vector, 
				 len_Theta = dat_new$len_Theta, Theta_idx = dat_new$Theta_idx, Q_idx = dat_new$Q_idx_Theta,
                 K_idx=dat_new$K_idx_Theta) 

fit = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains, cores=Ncores)
save(fit, dat_new, Nchains, Nsamples, file=paste0(dir, '/results_est_all/3D_est_all_R_v2_nsims', i, '.Rdata'))
