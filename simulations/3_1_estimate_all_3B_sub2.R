############################################################
###### estimate all parameters ########
############################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1] # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'

library(parallel)
library(rstan)
options(mc.cores = parallel::detectCores())
Nsamples = 1000
Nchains = 3
Ncores=Nchains

source("sim_settings_v2.R") # for reformat simulated data

load(paste0(dir, '/sim_data_NB/sim_3B_file_names.Rdata'))

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(i)

for (b in 19:21){
    print(paste0('file_idx: ', b))
    file = file_names[[b]]
    load(file, dat <- new.env())
    nsims = dat$nsims
    N = dat$sim_param$N
    theta_mu_true = dat$sim_param$params[[7]]
    sigma_eps_true = dat$sim_param$params[[5]]
    Theta_true = dat$sim_param$params[[6]]
    cov_true = dat$sim_param$cov_alpha_true

    model_file = "mfpca_weakPrior.stan"
    smod = stan_model(model_file)

    print(paste0('i: ', i))
    
    error_catch = try(sim_data_reformat(dat = dat, nsim=i), TRUE) # check whether error occurs
    if (class(error_catch) == 'try-error'){
        print('NA')
        fit = dat_new = Nchains = Nsamples = diff_time = NA
    }else{
        dat_new = sim_data_reformat(dat = dat, nsim=i)
        pca_data <- list(N = dat_new$num_subjects, P = dat_new$num_blocks, K = dat_new$num_PCs, Q = dat_new$num_basis,
                 M = dat_new$num_visits_all, subject_starts = dat_new$subject_starts, subject_stops = dat_new$subject_stops,
                 subject_starts_each = dat_new$subject_starts_each, subject_stops_each = dat_new$subject_stops_each,
                 V = dat_new$visits_matrix, B = dat_new$spline_basis, Y = dat_new$response_vector, 
                 len_Theta = dat_new$len_Theta, Theta_idx = dat_new$Theta_idx, Q_idx = dat_new$Q_idx_Theta,
                 K_idx=dat_new$K_idx_Theta) 

        start_time = Sys.time()
        fit = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains, cores=Ncores)
        end_time = Sys.time()
        diff_time = end_time - start_time # can use as.numeric(diff_time) to extract values
     }
            
    save(fit, dat_new, Nchains, Nsamples, diff_time, 
         file=paste0(dir, '/results_sim_NB/', i, '_', substr(file, start=81, stop=1000)))
}    



