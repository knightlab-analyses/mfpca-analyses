# calculate computation time for 5 modalities based on 1_4_mfpca_t2d_ALL_model_selection.R 
# single model fitting, no model selection


args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../code/mfpca_apply_fcns.R')
#--------------------------
# 1. prepare data
#--------------------------
df = read.csv('../data/HMP2_t2d.txt', sep='\t')
colnames(df)[2] = 'SubjectUnique'

microb = df[, c('SampleID', 'SubjectUnique', 'CollectionDate', 'shannon_microb', 'CL4')]
metab = df[, c('SubjectUnique', 'CollectionDate', 'shannon_metab', 'CL4')]
protein = df[, c('SubjectUnique', 'CollectionDate', 'shannon_protein', 'CL4')]
cytok = df[, c('SubjectUnique', 'CollectionDate', 'shannon_cytok', 'CL4')]
a1c = df[, c('SubjectUnique', 'CollectionDate', 'A1C', 'CL4')]

# choose 5 measurements
dat_mfpca <- data_mfpca_transformed(df_var=list(microb, metab, protein, cytok, a1c), 
                                    shared_id_unique = 'SubjectUnique', 
									shared_time_var = 'CollectionDate', 
                                    y_var = c('shannon_microb', 'shannon_metab', 'shannon_protein', 
                                              'shannon_cytok', 'A1C'), 
									covariates='CL4', keep_all_obs=FALSE, min_visit=0, nseed=31)


#--------------------------
# 2. model selection
#--------------------------
mdl_list <- list(list(num_PCs=c(2,2,2,2,2), nknots=c(1,1,1,1,1)))

names(mdl_list) = paste0('mdl', seq(1, length(mdl_list), 1))

#i <- as.numeric(Sys.getenv("PBS_ARRAYID"))
i = 1 # since only one set of model parameters
set.seed(i)
param_stan <- input_mfpca_stan(dat_mfpca=dat_mfpca, num_PCs=mdl_list[[i]]$num_PCs, 
							   nknots=mdl_list[[i]]$nknots, nseed=i)

pca_data <- list(N = param_stan$num_subjects, P = param_stan$num_blocks, K = param_stan$num_PCs, Q = param_stan$num_basis,
				 cov_size = param_stan$cov_size, cov_starts = param_stan$cov_starts, cov_stops = param_stan$cov_stops,
				 M = param_stan$num_visits_all, subject_starts = param_stan$subject_starts, subject_stops =
                 param_stan$subject_stops,
				 subject_starts_each = param_stan$subject_starts_each, subject_stops_each = param_stan$subject_stops_each,
				 V = param_stan$visits_matrix, B = param_stan$spline_basis, Y = param_stan$response_vector, 
				 len_Theta = param_stan$len_Theta, Theta_idx = param_stan$Theta_idx, Q_idx = param_stan$Q_idx_Theta,
                 K_idx=param_stan$K_idx_Theta) 


#--------------------------------------------
# MFPCA modeling
#--------------------------------------------
library(parallel)
library(rstan)
options(mc.cores = parallel::detectCores())
Nsamples = 1000
Nchains = 3
Ncores=Nchains

model_file = "../code/mfpca_weakPrior_loo.stan"
smod = stan_model(model_file)
start_time = Sys.time()
fit = sampling(smod, data = pca_data, iter=Nsamples, chains=Nchains, cores=Ncores)
end_time = Sys.time()
diff_time = end_time - start_time

save(diff_time, mdl_list, dat_mfpca,
     file=paste0(dir, '/results/mfpca_t2d_compute_time_mod5.RData'))



