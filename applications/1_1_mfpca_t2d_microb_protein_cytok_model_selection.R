#-------------------------------------------
# run on barnacle for model selection
#-------------------------------------------
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
#metab = df[, c('SampleID', 'SubjectUnique', 'CollectionDate', 'shannon_metab', 'CL4')]
protein = df[, c('SampleID', 'SubjectUnique', 'CollectionDate', 'shannon_protein', 'CL4')]
cytok = df[, c('SampleID', 'SubjectUnique', 'CollectionDate', 'shannon_cytok', 'CL4')]

dat_mfpca <- data_mfpca_transformed(df_var=list(microb, protein, cytok),
                                    shared_id_unique = 'SubjectUnique', 
									shared_time_var = 'CollectionDate', 
                                    y_var = c('shannon_microb', 'shannon_protein', 'shannon_cytok'),
									covariates = 'CL4', keep_all_obs=FALSE, min_visit=0, nseed=31)

#--------------------------
# 2. model selection
#--------------------------
# models for microb + protein + cytok
mdl_list <- list(list(num_PCs=c(3,4,3), nknots=c(1,1,1)),
                 list(num_PCs=c(2,4,3), nknots=c(1,1,1)),
                 list(num_PCs=c(3,3,3), nknots=c(1,1,1)),
                 list(num_PCs=c(3,4,4), nknots=c(1,1,1)))
names(mdl_list) = paste0('mdl', seq(1, length(mdl_list), 1))

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(i)
param_stan <- input_mfpca_stan(dat_mfpca=dat_mfpca, num_PCs=mdl_list[[i]]$num_PCs, 
							   nknots=mdl_list[[i]]$nknots, nseed=i)

pca_data <- list(N = param_stan$num_subjects, P = param_stan$num_blocks, K = param_stan$num_PCs, Q = param_stan$num_basis,
				 cov_size = param_stan$cov_size, cov_starts = param_stan$cov_starts, cov_stops = param_stan$cov_stops,
				 M = param_stan$num_visits_all, subject_starts = param_stan$subject_starts, subject_stops = param_stan$subject_stops,
				 subject_starts_each = param_stan$subject_starts_each, subject_stops_each = param_stan$subject_stops_each,
				 V = param_stan$visits_matrix, B = param_stan$spline_basis, Y = param_stan$response_vector, 
				 len_Theta = param_stan$len_Theta, Theta_idx = param_stan$Theta_idx, Q_idx = param_stan$Q_idx_Theta, K_idx=param_stan$K_idx_Theta) 


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
fit = sampling(smod, data = pca_data, iter=Nsamples, chains=Nchains, cores=Ncores)
save(fit, mdl_list, Nsamples, Nchains, pca_data, param_stan, dat_mfpca, 
	 file=paste0(dir, '/results/mfpca_t2d_microb_protein_cytok_mdl_', i, '.RData'))




