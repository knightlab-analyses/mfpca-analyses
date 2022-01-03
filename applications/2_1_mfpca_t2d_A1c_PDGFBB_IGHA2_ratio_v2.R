#-------------------------------------------------------------
# run on barnacle for model selection
# include all Omics measurements & no removal of NA's in A1C
# use less extreme pseudo number to deal with zero: 0.001
#-------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../data/mfpca_apply_fcns.R')
#--------------------------
# 1. prepare data
#--------------------------
df = read.csv('../data/HMP2_t2d_update120420.txt', sep='\t')
colnames(df)[2] = 'SubjectUnique'
df$ratio = log((df$genus_unclassified_Lachnospiraceae/df$genus_Bacteroides) + 0.001)
df$PDGFBB_log = log(df$PDGFBB)

microb = df[, c('SampleID', 'SubjectUnique', 'CollectionDate', 'ratio', 'CL4')]
microb = microb[complete.cases(microb$ratio), ]

protein = df[, c('SubjectUnique', 'CollectionDate', 'IGHA2', 'CL4')]
protein = protein[complete.cases(protein$IGHA2), ]

cytok = df[, c('SubjectUnique', 'CollectionDate', 'PDGFBB_log', 'CL4')]
cytok = cytok[complete.cases(cytok$PDGFBB_log), ]

a1c = df[, c('SubjectUnique', 'CollectionDate', 'A1C', 'CL4')]
a1c = a1c[complete.cases(a1c$A1C), ]

dat_mfpca <- data_mfpca_transformed(df_var=list(microb, protein, cytok, a1c), 
                                    shared_id_unique = 'SubjectUnique', 
									shared_time_var = 'CollectionDate', 
                                    y_var = c('ratio', 'IGHA2', 
                                              'PDGFBB_log', 'A1C'), 
									covariates='CL4', keep_all_obs=FALSE, min_visit=0, nseed=31)


#--------------------------
# 2. model selection
#--------------------------
# order of models: microb, protein, cytok, a1c
mdl_list <- list(list(num_PCs=c(2,2,4,4), nknots=c(1,2,1,1)),
                 list(num_PCs=c(2,3,4,4), nknots=c(1,2,1,1)),
                 list(num_PCs=c(2,2,4,4), nknots=c(1,2,2,1)),
                 list(num_PCs=c(2,2,3,3), nknots=c(1,2,2,1)),
                 list(num_PCs=c(2,2,2,3), nknots=c(1,2,2,1)),
                 list(num_PCs=c(2,2,2,2), nknots=c(1,2,2,1)),
                 list(num_PCs=c(2,2,3,3), nknots=c(1,1,1,1)))

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
	 file=paste0(dir, '/results/mfpca_t2d_A1c_PDGFBB_IGHA2_ratio_v2_', i, '.RData'))





