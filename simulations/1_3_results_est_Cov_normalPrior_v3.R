#-----------------------------
# results from 1k simulated data
#-----------------------------
# module load R_3.5.3-9.3.0 binutils_2.30
dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'
b = 0
load(paste0(dir, "/sim_data/sim_dat_3B_N_100_R_corr_", b, "_Complete_0.8.Rdata"), dat <- new.env())

nsims=dat$nsims
cov_true = dat$sim_param$cov_alpha_true
cov_params = rep(0, dim(cov_true)[1] * dim(cov_true)[2])

idx = 0
for (nrow in 1:dim(cov_true)[1]){
	for (ncol in 1:dim(cov_true)[2]){
		idx = idx + 1
		cov_params[idx] = paste('cov_alpha[', nrow, ',', ncol, ']', sep='')
	}
}

cov_mean_list = cov_q025_list = cov_q975_list = cov_cov_95 = list()
cov_coverage_95 = matrix(rep(0, nsims*length(cov_params)), nrow=nsims)
colnames(cov_coverage_95) = cov_params
error_nsim = NULL
for (i in 1:nsims){
	tmp = try(load(paste0(dir, '/results_test_cov/R_corr_', dat$sim_param$R_corr, 
                          '_test_cov_normalPrior_v3_', i, '.Rdata')), TRUE)
    if(isTRUE(class(tmp)=="try-error")) { 
        error_nsim = c(error_nsim, i);
        cov_coverage_95[i,] = NA
        next } 

	cov_mean_list[[i]] = results_cov[, 'mean']
	cov_q025_list[[i]] = results_cov[, '2.5%']
	cov_q975_list[[i]] = results_cov[, '97.5%']

	idx_cov = 0
	for (nrow in 1:dim(cov_true)[1]){
		for (ncol in 1:dim(cov_true)[2]){
			idx_cov = idx_cov + 1
			if (cov_true[nrow, ncol] <= cov_q975_list[[i]][idx_cov] & 
				cov_true[nrow, ncol] >= cov_q025_list[[i]][idx_cov]){
			cov_coverage_95[i, idx_cov] = 1
			} 
		}	
	}

}

# save results
cov_coverage_95_prob = colMeans(cov_coverage_95, na.rm=T)
# cov_est = apply(simplify2array(cov_mean_list), 1:2, mean, na.rm=T) 
# cov_q025_est = apply(simplify2array(cov_q025_list), 1:2, mean, na.rm=T)
# cov_q975_est = apply(simplify2array(cov_q975_list), 1:2, mean, na.rm=T)

save(nsims, cov_true, cov_coverage_95_prob, error_nsim,
	file=paste0('/home/lij014/Oct_2020_mfpca/simulations/results_summary/results_R_corr_', 
                b, '_test_cov_normalPrior_v3.RData'))
