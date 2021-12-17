####################################################
### perform post-hoc rotations on MFPCA outputs ####
####################################################
library(Matrix)
library(rstan)
library(corpcor) # covariance decomposition 

post_hoc_rotation = function(fit, Nchains, Nsamples, N, P, K, Q){
	###################################
	##### parameters list ############
	# fit: model returns from mfpca
	# Nchains: #chains used in stan
	# Nsamples: #samples used in stan

	# N = dat_new$num_subjects
	# P = dat_new$num_blocks
	# K = dat_new$num_PCs
	# Q = dat_new$num_basis

	Theta = extract(fit, "Theta", permuted=FALSE)
	cov = extract(fit, "cov_alpha", permuted=FALSE)
	theta_mu = extract(fit, "theta_mu", permuted=FALSE)
	alpha = extract(fit, "alpha", permuted=FALSE)
	sigma_eps = extract(fit, "sigma_eps", permuted=FALSE)

	Theta_old = Theta_new = array(0, dim=c(sum(Q), sum(K), Nchains*Nsamples/2))
    cov_old = cov_new = corr_new = array(0, dim=c(sum(K), sum(K), Nchains*Nsamples/2))
	theta_mu_new = array(0, dim=c(sum(Q), Nchains*Nsamples/2))
	alpha_old = alpha_new = array(0, dim=c(sum(K), N, Nchains*Nsamples/2)) 

	ind = 0
	for (i in 1:dim(Theta)[1]){
		for (j in 1:dim(Theta)[2]){
			ind = ind + 1 # loop through all posterior draws
			Theta_old[,,ind] = array(Theta[i,j,],dim=c(sum(Q), sum(K)))
			cov_old[,,ind] = array(cov[i,j,], dim=c(sum(K), sum(K)))
			theta_mu_new[,ind] = array(theta_mu[i,j,])
			alpha_old[,,ind] = array(alpha[i,j,],dim=c(sum(K), N))

			indq = 1;
	    	indk = 1;
	    	for (p in 1:P){
	    	poolvar = as.matrix(cov_old[,,ind][indk:sum(K[1:p]),indk:sum(K[1:p])]) 
	    	temp_sigma = Theta_old[,,ind][indq:sum(Q[1:p]),indk:sum(K[1:p])] %*% poolvar %*% t(Theta_old[,,ind][indq:sum(Q[1:p]),indk:sum(K[1:p])])
	    	eigen_temp_sigma=eigen(temp_sigma)
			v_temp=Re(eigen_temp_sigma$vectors)
			d_temp=Re(eigen_temp_sigma$values)

			# rotate Theta
			Theta_new[,,ind][indq:sum(Q[1:p]),indk:sum(K[1:p])] = v_temp[, 1:K[p]]
			for(k in 1:K[p]){
	        	Theta_new[,,ind][,(k+indk-1)] = sign(Theta_new[,,ind][indq, (k+indk-1)]) * Theta_new[,,ind][,(k+indk-1)];
	      	}

			indk=indk+K[p]
			indq=indq+Q[p]
			} # end p loop

			# rotate cov_alpha
			cov_new[,,ind] = t(Theta_new[,,ind]) %*% Theta_old[,,ind] %*% cov_old[,,ind]%*%
			 			     t(Theta_old[,,ind]) %*% Theta_new[,,ind] 

            # obtain correlation matrix
			corr_new[,,ind] = decompose.cov(cov_new[,,ind])$r 
            
			# rotate alpha
			alpha_new[,, ind] = t(Theta_new[,,ind]) %*% Theta_old[,,ind] %*% alpha_old[,,ind] 			   
		}
	}

	ALPHA_array = alpha_new
	MU_array = theta_mu_new
	THETA_array = Theta_new
	COV_array = cov_new
    COR_array = corr_new

	return(list('ALPHA_array'=ALPHA_array, 'MU_array'=MU_array, 'THETA_array'=THETA_array, 
                'COV_array'=COV_array, 'COR_array'=COR_array))
}


#----------------------------------------------------------------------------------------------------
# general function to calculate mutual information
# file_summary_output: file directory for output from output_sim function within sim_setting_v2.R
# K: dimension of PC for each block
#----------------------------------------------------------------------------------------------------
H_est = function(sub_corr){
	# calculate entropy

	if (is.null(dim(sub_corr))){
		H_est = 1/2 * log(sub_corr)
	}else{
		H_est = 1/2 * log(det(sub_corr))
	}
	return(H_est)
}

MI_norm = function(MI){
	# normalize MI to be [0, 1]

	MI_norm = (1 - exp(-MI * 2 )) ^ (1/2) 
	return(MI_norm)
}

# update MI calculation 
MI_est = function(corr_matrix, K){
	# calculate mutual informaion

	R = corr_matrix
	P = length(K)

	# index for each block (R_out_each: R after take out each variable; R_in_each: R of keep one variable = take out two if 3 in total)
	idx_end = idx_start = rep(0, P)
	idx_each = R_out_each = R_in_each = list()
	for (p in 1:P){
		idx_end[p] = sum(K[1:p])
		idx_start[p] = idx_end[p] - K[p] + 1 
		idx_each[[p]] = seq(idx_start[p], idx_end[p], 1)
		R_out_each[[p]] = R[-idx_each[[p]], -idx_each[[p]]] 
		R_in_each[[p]] = R[idx_each[[p]], idx_each[[p]]]
	}

    # mutual information
    MI_12 = - H_est(R_out_each[[3]])
    MI_13 = - H_est(R_out_each[[2]])
    MI_23 = - H_est(R_out_each[[1]])    
    
	# conditional mutual information
	CMI_12 = H_est(R_out_each[[2]]) + H_est(R_out_each[[1]]) - H_est(R_in_each[[3]]) - H_est(R)
	CMI_13 = H_est(R_out_each[[3]]) + H_est(R_out_each[[1]]) - H_est(R_in_each[[2]]) - H_est(R)
	CMI_23 = H_est(R_out_each[[3]]) + H_est(R_out_each[[2]]) - H_est(R_in_each[[1]]) - H_est(R)

	# normalized mutual information
	norm_MI_12 = MI_norm(MI_12)
	norm_MI_13 = MI_norm(MI_13)
	norm_MI_23 = MI_norm(MI_23)
    
	norm_CMI_12 = MI_norm(CMI_12)
	norm_CMI_13 = MI_norm(CMI_13)
	norm_CMI_23 = MI_norm(CMI_23)    

	return(list('MI_12'= MI_12, 'MI_13' = MI_13, 'MI_23' = MI_23,
                'CMI_12'= CMI_12, 'CMI_13' = CMI_13, 'CMI_23' = CMI_23,
                'norm_MI_12' = norm_MI_12, 'norm_MI_13' = norm_MI_13, 'norm_MI_23' = norm_MI_23,
                'norm_CMI_12' = norm_CMI_12, 'norm_CMI_13' = norm_CMI_13, 'norm_CMI_23' = norm_CMI_23))
}

MI_est_sim = function(post_rotation_results, K, R){
	# R: true correlation matrix from simulated data
	
	R_array = post_rotation_results$COR_array
	MI_true_results = MI_est(corr_matrix=R, K=K) 

	nloop = dim(R_array)[3]
	MI_12 = MI_13 = MI_23 = norm_MI_12 = norm_MI_13 = norm_MI_23 = 
    CMI_12 = CMI_13 = CMI_23 = norm_CMI_12 = norm_CMI_13 = norm_CMI_23 = rep(0, nloop)
    for (i in 1:nloop){
    	MI_results = MI_est(corr_matrix=R_array[,,i], K=K)
		
		# store original mutual info
		MI_12[i] = MI_results$MI_12
		MI_13[i] = MI_results$MI_13
		MI_23[i] = MI_results$MI_23
		CMI_12[i] = MI_results$CMI_12
		CMI_13[i] = MI_results$CMI_13
		CMI_23[i] = MI_results$CMI_23        

		# store normalized mutual info
		norm_MI_12[i] = MI_results$norm_MI_12
		norm_MI_13[i] = MI_results$norm_MI_13
		norm_MI_23[i] = MI_results$norm_MI_23
		norm_CMI_12[i] = MI_results$norm_CMI_12
		norm_CMI_13[i] = MI_results$norm_CMI_13
		norm_CMI_23[i] = MI_results$norm_CMI_23        
	}

	MI_12_list = c(mean(MI_12), quantile(MI_12, 0.5), quantile(MI_12, 0.025), quantile(MI_12, 0.975))
	MI_13_list = c(mean(MI_13), quantile(MI_13, 0.5), quantile(MI_13, 0.025), quantile(MI_13, 0.975))
	MI_23_list = c(mean(MI_23), quantile(MI_23, 0.5), quantile(MI_23, 0.025), quantile(MI_23, 0.975))
	CMI_12_list = c(mean(CMI_12), quantile(CMI_12, 0.5), quantile(CMI_12, 0.025), quantile(CMI_12, 0.975))
	CMI_13_list = c(mean(CMI_13), quantile(CMI_13, 0.5), quantile(CMI_13, 0.025), quantile(CMI_13, 0.975))
	CMI_23_list = c(mean(CMI_23), quantile(CMI_23, 0.5), quantile(CMI_23, 0.025), quantile(CMI_23, 0.975))    

	norm_MI_12_list = c(mean(norm_MI_12), quantile(norm_MI_12, 0.5), 
						quantile(norm_MI_12, 0.025), quantile(norm_MI_12, 0.975))
	norm_MI_13_list = c(mean(norm_MI_13), quantile(norm_MI_13, 0.5), 
					    quantile(norm_MI_13, 0.025), quantile(norm_MI_13, 0.975))
	norm_MI_23_list = c(mean(norm_MI_23), quantile(norm_MI_23, 0.5), 
						quantile(norm_MI_23, 0.025), quantile(norm_MI_23, 0.975))
	norm_CMI_12_list = c(mean(norm_CMI_12), quantile(norm_CMI_12, 0.5), 
						quantile(norm_CMI_12, 0.025), quantile(norm_CMI_12, 0.975))
	norm_CMI_13_list = c(mean(norm_CMI_13), quantile(norm_CMI_13, 0.5), 
					    quantile(norm_CMI_13, 0.025), quantile(norm_CMI_13, 0.975))
	norm_CMI_23_list = c(mean(norm_CMI_23), quantile(norm_CMI_23, 0.5), 
						quantile(norm_CMI_23, 0.025), quantile(norm_CMI_23, 0.975))
    
	return(list('MI_true_results'=MI_true_results, 
                'MI_12_list'=MI_12_list, 'MI_13_list'=MI_13_list, 'MI_23_list'=MI_23_list, 
                'CMI_12_list'=CMI_12_list, 'CMI_13_list'=CMI_13_list, 'CMI_23_list'=CMI_23_list,                 
				'norm_MI_12_list'=norm_MI_12_list, 'norm_MI_13_list'=norm_MI_13_list, 'norm_MI_23_list'=norm_MI_23_list,
                'norm_CMI_12_list'=norm_CMI_12_list, 'norm_CMI_13_list'=norm_CMI_13_list, 'norm_CMI_23_list'=norm_CMI_23_list))

}

#---------------------------------------
# comment out original MI calculation
#---------------------------------------
# MI_est = function(corr_matrix, K){
# 	# calculate mutual informaion

# 	R = corr_matrix
# 	P = length(K)

# 	# index for each block (R_out_each: R after take out each variable; R_in_each: R of keep one variable = take out two if 3 in total)
# 	idx_end = idx_start = rep(0, P)
# 	idx_each = R_out_each = R_in_each = list()
# 	for (p in 1:P){
# 		idx_end[p] = sum(K[1:p])
# 		idx_start[p] = idx_end[p] - K[p] + 1 
# 		idx_each[[p]] = seq(idx_start[p], idx_end[p], 1)
# 		R_out_each[[p]] = R[-idx_each[[p]], -idx_each[[p]]] 
# 		R_in_each[[p]] = R[idx_each[[p]], idx_each[[p]]]
# 	}

# 	# mutual information
# 	MI_all = -H_est(R)
# 	MI_12 = H_est(R_out_each[[2]]) + H_est(R_out_each[[1]]) - H_est(R_in_each[[3]]) - H_est(R)
# 	MI_13 = H_est(R_out_each[[3]]) + H_est(R_out_each[[1]]) - H_est(R_in_each[[2]]) - H_est(R)
# 	MI_23 = H_est(R_out_each[[3]]) + H_est(R_out_each[[2]]) - H_est(R_in_each[[1]]) - H_est(R)

# 	# normalized mutual information
# 	norm_MI_all = MI_norm(MI_all)
# 	norm_MI_12 = MI_norm(MI_12)
# 	norm_MI_13 = MI_norm(MI_13)
# 	norm_MI_23 = MI_norm(MI_23)

# 	return(list('MI_all'=MI_all, 'MI_12'=MI_12, 'MI_13'=MI_13, 'MI_23'=MI_23,
# 				'norm_MI_all'=norm_MI_all, 'norm_MI_12'=norm_MI_12,
# 				'norm_MI_13'=norm_MI_13, 'norm_MI_23'=norm_MI_23))
# }

# MI_est_sim = function(post_rotation_results, K, R){
# 	# R: true correlation matrix from simulated data
	
# 	R_array = post_rotation_results$COR_array
# 	MI_true_results = MI_est(corr_matrix=R, K=K) 

# 	nloop = dim(R_array)[3]
# 	MI_all = MI_12 = MI_13 = MI_23 = norm_MI_all = norm_MI_12 = norm_MI_13 = norm_MI_23 = rep(0, nloop)
#     for (i in 1:nloop){
#     	MI_results = MI_est(corr_matrix=R_array[,,i], K=K)
		
# 		# store original mutual info
# 		MI_all[i] = MI_results$MI_all
# 		MI_12[i] = MI_results$MI_12
# 		MI_13[i] = MI_results$MI_13
# 		MI_23[i] = MI_results$MI_23

# 		# store normalized mutual info
# 		norm_MI_all[i] = MI_results$norm_MI_all
# 		norm_MI_12[i] = MI_results$norm_MI_12
# 		norm_MI_13[i] = MI_results$norm_MI_13
# 		norm_MI_23[i] = MI_results$norm_MI_23
# 	}

# 	MI_all_list = c(mean(MI_all), quantile(MI_all, 0.5), quantile(MI_all, 0.025), quantile(MI_all, 0.975))
# 	MI_12_list = c(mean(MI_12), quantile(MI_12, 0.5), quantile(MI_12, 0.025), quantile(MI_12, 0.975))
# 	MI_13_list = c(mean(MI_13), quantile(MI_13, 0.5), quantile(MI_13, 0.025), quantile(MI_13, 0.975))
# 	MI_23_list = c(mean(MI_23), quantile(MI_23, 0.5), quantile(MI_23, 0.025), quantile(MI_23, 0.975))

# 	norm_MI_all_list = c(mean(norm_MI_all), quantile(norm_MI_all, 0.5), 
# 						 quantile(norm_MI_all, 0.025), quantile(norm_MI_all, 0.975))
# 	norm_MI_12_list = c(mean(norm_MI_12), quantile(norm_MI_12, 0.5), 
# 						quantile(norm_MI_12, 0.025), quantile(norm_MI_12, 0.975))
# 	norm_MI_13_list = c(mean(norm_MI_13), quantile(norm_MI_13, 0.5), 
# 					    quantile(norm_MI_13, 0.025), quantile(norm_MI_13, 0.975))
# 	norm_MI_23_list = c(mean(norm_MI_23), quantile(norm_MI_23, 0.5), 
# 						quantile(norm_MI_23, 0.025), quantile(norm_MI_23, 0.975))

# 	return(list('MI_true_results'=MI_true_results, 'MI_all_list'=MI_all_list, 'MI_12_list'=MI_12_list,
# 				'MI_13_list'=MI_13_list, 'MI_23_list'=MI_23_list, 'norm_MI_all_list'=norm_MI_all_list,
# 				'norm_MI_12_list'=norm_MI_12_list, 'norm_MI_13_list'=norm_MI_13_list,
# 				'norm_MI_23_list'=norm_MI_23_list))

# }



