##########################################################
## various functions for simulation settings ##
##########################################################
# source("basis_setup_sparse.R")
library(Matrix)
library(rstan)

##########################################################################
#### simulated data based on DIABIMMUNE Type I diabetes ##################
##########################################################################
## real values of nknots, Q_true, MU_true, THETA_true (orthornormal) come from 16s_shannon, shotgun_shannonEquilibility, metab_shannon
sim_setting_3B_T1D = function(T_min=0, T_max=1, N_T_max=10, N_T_mean=8, P=3, N=100, nknots=c(2,1,1), orth=TRUE, noise=0.01,
	                         Q_true=c(2,2,1), D_sd=c(2, 1, 2, 1, 1), R_corr=0,
	                         MU_true=c(c(-0.4203558,-0.2586524,0.3162561,0.3769931,0.2210936,0.1521339), # 16s_mean
	                         	       c(-0.48087048,0.09997186,0.14312583,0.02459493,-0.11323411), # shotgun_mean
	                       	           c(-0.19646573,-0.04440552,0.07264066,0.08270550,-0.18198343)), # metab_mean
	                         THETA_true=cbind(c(0.5721651, 0.2388504, 0.4930518, 0.1558751, -0.2376850, -0.5400800, rep(0,10)), # 16s_PC1
                 							  c(0.3335056, -0.2884843, -0.6651972, -0.3244344, -0.4150059, -0.2925336, rep(0, 10)), # 16s_PC2
                 							  c(rep(0,6), 0.08243401, -0.11683765, -0.13970818, -0.35273529, 0.91411872, rep(0,5)), # shotgun_PC1
                 							  c(rep(0,6), 0.0502734,  0.4417559,  0.2677603, 0.7621657, 0.3869525, rep(0,5)), # shotgun_PC2
                 							  c(rep(0,11), 0.46726222, -0.03998188, 0.54589527, 0.60428104, 0.34192140)), # metab_PC1
			                 plotit=FALSE, fig_name='FPC_true.tiff'){ 
	################## description of parameters ######
	## P: number of blocks
	## N: number of subjects
	## Q_true: number of PCs in each block
	## R_corr: number of correlations across PCs
	## D_sd: standard devations of FPC scores
	######################################################

	## set up covariance matrix Sigma = D * R * D
	R_true=diag(rep(1,sum(Q_true)))
	if (R_corr == 1){
		R_true[3,5] = R_true[5,3] = 0.75
	} else if (R_corr == 2){
		R_true[3,5] = R_true[5,3] = 0.75
		R_true[2,3] = R_true[3,2] = 0.5
	} else if (R_corr == 3){
		R_true[3,5] = R_true[5,3] = 0.75
		R_true[2,3] = R_true[3,2] = 0.5
		R_true[1,5] = R_true[5,1] = 0.25
	} else if (R_corr == 4){
		R_true[3,5] = R_true[5,3] = 0.8
		R_true[2,3] = R_true[3,2] = 0.5
		R_true[1,5] = R_true[5,1] = 0.25
	}

	D_true=diag(D_sd)

	SIGMA_ALPHA_true = D_true%*%R_true%*%D_true

	SIGMA_OMEGA_true = noise*rep(1,P)

	## set up FPC_true
	T_range=c(T_min,T_max)
	time=seq(T_range[1],T_range[2],(T_range[2]-T_range[1])/(N_T_max-1))
	time_sparse=list()
	for(i in 1:N){
		time_sparse[[i]]=list()
		for(p in 1:P){
			time_sparse[[i]][[p]]=time
		}
	}

	basis_stuff=basis_setup_sparse(time_sparse,nknots,orth=orth,plotit=FALSE)
	# index below is fixed
	knots=basis_stuff[[1]]
	time_sparse=basis_stuff[[2]]
	time_sparse_combined=basis_stuff[[6]]
	phi_t=basis_stuff[[3]]
	time_cont=basis_stuff[[4]]
	phi_t_cont=	basis_stuff[[5]]
	Phi_D=basis_stuff[[7]]
	Phi_D_T=basis_stuff[[8]]

	# Set up FPCA loadings
	FPC_true=array(0,dim=c(N_T_max,sum(Q_true)))
	for(t in 1:N_T_max){
		ind=0
		for(q in 1:length(Q_true)){
			FPC_true[t,(ind+1):(ind+Q_true[q])]=(Phi_D[[1]][[t]]%*%THETA_true)[q,(ind+1):(ind+Q_true[q])]
			ind=ind+Q_true[q]
		}	
	}

	# plot FPC functions
	if (plotit == TRUE){
		tiff(fig_name, units="in", width=7, height=5, res=300)
		par(mfrow=c(2,2))
		plot(time, FPC_true[,1], type='l', lwd=2, ylim=c(min(FPC_true), max(FPC_true)))
		lines(time, FPC_true[,2], type='l', lwd=1, col='red')
		title(main='FPCs for Variable 1')
		plot(time, FPC_true[,3], type='l', lwd=2, ylim=c(min(FPC_true), max(FPC_true)))
		lines(time, FPC_true[,4], type='l', lwd=1, col='red')
		title(main='FPCs for Variable 2')
		plot(time, FPC_true[,5], type='l', lwd=2, ylim=c(min(FPC_true), max(FPC_true)))
		title(main='FPCs for Variable 3')
		dev.off()
	}


	## Set up param list
	params=list()
	params[[1]]=nknots
	params[[2]]=Q_true
	params[[3]]=R_true
	params[[4]]=D_true
	params[[5]]=SIGMA_OMEGA_true # error variance (sigma_eps in stan is SD)
	params[[6]]=THETA_true
	params[[7]]=MU_true # theta_mu
	params[[8]]=FPC_true

	return(list('T_range'=T_range, 'N_T_max'=N_T_max, 'N_T_mean'=N_T_mean, 'time'=time, 'complete'=N_T_mean/N_T_max, 
		        'P'=P, 'N'=N, 'R_corr'=R_corr, 'D_true'=D_true, 
		        'cov_alpha_true' = SIGMA_ALPHA_true,
		        'params'=params, 'orth'=orth, 'Q'=Q_true, 'nknots'=nknots))
}

sim_data_reformat = function(dat, nsim){
    ################## description of parameters ######
	## dat: input simulated data
	## nsim: index the number of repition in simulations
	######################################################

	N = dat$sim_param$N
	P = dat$sim_param$P
	K = dat$sim_param$Q
	Q = (dat$sim_param$nknots + 4)

	## below has to be generated in R; due to int/real type issue in Stan
	len_Theta = sum(K * Q) 
	Theta_idx = cumsum(K*Q)
	Q_idx_Theta = cumsum(Q)
	K_idx_Theta = cumsum(K)

	phi_t_stacked = NULL # basis
	V = matrix(0, N, P) # record the total number of visits for each subject for pth block
	visits_stop = cov_stop = rep(0, N) # stopping index for each subject in visits and covariance matrix
	time = Y = list()

	for (n in 1:N){
		Phi_D_t_n = t(dat$phi_t[[nsim]][[1]][[n]]) # transposed basis for 1st block
		for (p in 2:P){
			Phi_D_t_n = as.matrix(Matrix::bdiag(Phi_D_t_n, t(dat$phi_t[[nsim]][[p]][[n]]))) # sparse matrix (v1+v2)*(Q1+Q2)
		} # add as.matrix() to fill in zeros so that Stan can read in
	    phi_t_stacked = rbind(phi_t_stacked, Phi_D_t_n)	

		# create matrix to record number of visits
		time_n = dat$TIME_SPARSE[[nsim]][[n]]
		time[[n]] = unlist(time_n)

		for (p in 1:P){
			V[n, p] = length(time_n[[p]]) 
		}

		# stopping index
		visits_stop[n] = sum(V) # for each subject across blocks
		cov_stop[n] = sum(rowSums(V)^2) 
		visits_stop_each = cumsum(as.vector(t(V))) # for each subject in each block

		Y[[n]] = unlist(dat$Y_SPARSE[[nsim]][[n]])
	}

	visits_start = c(1, visits_stop[-N] + 1) 
	cov_start = c(1, cov_stop[-N] + 1) 
	cov_size = cov_stop[N]
	visits_start_each = c(1, visits_stop_each[-length(visits_stop_each)] + 1)

	M = sum(V) # M is the total number of records 
	#W = rowSums(V) # total number of visit for each subject combining all blocks
	#V1 = V[, 1] # visit vector for block 1
	#V2 = V[, 2] # visit vector for block 2

	return(list('num_subjects'=N, 'num_blocks'=P, 'num_PCs'=K, 'num_basis'=Q, 'response_vector'=unlist(Y), 'num_visits_all'=M,
		        'subject_starts' = visits_start, 'subject_stops' = visits_stop, 'cov_starts' = cov_start, 
                'cov_stops' = cov_stop, 
		        'subject_starts_each' = visits_start_each, 'subject_stops_each' = visits_stop_each,
		        'visits_matrix' = V, 'spline_basis' = phi_t_stacked, 'len_Theta' = len_Theta, 
		        'Theta_idx' = Theta_idx, 'Q_idx_Theta' = Q_idx_Theta, 'K_idx_Theta' = K_idx_Theta))
		        #'len_R' = len_R, 'idx_m_R' = idx_m_R, 'idx_n_R' = idx_n_R))
	# 'visits_block1' = V1, 'visits_block2' = V2, 'num_visits_each'=W, 
}

plot_sim_results = function(dat_sim, idx_sim, post_rotation_results,
							fig_name, mean_lim, FPC_lim){
	############ parameter description ######
	# dat_sim: simulated data with 100 replicates
	# idx_sim: index for replicated simulation
	# post_rotation_results: output after posthoc rotation
	# fig_name: desired figure path and name
	# mean_lim: range of y-axis for mean curve
	# FPC_lim: range of y-axis for FPC curve
	###########################################

	dat = dat_sim
	i = idx_sim
	Y_sparse = dat$Y_SPARSE[[i]]
	time_sparse = dat$TIME_SPARSE[[i]]
	results_basis=basis_setup_sparse(dat$TIME_SPARSE[[i]], nknots=dat$sim_param$nknots, orth=TRUE, plotit=FALSE)
	phi_t_cont = results_basis[[5]]
	phi_t = results_basis[[3]]
	time_cont = results_basis[[4]]
	ALPHA_array = post_rotation_results$ALPHA_array
	MU_array = post_rotation_results$MU_array
	THETA_array = post_rotation_results$THETA_array
	R_array = post_rotation_results$R_array 
	D_array = post_rotation_results$D_array

	nloop=dim(ALPHA_array)[3]
	first=1
	last=nloop

	MU_mean = MU_array[, first] #mean function across sampling sessions
	ALPHA_mean = ALPHA_array[,,first] # mean factor scores
	THETA_mean = THETA_array[,,first] # mean factor loading
	R_mean = R_array[,,first]
	D_mean = D_array[,first]

	for(iter in 2:nloop){
		MU_mean = MU_mean + MU_array[, iter]
		ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
		THETA_mean = THETA_mean + THETA_array[,,iter]
		R_mean = R_mean + R_array[,,iter]
		D_mean = D_mean + D_array[,iter]
	}

	MU_mean=cbind(MU_mean/(last-first+1)) 
	ALPHA_mean=cbind(ALPHA_mean/(last-first+1))
	THETA_mean=cbind(THETA_mean/(last-first+1)) 
	R_mean=cbind(R_mean/(last-first+1)) 
	D_mean=cbind(D_mean/(last-first+1)) 

	MU_true = dat$sim_param$params[[7]]
	Mu_true_functions=t(bdiag(cbind(phi_t_cont[[1]]),cbind(phi_t_cont[[2]]), cbind(phi_t_cont[[3]])))%*%MU_true
	Mu_functions=t(bdiag(cbind(phi_t_cont[[1]]),cbind(phi_t_cont[[2]]), cbind(phi_t_cont[[3]])))%*%MU_mean
	Mu1_true=Mu_true_functions[1:length(time_cont)]
	Mu2_true=Mu_true_functions[(length(time_cont)+1):(2*length(time_cont))]
	Mu3_true=Mu_true_functions[(2*length(time_cont)+1):(3*length(time_cont))]
	Mu1=Mu_functions[1:length(time_cont)]
	Mu2=Mu_functions[(length(time_cont)+1):(2*length(time_cont))]
	Mu3=Mu_functions[(2*length(time_cont)+1):(3*length(time_cont))]

	FPC_true = dat$sim_param$params[[8]]

	K = dat$sim_param$Q
	Q = (dat$sim_param$nknots + 4)
	N = dat$sim_param$N
	FPC1_mean=t(phi_t_cont[[1]])%*%THETA_mean[1:Q[1],1:K[1]]
	FPC2_mean=t(phi_t_cont[[2]])%*%THETA_mean[(Q[1]+1):(Q[1]+Q[2]),(K[1]+1):(K[1]+K[2])]
	FPC3_mean=t(phi_t_cont[[3]])%*%THETA_mean[(Q[1]+Q[2]+1):sum(Q),sum(K)]
	time = dat$sim_param$time

	tiff(fig_name, units="in", width=7, height=5, res=300)
	### figure 1: mean curve on individual spagetti plots
	par(mfrow=c(2,3))
	plot(time_cont,Mu1_true,type="l",ylim=mean_lim,xlim=c(0,1),lwd=1.5,col=4, ylab='response',
		font.lab=2, cex.lab=1, font.axis=1, xlab='time')
	for(i in 1:N){
		lines(time_sparse[[i]][[1]],Y_sparse[[i]][[1]],type="l",lwd=.25)
	}
	lines(time_cont,Mu1,type="l",col=2,lwd=1.5)
	title(main="Outcome for Variable 1")

	plot(time_cont,Mu2_true,type="l",ylim=mean_lim,xlim=c(0,1),lwd=1.5,col=4, ylab='response',
		font.lab=2, cex.lab=1, font.axis=1, xlab='time')
	for(i in 1:N){
		lines(time_sparse[[i]][[2]],Y_sparse[[i]][[2]],type="l",lwd=.25)
	}
	lines(time_cont,Mu2,type="l",col=2,lwd=1.5)
	title(main="Outcome for Variable 2")
	legend('bottomleft', c('True mean', 'Estimated'), col=c(4,2), lty=c(1,1), lwd=c(2,2), bty='n')

	plot(time_cont,Mu3_true,type="l",ylim=mean_lim,xlim=c(0,1),lwd=1.5,col=4, ylab='response',
		font.lab=2, cex.lab=1, font.axis=1, xlab='time')
	for(i in 1:N){
		lines(time_sparse[[i]][[3]],Y_sparse[[i]][[3]],type="l",lwd=.25)
	}
	lines(time_cont,Mu3,type="l",col=2,lwd=1.5)
	title(main="Outcome for Variable 3")	

	plot(time,FPC_true[,1],type="l",lwd=2,ylim=FPC_lim, ylab='FPC curve',
		font.lab=2, cex.lab=1, font.axis=1)
	lines(time,FPC_true[,2],type="l",lwd=1)
	lines(time_cont,FPC1_mean[,1],type="l",lwd=2,col=2)
	lines(time_cont,FPC1_mean[,2],type="l",lwd=1,col=2)
	title(main="PFCs for Variable 1")

	plot(time,FPC_true[,3],type="l",lwd=2,ylim=FPC_lim, ylab='FPC curve',
		font.lab=2, cex.lab=1, font.axis=1)
	lines(time,FPC_true[,4],type="l",lwd=1)
	lines(time_cont,FPC2_mean[,1],type="l",lwd=2,col=2)
	lines(time_cont,FPC2_mean[,2],type="l",lwd=1,col=2)
	title(main="PFCs for Variable 2")
	legend('bottomleft', c('True PC', 'Estimated'), col=c(1,2), lty=c(1,1), lwd=c(2,2), bty='n')

	plot(time,FPC_true[,5],type="l",lwd=2,ylim=FPC_lim, ylab='FPC curve',
		font.lab=2, cex.lab=1, font.axis=1)
	lines(time_cont,FPC3_mean[,1],type="l",lwd=2,col=2)
	title(main="PFCs for Variable 3")

	dev.off()

	# compute 95% CI for estimating R
	R_median = R_q025 = R_q975 = array(0,dim=c(dim(R_array)[1],dim(R_array)[2]))
	for (i in 1:dim(R_array)[1]){
	    for (j in 1:dim(R_array)[2]){
	    	R_median[i, j] = quantile(R_array[i, j, ], 0.5)
	        R_q025[i, j] = quantile(R_array[i, j, ], 0.025)
	        R_q975[i, j] = quantile(R_array[i, j, ], 0.975)
	    }

	}

	# 95% CI for Theta
	THETA_median = THETA_q025 = THETA_q975 = array(0,dim=c(dim(THETA_array)[1],dim(THETA_array)[2]))
	for (i in 1:dim(THETA_array)[1]){
	    for (j in 1:dim(THETA_array)[2]){
	    	THETA_median[i, j] = quantile(THETA_array[i, j, ], 0.5)
	        THETA_q025[i, j] = quantile(THETA_array[i, j, ], 0.025)
	        THETA_q975[i, j] = quantile(THETA_array[i, j, ], 0.975)
	    }

	}	

	# 95% CI for D
	D_median = D_q025 = D_q975 = array(0,dim=dim(D_array)[1])
	for (i in 1:dim(D_array)[1]){
    	D_median[i] = quantile(D_array[i, ], 0.5)
        D_q025[i] = quantile(D_array[i, ], 0.025)
        D_q975[i] = quantile(D_array[i, ], 0.975)
	}

	# 95% CI for theta_mu
	MU_median = MU_q025 = MU_q975 = array(0,dim=dim(MU_array)[1])
	for (i in 1:dim(MU_array)[1]){
    	MU_median[i] = quantile(MU_array[i, ], 0.5)
        MU_q025[i] = quantile(MU_array[i, ], 0.025)
        MU_q975[i] = quantile(MU_array[i, ], 0.975)
	}        

	return(list('ALPHA_mean'=ALPHA_mean, 'theta_mu_mean'=MU_mean, 'THETA_mean'=THETA_mean, 'D_mean'=D_mean, 'R_mean'=R_mean, 
		        'theta_mu_median'=MU_median, 'theta_mu_q025'=MU_q025, 'theta_mu_q975'=MU_q975,
		        'THETA_median'=THETA_median, 'THETA_q025'=THETA_q025, 'THETA_q975'=THETA_q975,
				'D_median'=D_median, 'D_q025'=D_q025, 'D_q975'=D_q975,
		        'R_median'=R_median, 'R_q025'=R_q025, 'R_q975'=R_q975))
}

output_sim = function(mfpca_output, idx_sim, fig_name, mean_lim, FPC_lim){
	load(mfpca_output)
	post_rotation_results = post_hoc_rotation(fit=fit, Nchains=Nchains, Nsamples=Nsamples, N=dat$sim_param$N, 
		                                      P=dat$sim_param$P, K=dat$sim_param$Q, Q=(dat$sim_param$nknots + 4))
	results = plot_sim_results(dat_sim=dat, idx_sim=1, post_rotation_results=post_rotation_results, 
				 	 mean_lim=mean_lim, FPC_lim=FPC_lim, fig_name=fig_name)
	return(list('results'=results, 'R_true'=dat$sim_param$params[[3]]))
}


