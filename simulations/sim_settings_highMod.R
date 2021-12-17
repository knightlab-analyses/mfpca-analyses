#-------------------------------------
# simulate data with high modalities
#-------------------------------------


sim_setting_5B = function(T_min=0, T_max=1, N_T_max=10, N_T_mean=8, P=5, N=100, nknots=c(2,1,1,1,1), orth=TRUE, noise=0.01,
	                         Q_true=c(2,2,1,2,2), D_sd=c(2,1,2,1,1,1,1,1,1), R_corr=3,
	                         MU_true=c(c(-0.4203558,-0.2586524,0.3162561,0.3769931,0.2210936,0.1521339), # 16s_mean
	                         	       c(-0.48087048,0.09997186,0.14312583,0.02459493,-0.11323411), # shotgun_mean
	                       	           c(-0.19646573,-0.04440552,0.07264066,0.08270550,-0.18198343), # metab_mean
                                       c(0.15332, -0.03228, -0.00979, -0.03137, -0.04282), # HMP: 10_microbial ratio
                                       c(0.01679, -0.08319, 0.09498, -0.0296, -0.02048)), # HMP: 8_CFD
	                         THETA_true=cbind(c(0.5721651, 0.2388504, 0.4930518, 0.1558751, -0.2376850, -0.5400800, rep(0,20)), # 16s_PC1
                 							  c(0.3335056, -0.2884843, -0.6651972, -0.3244344, -0.4150059, -0.2925336, rep(0, 20)), # 16s_PC2
                 							  c(rep(0,6), 0.08243401, -0.11683765, -0.13970818, -0.35273529, 0.91411872, rep(0,15)), # shotgun_PC1
                 							  c(rep(0,6), 0.0502734,  0.4417559,  0.2677603, 0.7621657, 0.3869525, rep(0,15)), # shotgun_PC2
                 							  c(rep(0,11), 0.46726222, -0.03998188, 0.54589527, 0.60428104, 0.34192140, rep(0,10)),# metab_PC1
                                              c(rep(0,16), 0.68736, 0.3049, 0.41587, 0.09523, -0.08018, rep(0,5)), # microbial ratio PC1
                                              c(rep(0,16), 0.31981, -0.29963, -0.21997, -0.14709, -0.07309, rep(0,5)), # microbial ratio PC2
                                              c(rep(0, 21), 0.52136, 0.06733, 0.38355, -0.20622, 0.2398), # CFD PC 1
                                              c(rep(0, 21), 0.27095, 0.05542, 0.13023, -0.03162, 0.0671)), # CFD PC 2  
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
	} else if (R_corr == 5){ # for 5 modalities
		R_true[3,5] = R_true[5,3] = 0.75
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




sim_setting_10B = function(T_min=0, T_max=1, N_T_max=10, N_T_mean=8, P=10, N=100, 
                           nknots=c(2,1,1,1,1,2,2,1,1,1), orth=TRUE, noise=0.01,
	                       Q_true=c(2,2,1,rep(2,7)), D_sd=c(2,1,2, rep(1,16)), R_corr=3,
                           MU_true=c(c(-0.4203558,-0.2586524,0.3162561,0.3769931,0.2210936,0.1521339), # 16s_mean
                                   c(-0.48087048,0.09997186,0.14312583,0.02459493,-0.11323411), # shotgun_mean
                                   c(-0.19646573,-0.04440552,0.07264066,0.08270550,-0.18198343), # metab_mean
                                   c(0.15332, -0.03228, -0.00979, -0.03137, -0.04282), # HMP: 10_microbial ratio
                                   c(0.01679, -0.08319, 0.09498, -0.0296, -0.02048), # HMP: 8_CFD
                                   c(-0.10397, 0.00971, 0.07276, -0.02548, -0.0226, 0.00216), # HMP: 9_IGHA2
                                   c(-0.25823, -0.0943, 0.06014, -0.02679, 0.1175, 0.02344), # HMP: 7_VTN
                                   c(-0.03123, 0.1959, -0.06728, -0.02116, 0.01609), # HMP: IP10
                                   c(0.00927, 0.02515, -0.06815, -0.03329, 0.00172), # HMP: MIG
                                   c(-0.02319, 0.10163, -0.12023, 0.01839, -0.00969)), # HMP: TGFA
	                       THETA_true=cbind(c(0.5721651, 0.2388504, 0.4930518, 0.1558751, -0.2376850, -0.5400800, 
                                              rep(0,20), rep(0, 27)), # 16s_PC1
                 							c(0.3335056, -0.2884843, -0.6651972, -0.3244344, -0.4150059, -0.2925336, 
                                              rep(0, 20), rep(0, 27)), # 16s_PC2
                 							c(rep(0,6), 0.08243401, -0.11683765, -0.13970818, -0.35273529, 0.91411872, 
                                              rep(0,15), rep(0, 27)), # shotgun_PC1
                 							c(rep(0,6), 0.0502734,  0.4417559,  0.2677603, 0.7621657, 0.3869525, 
                                              rep(0,15), rep(0, 27)), # shotgun_PC2
                 							c(rep(0,11), 0.46726222, -0.03998188, 0.54589527, 0.60428104, 0.34192140, 
                                              rep(0,10), rep(0, 27)),# metab_PC1
                                            c(rep(0,16), 0.68736, 0.3049, 0.41587, 0.09523, -0.08018, rep(0,5), 
                                              rep(0, 27)), # microbial ratio PC1
                                            c(rep(0,16), 0.31981, -0.29963, -0.21997, -0.14709, -0.07309, 
                                              rep(0,5), rep(0, 27)), # microbial ratio PC2
                                            c(rep(0,21), 0.52136, 0.06733, 0.38355, -0.20622, 0.2398, 
                                              rep(0, 27)), # CFD PC 1
                                            c(rep(0,21), 0.27095, 0.05542, 0.13023, -0.03162, 0.0671, 
                                              rep(0, 27)), # CFD PC 2 
                                            c(rep(0,26), 0.44775, 0.06591, 0.677, 0.34474, 0.29558, 0.15061, 
                                              rep(0, 21)), # IGHA2_pc1
                                            c(rep(0,26), 0.34735, -0.79692, -0.05098, -0.14884, 0.00307, -0.06748, 
                                              rep(0, 21)), # IGHA2_pc2
                                            c(rep(0,32), 0.39727, 0.48008, 0.58251, 0.31658, 0.22702, 0.13085, 
                                              rep(0,15)), # VTN: pc 1
                                            c(rep(0,32), 0.36821, 0.33103, -0.57105, 0.22785, -0.07472, 0.10843,
                                              rep(0,15)), # VTN: pc 2
                                            c(rep(0,38), 0.45792, 0.59746, 0.43786, 0.44388, 0.19057, 
                                              rep(0, 10)), # IP10: PC1
                                            c(rep(0,38), 0.26966, 0.01481, -0.11301, -0.16712, -0.04292, 
                                              rep(0, 10)), # IP10: PC2
                                            c(rep(0,43), 0.53779, 0.59574, 0.44162, 0.33604, 0.20985, rep(0,5)), # MIG: pc1
                                            c(rep(0,43), 0.24638, 0.21735, -0.51835, -0.44129, 0.54311, rep(0,5)), # MIG: pc2
                                            c(rep(0,48), 0.48826, 0.69479, 0.37806, 0.32675, 0.13464), # IGFA: pc1
                                            c(rep(0,48), 0.22832, -0.27935, 0.21909, -0.02197, 0.06159)), # IGFA: pc2            
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
	} else if (R_corr == 5){ # for 5 modalities
		R_true[3,5] = R_true[5,3] = 0.75
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
