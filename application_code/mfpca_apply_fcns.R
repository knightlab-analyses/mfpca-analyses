#--------------------------------------------------------------------
# functions for MFPCA applications to microbiome multi-omics data
#--------------------------------------------------------------------


#---------------------------------------------------------------------
# transform data for MFPCA
# df_var: dataframe variable for different datasets
# shared_id_unique: same unique id identifier in two datasets
# shared_time_var: same time variable in two datasets
# y_var: response variable in each dataset (order matters!)
# covariates: interested covariates for FPC score analysis
# keep_all_obs: whether keep all obs from each dataset
# min_visit: minimum visit required for each subject
# output dataset: time scaled in [0,1]; response standardized
#----------------------------------------------------------------------

# update function to include all observations (within minimum requirement of at least one obs/block for each subject)
# keep_all_obs: if FALSE, meaning only observations (subjects + time) with measurements at all blocks are retained
data_mfpca_transformed <- function(df_var, shared_id_unique, shared_time_var, y_var, 
                                   covariates, keep_all_obs = TRUE, min_visit=0, nseed=31){

	set.seed(nseed)
	dat_list = df_var
	P = length(df_var) 
	for (i in 1:P){
		tmp = dat_list[[i]]
		dat_list[[i]]$ID_merge = paste(tmp[, shared_id_unique], tmp[, shared_time_var])	
	}

    # merge multiple data frame; keep all obs
    dat = Reduce(function(...) merge(..., 
                 by=c('ID_merge', shared_id_unique, shared_time_var, covariates), 
                 all = keep_all_obs), dat_list) 
	dat$ID = as.numeric(as.factor(dat[, shared_id_unique]))
	P_id.old=unique(dat$ID); N=length(P_id.old)
	dat$P_id.old=dat$ID
	dat$time=dat[, shared_time_var]
	dat=dat[order(dat$ID,dat$time),]

	# create visit and n_visits variable
	dat$visit = 1
	dat$n_visits = 1
	for(i in 1:N){
		dat_i=dat[dat$P_id.old==P_id.old[i],]
		dat_i$ID=i
		dat_i$visit = 1:dim(dat_i)[1]
		dat_i$n_visits = dim(dat_i)[1]	
		dat[dat$P_id.old==P_id.old[i],]=dat_i
		rm(dat_i)
	}

	# delete subjects with fewer than k visits
	k = min_visit
	dat = dat[dat$n_visits >= k, ]

	P_id.old = unique(dat$ID); N = length(P_id.old)
	dat$P_id.old=dat$ID; 
	dat$ID=0
	for(i in 1:N){
		dat_i=dat[dat$P_id.old==P_id.old[i],]
		dat_i$ID=i
		dat[dat$P_id.old==P_id.old[i],]=dat_i
		rm(dat_i)
	}

	# re-scale time to be in [0,1]
	dat$time=(dat$time-min(dat$time))/(max(dat$time-min(dat$time)))
	N=max(dat$ID)

	# standardized response variables
	dat$y = scale(dat[, y_var])

	# formatting data
	Y_sparse=list()
	time_sparse=list()
	for(i in 1:N){
		Y_sparse[[i]]=list()
		time_sparse[[i]]=list()
		for(p in 1:P){
		  Y_sparse[[i]][[p]]=dat$y[dat$ID==i & !is.na(dat$y[,p]),p]
		  time_sparse[[i]][[p]]=dat$time[dat$ID==i & !is.na(dat$y[,p])]
		}
	}
    
    # remove subjects with no observation across all time for each block             
    idx_subject_remove = NULL
    for (i in 1:length(Y_sparse)){
        if (sum(sapply(Y_sparse[[i]], length) == 0) >= 1){
            idx_subject_remove = c(idx_subject_remove, i)
        }
    }             
    if (length(idx_subject_remove) > 0){
        Y_sparse = Y_sparse[-idx_subject_remove] 
        time_sparse = time_sparse[-idx_subject_remove] 
        dat = dat[dat[,shared_id_unique] %in% unique(dat[,shared_id_unique])[idx_subject_remove], ]
    }                     
	N = length(Y_sparse)

	# record original scale (on filtered data)
	time_origin = dat[, shared_time_var]
	mu_y_origin = colMeans(dat[, y_var], na.rm=T)
	sigma_y_origin = apply(dat[, y_var], 2, sd, na.rm=T)

	return(list('Y_sparse' = Y_sparse, 'time_sparse' = time_sparse, 'N' = N, 'P' = P, 
				'data' = dat, 'time_origin' = time_origin, 'mu_y_origin'= mu_y_origin, 
				'sigma_y_origin'= sigma_y_origin, 'shared_id_unique' = shared_id_unique, 
                'covariates'= covariates, 'idx_subject_remove' = idx_subject_remove))

}
                

#-----------------------------------------------------------------------------------------------------------------
# prepare input parameters for stan
# dat_mfpca: transformed data from data_mfpca_transformed() function
#-----------------------------------------------------------------------------------------------------------------
input_mfpca_stan = function(dat_mfpca, num_PCs, nknots, nseed=31){
	set.seed(nseed)

	N = dat_mfpca$N
	P = dat_mfpca$P
	K = num_PCs
	Q = (nknots + 4)
	time_sparse = dat_mfpca$time_sparse
	Y_sparse = dat_mfpca$Y_sparse


	# below has to be generated in R; due to int/real type issue in Stan
	len_Theta = sum(K * Q) 
	Theta_idx = cumsum(K*Q)
	Q_idx_Theta = cumsum(Q)
	K_idx_Theta = cumsum(K)  	

	# set up basis
	basis_stuff = basis_setup_sparse(time_sparse, nknots)
	knots = basis_stuff[[1]]
	phi_t = basis_stuff[[3]]
	time_cont = basis_stuff[[4]]
	phi_t_cont = basis_stuff[[5]]


	phi_t_stacked = NULL # basis
	V = matrix(0, N, P) # record the total number of visits for each subject for pth block
	visits_stop = cov_stop = rep(0, N) # stopping index for each subject in visits and covariance matrix
	time = Y = list()

	for (n in 1:N){
		Phi_D_t_n = t(phi_t[[1]][[n]]) # transposed basis for 1st block
		for (p in 2:P){
			Phi_D_t_n = as.matrix(Matrix::bdiag(Phi_D_t_n, t(phi_t[[p]][[n]]))) # sparse matrix (v1+v2)*(Q1+Q2)
		} # add as.matrix() to fill in zeros so that Stan can read in
	    phi_t_stacked = rbind(phi_t_stacked, Phi_D_t_n)	

		# create matrix to record number of visits
		time_n = time_sparse[[n]]
		time[[n]] = unlist(time_n)

		for (p in 1:P){
			V[n, p] = length(time_n[[p]]) 
		}

		# stopping index
		visits_stop[n] = sum(V) # for each subject across blocks
		cov_stop[n] = sum(rowSums(V)^2) 
		visits_stop_each = cumsum(as.vector(t(V))) # for each subject in each block

		Y[[n]] = unlist(Y_sparse[[n]])
	}

	visits_start = c(1, visits_stop[-N] + 1) 
	cov_start = c(1, cov_stop[-N] + 1) 
	cov_size = cov_stop[N]
	visits_start_each = c(1, visits_stop_each[-length(visits_stop_each)] + 1)

	M = sum(V)

	# check knots placement: knots
	return(list('num_subjects'=N, 'num_blocks'=P, 'num_PCs'=K, 'num_basis'=Q, 'response_vector'=unlist(Y), 'num_visits_all'=M,
		        'subject_starts' = visits_start, 'subject_stops' = visits_stop, 
		        'cov_starts' = cov_start, 'cov_stops' = cov_stop, 'cov_size' = cov_size, 
		        'subject_starts_each' = visits_start_each, 'subject_stops_each' = visits_stop_each,
		        'visits_matrix' = V, 'spline_basis' = phi_t_stacked, 'len_Theta' = len_Theta, 
		        'Theta_idx' = Theta_idx, 'Q_idx_Theta' = Q_idx_Theta, 'K_idx_Theta' = K_idx_Theta,
		        'nknots' = nknots, 'knots' = knots, 'phi_t_cont' = phi_t_cont, 'phi_t' = phi_t, 'time_cont' = time_cont))
}


#--------------------------------------------
# set up basis
#--------------------------------------------
basis_setup_sparse=function(time_sparse,nknots=rep(1,length(time_sparse[[1]])), orth=TRUE, delta=1/1000){

	library(splines)
	library(pracma)

	N=length(time_sparse)
	P=length(time_sparse[[1]])
	time_unique=list()
	for(p in 1:P){
		for(i in 1:N){
			time_sparse[[i]][[p]]=round(time_sparse[[i]][[p]]/delta)*delta
		}	
		time_unique[[p]]=time_sparse[[1]][[p]]
		for(i in 2:N){
			time_unique[[p]]=c(time_unique[[p]],time_sparse[[i]][[p]])
		}
		time_unique[[p]]=round(sort(unique(time_unique[[p]]))/delta)*delta
	}	
	time_unique_combined=time_unique[[1]]
	for(p in 2:P){
		time_unique_combined=c(time_unique_combined,time_unique[[p]])
	}
	time_unique_combined=sort(unique(time_unique_combined))
	
	T_min=min(time_unique_combined)
	T_max=max(time_unique_combined)
	time_cont=seq(T_min,T_max/delta)*delta
	time_cont=round(time_cont/delta)*delta
	
	knots=list()
	for(p in 1:P){
		K=nknots[p]+4
		qs=1/(nknots[p]+1)
		knots[[p]]=quantile(time_unique[[p]],qs)
		if(nknots[p]>1){
			for(q in 2:nknots[p]){
				knots[[p]]=c(knots[[p]],q*quantile(time_unique[[p]],qs))
			}
		}	
		knots[[p]]=as.vector(knots[[p]])
	}

	phi_t_cont=list()
	for(p in 1:P){
		phi_t_cont[[p]]=bs(time_cont,knots=knots[[p]],degree=3,intercept=TRUE)
		temp=phi_t_cont[[p]]
		for(k in 1:(nknots[p]+4)){
			if(orth==TRUE){
				if(k>1){
					for(q in 1:(k-1)){
						temp[,k]=temp[,k]-(sum(temp[,k]*temp[,k-q])/
							sum(temp[,k-q]^2))*temp[,k-q];
					}
				}
			}		
		    temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
		}
		phi_t_cont[[p]]=t(sqrt(1/delta)*temp)
	}	

	phi_t=list()
	for(p in 1:P){
		phi_t[[p]]=list()
		for(i in 1:N){
			phi_t[[p]][[i]]=array(0,dim=c((nknots[p]+4),length(time_sparse[[i]][[p]])))
			for(k in 1:(nknots[p]+4)){
				for(t in 1:length(time_sparse[[i]][[p]])){
					phi_t[[p]][[i]][k,t]=phi_t_cont[[p]][k,abs(time_cont-time_sparse[[i]][[p]][t])==
						min(abs(time_cont-time_sparse[[i]][[p]][t]))][1]
				}
			}
		}		
	}
			

	results=list()
	results[[1]]=knots
	results[[2]]=time_sparse
	results[[3]]=phi_t
	results[[4]]=time_cont
	results[[5]]=phi_t_cont

	return(results)	
}	


#-------------------------------------------------
# post hoc rotation on MFPCA results
# fit: model returns from mfpca
# Nchains: #chains used in stan
# Nsamples: #samples used in stan
#-------------------------------------------------

library(Matrix)
library(rstan)
library(corpcor) # for covariance decomposition

post_hoc_rotation = function(fit, Nchains, Nsamples, N, P, K, Q){
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
	prop_var = list()	
	for (i in 1:dim(Theta)[1]){
		for (j in 1:dim(Theta)[2]){
			ind = ind + 1 # loop through all posterior draws
			Theta_old[,,ind] = array(Theta[i,j,],dim=c(sum(Q), sum(K)))
			cov_old[,,ind] = array(cov[i,j,], dim=c(sum(K), sum(K)))
			theta_mu_new[,ind] = array(theta_mu[i,j,])
			alpha_old[,,ind] = array(alpha[i,j,],dim=c(sum(K), N))

			indq = 1;
	    	indk = 1;
	    	prop_var[[ind]] = list()
	    	for (p in 1:P){
		    	poolvar = as.matrix(cov_old[,,ind][indk:sum(K[1:p]),indk:sum(K[1:p])]) 
		    	temp_sigma = Theta_old[,,ind][indq:sum(Q[1:p]),indk:sum(K[1:p])] %*% poolvar %*% t(Theta_old[,,ind][indq:sum(Q[1:p]),indk:sum(K[1:p])])
		    	eigen_temp_sigma=eigen(temp_sigma)
				v_temp=Re(eigen_temp_sigma$vectors)
				d_temp=Re(eigen_temp_sigma$values)	
				prop_var[[ind]][[p]] = d_temp / sum(d_temp)		

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

	# compute averaged explained variance (average across draws)
	prop_var_avg = list()
	for (p in 1:P){
		tmp = sapply(prop_var, '[[', p)
		prop_var_avg[[p]] = rowMeans(tmp)[1:K[p]]
	}

	return(list('ALPHA_array' = ALPHA_array, 'MU_array'= MU_array, 'THETA_array'= THETA_array, 
                'COV_array' = COV_array, 'COR_array'= COR_array, 
                'prop_var' = prop_var, 'prop_var_avg'= prop_var_avg))
}

#-------------------------------------------------------------
# alternative plotting: save results + plot each individually
#-------------------------------------------------------------
output_mfpca = function(dat_mfpca, param_stan, post_rotation_results, title){
	K = param_stan$num_PCs
	Q = param_stan$num_basis
	N = param_stan$num_subjects
	P = param_stan$num_blocks
	nknots = param_stan$nknots
	time_origin = dat_mfpca$time_origin
	mu_y = dat_mfpca$mu_y_origin
	sigma_y = dat_mfpca$sigma_y_origin
	ids_origin = dat_mfpca$shared_id_unique
	covariates = dat_mfpca$covariates

	Y_sparse = dat_mfpca$Y_sparse
	time_sparse = dat_mfpca$time_sparse
	phi_t_cont = param_stan$phi_t_cont
	phi_t = param_stan$phi_t
	time_cont = param_stan$time_cont
	ALPHA_array = post_rotation_results$ALPHA_array
	MU_array = post_rotation_results$MU_array
	THETA_array = post_rotation_results$THETA_array
	COV_array = post_rotation_results$COV_array
	COR_array = post_rotation_results$COR_array	 
	prop_var_avg = post_rotation_results$prop_var_avg 


	nloop=dim(ALPHA_array)[3]
	first=1
	last=nloop
	MU_mean = MU_array[, first] # mean function across sampling sessions
	ALPHA_mean = ALPHA_array[,,first] # mean factor scores
	THETA_mean = THETA_array[,,first] # mean factor loading
	COV_mean = COV_array[,,first]
	COR_mean = COR_array[,,first]    

	for(iter in 2:nloop){
		MU_mean = MU_mean + MU_array[, iter]
		ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
		THETA_mean = THETA_mean + THETA_array[,,iter]
		COV_mean = COV_mean + COV_array[,,iter]
		COR_mean = COR_mean + COR_array[,,iter]        
	}

	MU_mean = cbind(MU_mean/(last-first+1)) 
	ALPHA_mean = cbind(ALPHA_mean/(last-first+1))
	THETA_mean = cbind(THETA_mean/(last-first+1)) 
	COV_mean = cbind(COV_mean/(last-first+1)) 
	COR_mean = cbind(COR_mean/(last-first+1))     

	tmp = bdiag(cbind(phi_t_cont[[1]]), cbind(phi_t_cont[[2]]))
	if(P > 2){
		for(p in 3:P){
			tmp = bdiag(tmp,cbind(phi_t_cont[[p]]))
		}
	}
	Mu_functions=t(tmp)%*%MU_mean

	Mu = list()
	for(p in 1:P){
		Mu[[p]] = Mu_functions[((p-1)*length(time_cont)+1):(p*length(time_cont))]
	}

	FPC_mean = list()
	ind1 = 1
	ind2 = 1
	for(p in 1:P){
		FPC_mean[[p]] = t(phi_t_cont[[p]])%*%THETA_mean[ind1:(ind1+nknots[p]+3),ind2:(ind2+K[p]-1)]
		ind1 = ind1 + nknots[p]+4
		ind2 = ind2 + K[p]
	}

	#setting the x- and y-axis limit for plotting
	Y_min = list()
	Y_max = list()
	for(p in 1:P){
		Y_min[[p]]=min(Y_sparse[[1]][[p]])	
		Y_max[[p]]=max(Y_sparse[[1]][[p]])	
		for(i in 2:N){
			Y_min[[p]]=min(Y_min[[p]],min(Y_sparse[[i]][[p]]))
			Y_max[[p]]=max(Y_max[[p]],max(Y_sparse[[i]][[p]]))
		}
	}

	#---------------------------------------
	# output FPC scores with original ID
	#---------------------------------------
	table_scores = as.data.frame(t(ALPHA_mean))
	names_scores = NULL
	names_covs = NULL
	for (p in 1:P){
		for (k in 1:K[p]){
			names_scores = c(names_scores, paste0(title[p], '_FPC', k))
		}
	}
	colnames(table_scores) = names_scores

	df = dat_mfpca$data[, c('ID', ids_origin, covariates)]
	df = df %>% distinct(!!as.symbol(ids_origin), .keep_all=TRUE)
	table_scores = cbind(df, table_scores)

	return(list('P' = P, 'N' = N, 'K' = K, 
		        'time_cont' = time_cont, 'time_origin' = time_origin, 'time_sparse' = time_sparse,
		        'Mu' = Mu, 'ALPHA_mean' = ALPHA_mean, 'FPC_mean' = FPC_mean, 'prop_var_avg' = prop_var_avg,
		        'Y_min' = Y_min, 'Y_max' = Y_max, 'Y_sparse' = Y_sparse, 
		        'mu_y' = mu_y, 'sigma_y' = sigma_y, 'FPC_scores' = table_scores, 'title' = title))
}


plots_mfpca = function(output_mfpca, plot_mean, mean_x_label, mean_y_label, 
					   mean_x_range, mean_y_ticks, mean_x_ticks, mean_title,
					   plot_fpc, fpc_y_lim, fpc_x_label, fpc_y_label, fpc_title, fpc_y_ticks, fpc_x_ticks, 
					   plot_fpc_scores, scores_selected=NULL, scores_group_compare=NULL, scores_y_label=NULL){

	P = output_mfpca$P
	N = output_mfpca$N
	K = output_mfpca$K
	time_cont = output_mfpca$time_cont
	time_origin = output_mfpca$time_origin
	time_sparse = output_mfpca$time_sparse
	Mu = output_mfpca$Mu
	sigma_y = output_mfpca$sigma_y
	mu_y = output_mfpca$mu_y
	Y_min = output_mfpca$Y_min
	Y_max = output_mfpca$Y_max
	Y_sparse = output_mfpca$Y_sparse
	FPC_mean = output_mfpca$FPC_mean
	table_scores = output_mfpca$FPC_scores

	if (plot_mean){
		figs_mean_list = list()
		for(p in 1:P){
			plot(time_cont*max(time_origin),Mu[[p]]*sigma_y[p] + mu_y[p],type="l", yaxt="n", xaxt="n",
				 ylim=c(Y_min[[p]]*sigma_y[p] + mu_y[p], Y_max[[p]]*sigma_y[p] + mu_y[p]),
				 xlim=mean_x_range[[p]],lwd=5,col=4, ylab=mean_y_label[p], xlab=mean_x_label, font.lab=2, cex.lab=1.2)
			for(i in 1:N){
				lines(time_sparse[[i]][[p]]*max(time_origin),Y_sparse[[i]][[p]]*sigma_y[p] + mu_y[p],type="l",lwd=.25)
			}
			title(main=mean_title[p])
			axis(2, at = mean_y_ticks[[p]], las = 2)
			axis(1, at = mean_x_ticks[[p]], las = 1)
        	assign(paste0('figure_mean', p), recordPlot())
        	figs_mean_list[[p]] = eval(as.name(paste0('figure_mean', p)))
		}
	}else{figs_mean_list = NULL}

	if (plot_fpc){
		figs_fpc_list = list()
		for(p in 1:P){
		    tmp = NULL
		    tmp = cbind(tmp, time_cont*max(time_origin))
		    for (k in 1:K[[p]]){
		        tmp = cbind(tmp, FPC_mean[[p]][,k]*sigma_y[p] + mu_y[p])
		    }
		    tmp = as.data.frame(tmp)
		    colnames(tmp) = c('time', paste0('fpc', seq(1, K[[p]], 1)))
		    
		    df = reshape::melt(tmp, 'time', paste0('fpc', seq(1, K[[p]], 1)))
		    
		    figs_fpc_list[[p]] <- ggplot() + geom_line(data = df, aes(x = time, y = value, 
		                                     colour = variable, linetype = variable), lwd = 1) +
		                          guides(linetype = F) + labs(colour = 'curve') + 
		                          labs(title = fpc_title[p], x = fpc_x_label, y = fpc_y_label[p]) +
		                          scale_x_continuous(breaks = fpc_x_ticks[[p]]) +
		                          scale_y_continuous(limits = fpc_y_lim[[p]], breaks = fpc_y_ticks[[p]]) +
		                          theme_classic() +
		                          theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
		                                axis.text.x = element_text(size = 10, face = "bold"),
		                                axis.text.y = element_text(size = 10, face = "bold"),
		                                axis.title.x = element_text(size = 12, face = "bold"),
		                                axis.title.y = element_text(size = 12, face = "bold"),
		                                legend.title = element_blank(), legend.position = 'top')   
    
		}  
	}else{figs_fpc_list = NULL}  		

if (plot_fpc_scores){
		figs_scores_list = list()
		for (p in 1:P){
			if (length(scores_group_compare) == 2){
				tmp <- ggplot(aes(y = eval(as.name(scores_selected[[p]])), 
		                                    x = eval(as.name(scores_group_compare[[1]])), 
		                                    fill = eval(as.name(scores_group_compare[[2]]))), 
		                                    data = table_scores) + 
		    						 geom_boxplot() + 
		                             ylab(scores_y_label[p]) + xlab(scores_group_compare[[1]]) +
		                             theme_classic() + theme(axis.text = element_text(size=10, face="bold", color='black'),
		                                  axis.title.x = element_text(size=12, face="bold"),
		                                  axis.title.y = element_text(size=12, face="bold"),
		                                  legend.position = "top") + scale_fill_discrete(name = scores_group_compare[[2]])
		    }else{
				tmp <- ggplot(aes(y = eval(as.name(scores_selected[[p]])), 
		                                    x = eval(as.name(scores_group_compare)), 
		                                    fill = eval(as.name(scores_group_compare))),
		                                    data = table_scores) + 
		    						 geom_boxplot() + 
		                             ylab(scores_y_label[p]) + xlab(scores_group_compare) +
		                             theme_classic() + theme(axis.text = element_text(size=10, face="bold", color='black'),
		                                  axis.title.x = element_text(size=12, face="bold"),
		                                  axis.title.y = element_text(size=12, face="bold"),
		                                  legend.position = "none") 
		                             # + geom_jitter(aes(colour=eval(as.name(scores_group_compare))), shape=16, 
		                             # 	         position=position_jitter(0.2)) # color weird
		    }
		     
		    figs_scores_list[[p]] <- ggplotGrob(tmp)
	    }  
	}else{figs_scores_list = NULL}                           	
	return(list('figs_mean_list' = figs_mean_list, 'figs_fpc_list' = figs_fpc_list, 'figs_scores_list' = figs_scores_list))
}

                 
#------------------------------------------------------------
# plot mfpca results including
# 1. standaridzed overal mean curve
# 2. original scaled mean curve
# 3. original scaled FPC curve (w/wo 1SD FPC scores)
#------------------------------------------------------------
library(dplyr) # for remove duplicated rows in data frame with distinct()
visu_mfpca = function(dat_mfpca, param_stan, post_rotation_results, title,
					  x_label, y_label, x_range, FPC_sd, fig_dir){

	K = param_stan$num_PCs
	Q = param_stan$num_basis
	N = param_stan$num_subjects
	P = param_stan$num_blocks
	nknots = param_stan$nknots
	time_origin = dat_mfpca$time_origin
	mu_y = dat_mfpca$mu_y_origin
	sigma_y = dat_mfpca$sigma_y_origin
	ids_origin = dat_mfpca$shared_id_unique
	covariates = dat_mfpca$covariates

	Y_sparse = dat_mfpca$Y_sparse
	time_sparse = dat_mfpca$time_sparse
	phi_t_cont = param_stan$phi_t_cont
	phi_t = param_stan$phi_t
	time_cont = param_stan$time_cont
	ALPHA_array = post_rotation_results$ALPHA_array
	MU_array = post_rotation_results$MU_array
	THETA_array = post_rotation_results$THETA_array
	COV_array = post_rotation_results$COV_array
	COR_array = post_rotation_results$COR_array	 
	prop_var_avg = post_rotation_results$prop_var_avg 


	nloop=dim(ALPHA_array)[3]
	first=1
	last=nloop
	MU_mean = MU_array[, first] # mean function across sampling sessions
	ALPHA_mean = ALPHA_array[,,first] # mean factor scores
	THETA_mean = THETA_array[,,first] # mean factor loading
	COV_mean = COV_array[,,first]
	COR_mean = COR_array[,,first]    

	for(iter in 2:nloop){
		MU_mean = MU_mean + MU_array[, iter]
		ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
		THETA_mean = THETA_mean + THETA_array[,,iter]
		COV_mean = COV_mean + COV_array[,,iter]
		COR_mean = COR_mean + COR_array[,,iter]        
	}

	MU_mean = cbind(MU_mean/(last-first+1)) 
	ALPHA_mean = cbind(ALPHA_mean/(last-first+1))
	THETA_mean = cbind(THETA_mean/(last-first+1)) 
	COV_mean = cbind(COV_mean/(last-first+1)) 
	COR_mean = cbind(COR_mean/(last-first+1))     

	tmp = bdiag(cbind(phi_t_cont[[1]]), cbind(phi_t_cont[[2]]))
	if(P > 2){
		for(p in 3:P){
			tmp = bdiag(tmp,cbind(phi_t_cont[[p]]))
		}
	}
	Mu_functions=t(tmp)%*%MU_mean

	Mu = list()
	for(p in 1:P){
		Mu[[p]] = Mu_functions[((p-1)*length(time_cont)+1):(p*length(time_cont))]
	}

	FPC_mean = list()
	ind1 = 1
	ind2 = 1
	for(p in 1:P){
		FPC_mean[[p]] = t(phi_t_cont[[p]])%*%THETA_mean[ind1:(ind1+nknots[p]+3),ind2:(ind2+K[p]-1)]
		ind1 = ind1 + nknots[p]+4
		ind2 = ind2 + K[p]
	}

	#setting the x- and y-axis limit for plotting
	Y_min = list()
	Y_max = list()
	for(p in 1:P){
		Y_min[[p]]=min(Y_sparse[[1]][[p]])	
		Y_max[[p]]=max(Y_sparse[[1]][[p]])	
		for(i in 2:N){
			Y_min[[p]]=min(Y_min[[p]],min(Y_sparse[[i]][[p]]))
			Y_max[[p]]=max(Y_max[[p]],max(Y_sparse[[i]][[p]]))
		}
	}

	#----------------------------------------------------
	# fig: spaghetti + mean curve (standardized scale)
	#----------------------------------------------------
    #capabilities()
    print('figure mean at transformed scale')
	png(paste0(fig_dir, 'mean_spaghetti_scaled.png'), units="in", width=10, height=7, res=300)
    par(mfrow=c(ceiling(P/2),2))
	for(p in 1:P){
		plot(time_cont*max(time_origin),Mu[[p]], type="l", ylim=c(Y_min[[p]],Y_max[[p]]),
			 xlim=x_range, lwd=5, col=4, ylab=y_label, xlab=x_label)
		for(i in 1:N){
			lines(time_sparse[[i]][[p]]*max(time_origin), Y_sparse[[i]][[p]], type="l", lwd=.25)
		}
		#title(main=names(dat)[ind_y[p]])
		title(main=title[p])
	}
	dev.off()

	#----------------------------------------------------
	# fig: spaghetti + mean curve (transformed scale)
	#----------------------------------------------------
    print('figure mean at original scale')
	png(paste0(fig_dir,'mean_spaghetti_original.png'), units="in", width=10, height=7, res=300)
    par(mfrow=c(ceiling(P/2),2))
	for(p in 1:P){
		plot(time_cont*max(time_origin),Mu[[p]]*sigma_y[p] + mu_y[p],type="l",
			 ylim=c(Y_min[[p]]*sigma_y[p] + mu_y[p], Y_max[[p]]*sigma_y[p] + mu_y[p]),
			 xlim=x_range,lwd=5,col=4, ylab=y_label, xlab=x_label, font.lab=2, cex.lab=1.2)
		for(i in 1:N){
			lines(time_sparse[[i]][[p]]*max(time_origin),Y_sparse[[i]][[p]]*sigma_y[p] + mu_y[p],type="l",lwd=.25)
		}
		title(main=title[p])
	}
	dev.off()

	#-------------------------------------------------------------------
	# fig: FPC curves +- 1SD (FPC scores) (transformed scale)
	#-------------------------------------------------------------------
    print('figure fpc curves')
	idx_alpha = 0
	for(p in 1:P){
		png(paste0(fig_dir, 'PCs_original_', title[p], '_1sd_', FPC_sd, '.png'), units="in", width=7, height=5, res=300)
        
        if (p <= 2){
            par(mfrow=c(2,2))
        }else{
            par(mfrow=c(ceiling(p/2),2))
        }    
		
		for (k in 1:K[[p]]){
			idx_alpha = idx_alpha + 1
			prop_var_approx = paste0(round(prop_var_avg[[p]][[k]]*100, 2), '%')
			alpha_1sd = sd(ALPHA_mean[idx_alpha, ])

			plot(time_cont*max(time_origin), Mu[[p]]*sigma_y[p] + mu_y[p], type="l", 
				ylim=c(Y_min[[p]]*sigma_y[p] + mu_y[p], Y_max[[p]]*sigma_y[p] + mu_y[p]), 
		 		lwd=2,col=1, ylab=y_label, xlab=x_label, font.lab=2, cex.lab=1.2)

			if (FPC_sd){ # +- SD of FPC scores
				lines(time_cont*max(time_origin), 
			      (Mu[[p]] + FPC_mean[[p]][,k])*alpha_1sd*sigma_y[p] + mu_y[p],type="l",lwd=3,lty=2,col=2) # red
				lines(time_cont*max(time_origin), 
			      (Mu[[p]] - FPC_mean[[p]][,k])*alpha_1sd*sigma_y[p] + mu_y[p],type="l",lwd=3,lty=2,col=3) # green
			}else{
				lines(time_cont*max(time_origin), 
			      (Mu[[p]] + FPC_mean[[p]][,k])*sigma_y[p] + mu_y[p],type="l",lwd=3,lty=2,col=2) # red
				lines(time_cont*max(time_origin), 
			      (Mu[[p]] - FPC_mean[[p]][,k])*sigma_y[p] + mu_y[p],type="l",lwd=3,lty=2,col=3) # green
			}

			title(main=paste(paste('PC', k, sep=' '), ' (', prop_var_approx, ' )', sep=''))
			#axis(1, font=2) # make x-axis ticks label bold
			legend('topright', c('+ pc', '- pc'), lty=c(2,2), lwd=c(3,3), col=c(2, 3), bty='n', cex=0.5)
		}	
		dev.off()
	}

	#---------------------------------------
	# output FPC scores with original ID
	#---------------------------------------
	table_scores = as.data.frame(t(ALPHA_mean))
	names_scores = NULL
	names_covs = NULL
	for (p in 1:P){
		for (k in 1:K[p]){
			names_scores = c(names_scores, paste0(title[p], '_FPC', k))
		}
	}
	colnames(table_scores) = names_scores

# 	df = dat_mfpca$data[, c('ID', 'ID_unique', covariates)]
# 	df = df %>% distinct(ID_unique, .keep_all=TRUE)
    df = dat_mfpca$data[, c('ID', ids_origin, covariates)]
    df = df %>% distinct(!!as.symbol(ids_origin), .keep_all=TRUE)
	table_scores = cbind(df, table_scores)

	return(list('FPC_scores' = table_scores))
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

# updated MI calculation 
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

    pairs = combn(seq(1,P, 1), 2)
    pair_names_MI = pair_names_CMI = pair_names_MI_norm = pair_names_CMI_norm = NULL
    for (p in 1:dim(pairs)[2]){
        pair_names_MI = c(pair_names_MI, paste0('MI_', pairs[1, p], pairs[2, p]))
        pair_names_CMI = c(pair_names_CMI, paste0('CMI_', pairs[1, p], pairs[2, p]))
        pair_names_MI_norm = c(pair_names_MI_norm, paste0('norm_MI_', pairs[1, p], pairs[2, p]))
        pair_names_CMI_norm = c(pair_names_CMI_norm, paste0('norm_CMI_', pairs[1, p], pairs[2, p]))

        # compute mutual information
        assign(paste0('MI_', pairs[1, p], pairs[2, p]), - H_est(R[unlist(idx_each[c(pairs[1, p], pairs[2, p])]), 
                   unlist(idx_each[c(pairs[1, p], pairs[2, p])])]))

        # compute conditional mutual information
        assign(paste0('CMI_', pairs[1, p], pairs[2, p]), 
               H_est(R_out_each[[pairs[1, p]]]) + H_est(R_out_each[[pairs[2, p]]]) - 
               H_est(R[-unlist(idx_each[c(pairs[1, p], pairs[2, p])]), 
                       -unlist(idx_each[c(pairs[1, p], pairs[2, p])])]) - H_est(R))
        
        # compute normalized MI
        assign(paste0('norm_MI_', pairs[1, p], pairs[2, p]), 
               MI_norm(eval(as.name(paste0('MI_', pairs[1, p], pairs[2, p])))))
        assign(paste0('norm_CMI_', pairs[1, p], pairs[2, p]), 
               MI_norm(eval(as.name(paste0('CMI_', pairs[1, p], pairs[2, p])))))
    }
    
    results_list = list()
    for (names in c(pair_names_MI, pair_names_CMI, pair_names_MI_norm, pair_names_CMI_norm)){
        results_list[[names]] = eval(as.name(names))
    }
    return(results_list)
}

MI_est_data = function(post_rotation_results, K){
	# R: true correlation matrix from simulated data
	
    R_array = post_rotation_results$COR_array
	nloop = dim(R_array)[3]
	P = length(K)

    MI_results = list()
	for (i in 1:nloop){
        MI_results[[i]] = MI_est(corr_matrix=R_array[,,i], K=K)
    }       
    
    list_names = MI_names = names(MI_results[[1]])
    summary_names = paste0(MI_names, '_list')
    results_list = list()
    for (m in 1:length(MI_names)){

        tmp = list()
        for (i in 1:nloop){
            tmp[i] = MI_results[[i]][MI_names[m]]
        }

        assign(list_names[m], tmp)

        results_tmp = as.numeric(eval(as.name(list_names[m])))
        assign(summary_names[m], c(mean(results_tmp), quantile(results_tmp, 0.5), 
                                   quantile(results_tmp, 0.025), quantile(results_tmp, 0.975)))
        results_list[[summary_names[m]]] = eval(as.name(summary_names[m]))
    }
    
    return(results_list)
}


#-----------------------------------
# diagnostic plot - LOOIC
#-----------------------------------
diagnostic_looic = function(loo_list, num_subjects, fig_dir=NULL){
	N = num_subjects
	pkdf <- data.frame(pk=loo_list$diagnostics$pareto_k, id=1:N)
	fig_looic <- ggplot(pkdf, aes(x=id, y=pk)) + geom_point(shape=3, color="blue") +
	    labs(x="Observation left out", y="Pareto shape k") +
	    geom_hline(yintercept = 0.7, linetype=2, color="red", size=0.2) +
	    ggtitle("PSIS-LOO diagnostics") + theme_classic() + 
	    theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"),
	    	  axis.text.x= element_text(size=10, face="bold"),
	    	  axis.text.y= element_text(size=10, face="bold"),
		 	  axis.title.x= element_text(size=12, face="bold"),
		      axis.title.y= element_text(size=12, face="bold")) 

	if (!is.null(fig_dir)){    
		#ggsave(paste0(fig_dir, 'diagnostic_looic_fig.pdf'), width=4, height=4, dpi=300)
		ggsave(paste0(fig_dir, 'diagnostic_looic_fig.png'), width=4, height=4, dpi=300)
	}	

	return(fig_looic)
}


#---------------------------------------------------------
# diagnostic plot - posterior predictive checking
#---------------------------------------------------------
library("bayesplot")
library("ggplot2")
diagnostic_posterior = function(mfpca_fit, Nsamples, Nchains, visits_matrix, response_observed, fig_title, fig_dir=NULL){
	Ynew = extract(mfpca_fit,"Ynew",permuted=FALSE)
	V = visits_matrix
	Ynew_transform = matrix(rep(0, Nsamples/2 * Nchains * sum(V)), ncol=sum(V))
	ind = 0
	for (i in 1:(Nsamples/2)){
		for (j in 1:Nchains){
			ind = ind + 1
			Ynew_transform[ind, ] = Ynew[i,j,]
		}
	}
	Ynew_mean = colMeans(Ynew_transform)

	color_scheme_set("brightblue")
	fig_posterior <- ppc_dens_overlay(response_observed, Ynew_transform) + 
					 ggtitle(fig_title) + labs(x="Standardized response", y="Kernel density") +
					 theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"),
					 	   axis.text.x= element_text(size=10, face="bold"),
				    	   axis.text.y= element_text(size=10, face="bold"),
					 	   axis.title.x= element_text(size=12, face="bold"),
					       axis.title.y= element_text(size=12, face="bold"),
					 	   legend.position = 'top')

	if (!is.null(fig_dir)){
		#ggsave(paste0(fig_dir, 'diagnostic_posterior_fig.pdf'), width=4, height=4, dpi=300)
		ggsave(paste0(fig_dir, 'diagnostic_posterior_fig.png'), width=4, height=4, dpi=300)		
	}

	return(fig_posterior)
}

