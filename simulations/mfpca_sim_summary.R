#--------------------------------------------------------
# calculate coverage probabilities for covariance matrix
#--------------------------------------------------------
coverage_cov_95 = function(dir_home, file_data, file_results, nsims){
    dir_data = paste0(dir_home, file_data)
    
    load(dir_data, dat <- new.env())
    cov_true = dat$sim_param$cov_alpha_true
    cov_params = rep(0, dim(cov_true)[1] * dim(cov_true)[2])
    idx = 0
    for (nrow in 1:dim(cov_true)[1]){
        for (ncol in 1:dim(cov_true)[2]){
            idx = idx + 1
            cov_params[idx] = paste('cov_alpha[', nrow, ',', ncol, ']', sep='')
        }
    }

    cov_median_list_weak = cov_q025_list_weak = cov_q975_list_weak = cov_cov_95_weak = list()
    cov_coverage_95_weak = matrix(rep(0, nsims*length(cov_params)), nrow=nsims)
    colnames(cov_coverage_95_weak) = cov_params
    error_nsim_weak = NULL
    for (i in 1:nsims){
        dir_results = paste0(dir_home, file_results, i, '.Rdata')
        tmp = try(load(dir_results), TRUE)
        if(isTRUE(class(tmp)=="try-error")) { 
            error_nsim_weak = c(error_nsim_weak, i);
            cov_coverage_95_weak[i,] = NA
            next } 

        cov_median_list_weak[[i]] = results_nsim$est_results$COV_median
        cov_q025_list_weak[[i]] = results_nsim$est_results$COV_q025
        cov_q975_list_weak[[i]] = results_nsim$est_results$COV_q975

        idx_cov = 0
        for (nrow in 1:dim(cov_true)[1]){
            for (ncol in 1:dim(cov_true)[2]){
                idx_cov = idx_cov + 1
                # rounding error with zero within-block covariance
                if (cov_true[nrow, ncol] <= round(cov_q975_list_weak[[i]][nrow, ncol], 5) & 
                    cov_true[nrow, ncol] >= round(cov_q025_list_weak[[i]][nrow, ncol], 5)){
                cov_coverage_95_weak[i, idx_cov] = 1
                } 
            }	
        }

    }
    
    cov_coverage_95_prob_weak = colMeans(cov_coverage_95_weak, na.rm=T)

    return(list('cov_true'=cov_true, 'cov_coverage_95_weak'=cov_coverage_95_weak, 
                'cov_coverage_95_prob_weak'=cov_coverage_95_prob_weak, 
                'cov_median_list_weak'=cov_median_list_weak,
                'cov_q025_list_weak'=cov_q025_list_weak, 
                'cov_q975_list_weak'=cov_q975_list_weak))
}


#--------------------------------------------------------
# calculate coverage probabilities for correlation matrix
#--------------------------------------------------------
coverage_cor_95 = function(dir_home, file_data, file_results, nsims){
    dir_data = paste0(dir_home, file_data)
    
    load(dir_data, dat <- new.env())
    cor_true = dat$sim_param$params[[3]]
    cor_params = rep(0, dim(cor_true)[1] * dim(cor_true)[2])
    idx = 0
    for (nrow in 1:dim(cor_true)[1]){
        for (ncol in 1:dim(cor_true)[2]){
            idx = idx + 1
            cor_params[idx] = paste('cor_alpha[', nrow, ',', ncol, ']', sep='')
        }
    }

    cor_median_list_weak = cor_q025_list_weak = cor_q975_list_weak = cor_cor_95_weak = list()
    cor_coverage_95_weak = matrix(rep(0, nsims*length(cor_params)), nrow=nsims)
    colnames(cor_coverage_95_weak) = cor_params
    error_nsim_weak = NULL
    for (i in 1:nsims){
        dir_results = paste0(dir_home, file_results, i, '.Rdata')
        tmp = try(load(dir_results), TRUE)
        if(isTRUE(class(tmp)=="try-error")) { 
            error_nsim_weak = c(error_nsim_weak, i);
            cor_coverage_95_weak[i,] = NA
            next } 

        cor_median_list_weak[[i]] = results_nsim$est_results$COR_median
        cor_q025_list_weak[[i]] = results_nsim$est_results$COR_q025
        cor_q975_list_weak[[i]] = results_nsim$est_results$COR_q975

        idx_cov = 0
        for (nrow in 1:dim(cor_true)[1]){
            for (ncol in 1:dim(cor_true)[2]){
                idx_cov = idx_cov + 1
                # rounding error with zero within-block covariance
                if (cor_true[nrow, ncol] <= round(cor_q975_list_weak[[i]][nrow, ncol], 5) & 
                    cor_true[nrow, ncol] >= round(cor_q025_list_weak[[i]][nrow, ncol], 5)){
                cor_coverage_95_weak[i, idx_cov] = 1
                } 
            }	
        }

    }
    
    cor_coverage_95_prob_weak = colMeans(cor_coverage_95_weak, na.rm=T)

    return(list('cor_true'=cor_true, 'cor_coverage_95_weak'=cor_coverage_95_weak, 
                'cor_coverage_95_prob_weak'=cor_coverage_95_prob_weak, 
                'cor_median_list_weak'=cor_median_list_weak,
                'cor_q025_list_weak'=cor_q025_list_weak, 
                'cor_q975_list_weak'=cor_q975_list_weak))
}

#------------------------------------------------------------------------
# calculate coverage probabilities for mutual information (updated
#------------------------------------------------------------------------
coverage_MI_95 = function(dir_home, file_results, nsims){
    MI_params = c('MI_12', 'MI_13', 'MI_23', 'CMI_12', 'CMI_13', 'CMI_23',
                  'norm_MI_12', 'norm_MI_13', 'norm_MI_23',
                  'norm_CMI_12', 'norm_CMI_13', 'norm_CMI_23')

    MI_median_list_weak = MI_q025_list_weak = MI_q975_list_weak = MI_MI_95_weak = list()
    MI_coverage_95_weak = matrix(rep(0, nsims*length(MI_params)), nrow=nsims)
    colnames(MI_coverage_95_weak) = MI_params
    error_nsim_weak = NULL

    for (i in 1:nsims){
        dir_results = paste0(dir_home, file_results, i, '.Rdata')
        tmp = try(load(dir_results), TRUE)
        if(isTRUE(class(tmp)=="try-error")) { 
            error_nsim_weak = c(error_nsim_weak, i);
            MI_coverage_95_weak[i,] = NA
            next } 

        MI_median_list_weak[[i]] = c(results_nsim$MI_results$MI_12_list[['50%']],
                                     results_nsim$MI_results$MI_13_list[['50%']],
                                     results_nsim$MI_results$MI_23_list[['50%']],
                                     results_nsim$MI_results$CMI_12_list[['50%']],
                                     results_nsim$MI_results$CMI_13_list[['50%']],
                                     results_nsim$MI_results$CMI_23_list[['50%']],
                                     results_nsim$MI_results$norm_MI_12_list[['50%']],
                                     results_nsim$MI_results$norm_MI_13_list[['50%']],
                                     results_nsim$MI_results$norm_MI_23_list[['50%']],
                                     results_nsim$MI_results$norm_CMI_12_list[['50%']],
                                     results_nsim$MI_results$norm_CMI_13_list[['50%']],
                                     results_nsim$MI_results$norm_CMI_23_list[['50%']])
        MI_q025_list_weak[[i]] = c(results_nsim$MI_results$MI_12_list[['2.5%']],
                                     results_nsim$MI_results$MI_13_list[['2.5%']],
                                     results_nsim$MI_results$MI_23_list[['2.5%']],
                                     results_nsim$MI_results$CMI_12_list[['2.5%']],
                                     results_nsim$MI_results$CMI_13_list[['2.5%']],
                                     results_nsim$MI_results$CMI_23_list[['2.5%']],
                                     results_nsim$MI_results$norm_MI_12_list[['2.5%']],
                                     results_nsim$MI_results$norm_MI_13_list[['2.5%']],
                                     results_nsim$MI_results$norm_MI_23_list[['2.5%']],
                                     results_nsim$MI_results$norm_CMI_12_list[['2.5%']],
                                     results_nsim$MI_results$norm_CMI_13_list[['2.5%']],
                                     results_nsim$MI_results$norm_CMI_23_list[['2.5%']])
        MI_q975_list_weak[[i]] = c(results_nsim$MI_results$MI_12_list[['97.5%']],
                                     results_nsim$MI_results$MI_13_list[['97.5%']],
                                     results_nsim$MI_results$MI_23_list[['97.5%']],
                                     results_nsim$MI_results$CMI_12_list[['97.5%']],
                                     results_nsim$MI_results$CMI_13_list[['97.5%']],
                                     results_nsim$MI_results$CMI_23_list[['97.5%']],
                                     results_nsim$MI_results$norm_MI_12_list[['97.5%']],
                                     results_nsim$MI_results$norm_MI_13_list[['97.5%']],
                                     results_nsim$MI_results$norm_MI_23_list[['97.5%']],
                                     results_nsim$MI_results$norm_CMI_12_list[['97.5%']],
                                     results_nsim$MI_results$norm_CMI_13_list[['97.5%']],
                                     results_nsim$MI_results$norm_CMI_23_list[['97.5%']])

        for (idx in 1:length(MI_params)){
            if (results_nsim$MI_results$MI_true_results[[idx]] == 0){
                tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 0)
                tmp_q925 = round(MI_q025_list_weak[[i]][[idx]], 0)
#                 tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 3)
#                 tmp_q925 = ifelse(MI_q025_list_weak[[i]][[idx]] < 0.1, 0, MI_q025_list_weak[[i]][[idx]]) 
                  # rounding error (needs to be different for original / standardized MI
            }else{
                tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 3) 
                tmp_q925 = round(MI_q025_list_weak[[i]][[idx]], 3)                
            }
            
            if (results_nsim$MI_results$MI_true_results[[idx]] <= tmp_q975 &
                results_nsim$MI_results$MI_true_results[[idx]] >= tmp_q925) {
                MI_coverage_95_weak[i, idx] = 1
            }
        }
    }
    
    MI_coverage_95_prob_weak = colMeans(MI_coverage_95_weak, na.rm=T)
    
    return(list('MI_true'= results_nsim$MI_results$MI_true_results, 
                'MI_coverage_95_weak'= MI_coverage_95_weak, 
                'MI_coverage_95_prob_weak'= MI_coverage_95_prob_weak, 
                'MI_median_list_weak'= MI_median_list_weak,
                'MI_q025_list_weak' = MI_q025_list_weak, 
                'MI_q975_list_weak' = MI_q975_list_weak))
}    

# #--------------------------------------------------------
# # calculate coverage probabilities for mutual information
# #--------------------------------------------------------
# coverage_MI_95 = function(dir_home, file_results, nsims){
#     MI_params = c('MI_all', 'MI_12', 'MI_13', 'MI_23', 'norm_MI_all', 'norm_MI_12', 'norm_MI_13', 'norm_MI_23')

#     MI_median_list_weak = MI_q025_list_weak = MI_q975_list_weak = MI_MI_95_weak = list()
#     MI_coverage_95_weak = matrix(rep(0, nsims*length(MI_params)), nrow=nsims)
#     colnames(MI_coverage_95_weak) = MI_params
#     error_nsim_weak = NULL

#     for (i in 1:nsims){
#         dir_results = paste0(dir_home, file_results, i, '.Rdata')
#         tmp = try(load(dir_results), TRUE)
#         if(isTRUE(class(tmp)=="try-error")) { 
#             error_nsim_weak = c(error_nsim_weak, i);
#             MI_coverage_95_weak[i,] = NA
#             next } 

#         MI_median_list_weak[[i]] = c(results_nsim$MI_results$MI_all_list[['50%']],
#                                      results_nsim$MI_results$MI_12_list[['50%']],
#                                      results_nsim$MI_results$MI_13_list[['50%']],
#                                      results_nsim$MI_results$MI_23_list[['50%']],
#                                      results_nsim$MI_results$norm_MI_all_list[['50%']],
#                                      results_nsim$MI_results$norm_MI_12_list[['50%']],
#                                      results_nsim$MI_results$norm_MI_13_list[['50%']],
#                                      results_nsim$MI_results$norm_MI_23_list[['50%']])
#         MI_q025_list_weak[[i]] = c(results_nsim$MI_results$MI_all_list[['2.5%']],
#                                      results_nsim$MI_results$MI_12_list[['2.5%']],
#                                      results_nsim$MI_results$MI_13_list[['2.5%']],
#                                      results_nsim$MI_results$MI_23_list[['2.5%']],
#                                      results_nsim$MI_results$norm_MI_all_list[['2.5%']],
#                                      results_nsim$MI_results$norm_MI_12_list[['2.5%']],
#                                      results_nsim$MI_results$norm_MI_13_list[['2.5%']],
#                                      results_nsim$MI_results$norm_MI_23_list[['2.5%']])
#         MI_q975_list_weak[[i]] = c(results_nsim$MI_results$MI_all_list[['97.5%']],
#                                      results_nsim$MI_results$MI_12_list[['97.5%']],
#                                      results_nsim$MI_results$MI_13_list[['97.5%']],
#                                      results_nsim$MI_results$MI_23_list[['97.5%']],
#                                      results_nsim$MI_results$norm_MI_all_list[['97.5%']],
#                                      results_nsim$MI_results$norm_MI_12_list[['97.5%']],
#                                      results_nsim$MI_results$norm_MI_13_list[['97.5%']],
#                                      results_nsim$MI_results$norm_MI_23_list[['97.5%']])

#         for (idx in 1:length(MI_params)){
#             if (results_nsim$MI_results$MI_true_results[[idx]] == 0){
#                 tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 0)
#                 tmp_q925 = round(MI_q025_list_weak[[i]][[idx]], 0)
# #                 tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 3)
# #                 tmp_q925 = ifelse(MI_q025_list_weak[[i]][[idx]] < 0.1, 0, MI_q025_list_weak[[i]][[idx]]) 
#                   # rounding error (needs to be different for original / standardized MI
#             }else{
#                 tmp_q975 = round(MI_q975_list_weak[[i]][[idx]], 3) 
#                 tmp_q925 = round(MI_q025_list_weak[[i]][[idx]], 3)                
#             }
            
#             if (results_nsim$MI_results$MI_true_results[[idx]] <= tmp_q975 &
#                 results_nsim$MI_results$MI_true_results[[idx]] >= tmp_q925) {
#                 MI_coverage_95_weak[i, idx] = 1
#             }
#         }
#     }
    
#     MI_coverage_95_prob_weak = colMeans(MI_coverage_95_weak, na.rm=T)
    
#     return(list('MI_true'= results_nsim$MI_results$MI_true_results, 
#                 'MI_coverage_95_weak'= MI_coverage_95_weak, 
#                 'MI_coverage_95_prob_weak'= MI_coverage_95_prob_weak, 
#                 'MI_median_list_weak'= MI_median_list_weak,
#                 'MI_q025_list_weak' = MI_q025_list_weak, 
#                 'MI_q975_list_weak' = MI_q975_list_weak))
# }    

#----------------------------------------------------------------
# generate summarizing table for mutual information estimation
#----------------------------------------------------------------
summary_MI <- function(dir_home, file_dir, R_ver){
    # nsims: total number of simulation replicates
    # R_ver: which correlation version
    
    load(paste0(dir_home, file_dir))
    
    if (R_ver == 0){
        results = results_R_v0_MI
    }else if (R_ver == 1){
        results = results_R_v1_MI
    }else if (R_ver == 2){
        results = results_R_v2_MI
    }else if (R_ver == 3){
        results = results_R_v3_MI
    }else if (R_ver == 4){
        results = results_R_v4_MI
    }
    
    MI_true = unlist(results$MI_true)
    MI_median = matrix(unlist(results$MI_median_list_weak), byrow=T, 
                   nrow=length(results$MI_median_list_weak))
    MI_q025 = matrix(unlist(results$MI_q025_list_weak), byrow=T, 
                   nrow=length(results$MI_q025_list_weak))
    MI_q975 = matrix(unlist(results$MI_q975_list_weak), byrow=T, 
                   nrow=length(results$MI_q975_list_weak))
    MI_cov_prob = unname(unlist(results$MI_coverage_95_prob_weak))

    stat_median = colMeans(MI_median)
    stat_q025 = colMeans(MI_q025)
    stat_q975 = colMeans(MI_q975)
    
    table_MI = as.data.frame(cbind(MI_true, stat_median, stat_q025, stat_q975, MI_cov_prob))
    colnames(table_MI) = c('Truth', 'Median', 'Q025', 'Q975', 'Cov_prob')
    table_MI = round(table_MI, 3)
    
    return(table_MI)
}