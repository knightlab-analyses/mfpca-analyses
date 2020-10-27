#-------------------------------------------
# results on optimal model
#-------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../../application_code/mfpca_apply_fcns.R')
#----------------------------------------------
# model diagnostics
#----------------------------------------------
mdl_idx = as.numeric(Sys.getenv("PBS_ARRAYID"))
load(paste0(dir, '/results/mfpca_stan_mdl_', mdl_idx, '.RData'))
load(paste0(dir, '/results/mfpca_looic_models_compare.RData'))
          
mdl_optimal_param = mdl_list[[paste0('mdl', mdl_idx)]]


diagnostic_looic(loo_list=looic.list[[mdl_idx]], num_subjects=pca_data$N, 
                 fig_dir='../figures/') 

diagnostic_posterior(mfpca_fit=fit, Nsamples=Nsamples, Nchains=Nchains, 
                     visits_matrix=pca_data$V, response_observed=pca_data$Y,
				     fig_title='posterior predictive checking', fig_dir='../figures/') 

#------------------------------------------
# results checking & visualization
#------------------------------------------
post_rotation_results <- post_hoc_rotation(fit = fit, Nchains = Nchains, Nsamples = Nsamples, 
									       N = pca_data$N, P = pca_data$P, K = pca_data$K, Q = pca_data$Q)


# with 1SD of FPC scores on FPC curves
mfpca_results <- visu_mfpca(dat_mfpca, param_stan, post_rotation_results, title = c('bacterial', 'molecular'),
						 x_label='Timepoint (weeks)', y_label ='Shannon', 
						 x_range = c(0,9), FPC_sd = TRUE, fig_dir='../figures/') 
                            #fig_dir=paste0(dir, '/figures/'))

# without 1SD of FPC scores on FPC curves
mfpca_results <- visu_mfpca(dat_mfpca, param_stan, post_rotation_results, title = c('bacterial', 'molecular'),
						 x_label='Timepoint (weeks)', y_label ='Shannon', 
						 x_range = c(0,9), FPC_sd = FALSE, fig_dir='../figures/') 
                            #fig_dir=paste0(dir, '/figures/'))
     
# compute mutual information
MI_results <- MI_est_data(post_rotation_results=post_rotation_results, K=pca_data$K)
     
#----------------------------------
# boxplots on FPC scores
#----------------------------------
table_scores = mfpca_results$FPC_scores

#fig_dir=paste0(dir, '/figures/')
fig_dir='../figures/'
ggplot(aes(y = bacterial_FPC1, x = BodySite, fill = Sex), data = table_scores) + geom_boxplot() + 
theme_classic() + theme(axis.text = element_text(color='black')) +
ylab('Bacterial Diversity FPC1') +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"),
      legend.position = 'none')  
ggsave(paste0(fig_dir, 'fpc1_scores_bacterial_SiteSex.pdf'), width = 4, height = 3, dpi=300)

ggplot(aes(y = molecular_FPC1, x = BodySite, fill = Sex), data = table_scores) + geom_boxplot() + 
theme_classic() + theme(axis.text = element_text(color='black')) +
ylab('Bacterial Diversity FPC1') +
theme(axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"),
      legend.position = 'none')   
ggsave(paste0(fig_dir, 'fpc1_scores_molecular_SiteSex.pdf'), width = 4, height = 3, dpi=300)

#----------------------------
# export results
#-----------------------------
save(mdl_optimal_param, post_rotation_results, 
     MI_results, table_scores, file=paste0(dir, '/results/mfpca_optimal_mdl_results.RData'))     

