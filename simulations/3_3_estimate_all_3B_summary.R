source('mfpca_sim_summary.R')

args = commandArgs(trailingOnly=TRUE)
dir = args[1] # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'

load(paste0(dir, '/sim_data_NB/sim_3B_file_names.Rdata'))

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))

n_files = length(file_names)
cov_prob_table = as.data.frame(matrix(0, nrow=n_files, ncol=5*5)) # 25 covariance parameters (R_corr = 3)
colnames(cov_prob_table) = c('cov_alpha[1,1]', 'cov_alpha[1,2]', 'cov_alpha[1,3]', 'cov_alpha[1,4]', 'cov_alpha[1,5]',
                             'cov_alpha[2,1]', 'cov_alpha[2,2]', 'cov_alpha[2,3]', 'cov_alpha[2,4]', 'cov_alpha[2,5]',
                             'cov_alpha[3,1]', 'cov_alpha[3,2]', 'cov_alpha[3,3]', 'cov_alpha[3,4]', 'cov_alpha[3,5]',
                             'cov_alpha[4,1]', 'cov_alpha[4,2]', 'cov_alpha[4,3]', 'cov_alpha[4,4]', 'cov_alpha[4,5]', 
                             'cov_alpha[5,1]', 'cov_alpha[5,2]', 'cov_alpha[5,3]', 'cov_alpha[5,4]', 'cov_alpha[5,5]')
for (b in 1:n_files){
    print(paste0('file_idx: ', b))
    file_path = file_names[[b]]
    file_name = substr(file_path, start=81, stop=1000)
    fig_name = substr(file_name,1,nchar(file_name)-6)
    results_path = paste0('/summary_est_NB/summary_', fig_name, '_nsims')
    
    load(paste0(dir, '/results_sim_NB/', i, '_', file_name))
    if (is.na(diff_time)){
        cov_prob_table[b, ] = NA
    }else{
        results_sim_supp_cov = coverage_cov_95(dir_home = dir, 
                               file_data = paste0('/sim_data_NB/', file_name), 
                               file_results = results_path, 
                               nsims=i) # total number of simulation replicates
        cov_prob_table[b, ] = unname(unlist(results_sim_supp_cov$cov_coverage_95_prob_weak))
    }
    print(cov_prob_table[b,])
}    

save(cov_prob_table, 
     file=paste0(dir, '/cov_prob_est_NB/cov_prob_est_all_3B_totalSims', i, '.Rdata'))

