############################################################
###### estimate all parameters ########
############################################################
args = commandArgs(trailingOnly=TRUE)
dir = args[1] # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'

i <- as.numeric(Sys.getenv("PBS_ARRAYID"))


#----------------------
#check results
#----------------------
source('mfpca_rotation_v2.R')
source('sim_settings_v2.R') # for plotting figures
source('basis_setup_sparse.R')

load(paste0(dir, '/sim_data_NB/sim_3B_file_names.Rdata'))

for (b in 1:length(file_names)){
    print(paste0('file_idx: ', b))
    
    file_path = file_names[[b]]
    file_name = substr(file_path, start=81, stop=1000)
    fig_name = substr(file_name,1,nchar(file_name)-6)

    #error_catch = try(load(paste0(dir, '/results_sim_NB/', i, '_', file_name)), TRUE) # check whether mfpca results exist
    #if (class(error_catch) == 'try-error'){
    
    load(paste0(dir, '/results_sim_NB/', i, '_', file_name)) # updated ACT code 3_1_xxx
    if (is.na(diff_time)){ 
        print('NA')
        results_nsim = paste0('file_idx: ', b, ' nsim: ', i, ' = NA')
    }else{
        results_nsim = output_sim(mfpca_output = paste0(dir, '/results_sim_NB/', i, '_', file_name),
                               dat_sim = file_path,
                               idx_sim = i, mean_lim = c(-15,15), FPC_lim = c(-5,5),
                               fig_name = paste0('../figures_supp_test/', fig_name, '_nsims', i, '.pdf'),
                               fig_on = FALSE)
    }    
    
    save(results_nsim, file=paste0(dir, '/summary_est_NB/summary_', fig_name, '_nsims', i, '.Rdata'))
}    
