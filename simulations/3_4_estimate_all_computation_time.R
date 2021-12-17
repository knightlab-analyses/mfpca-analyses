args = commandArgs(trailingOnly=TRUE)
dir = args[1]  # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'

nsim = as.numeric(Sys.getenv("PBS_ARRAYID"))

load(paste0(dir, '/sim_data_NB/sim_3B_file_names.Rdata'))

n_files = length(file_names)
time_matrix = matrix(0, nrow=n_files, ncol= nsim) 
# each row: a simulation scenario; each column: simulation replicate

summary_matrix = as.data.frame(matrix(0, nrow=length(file_names), ncol = 6))
colnames(summary_matrix) = c('num_mod', 'num_sample', 'num_time', 'prop_miss', 'mean_t', 'sd_t')

for (b in 1:n_files){ # problem start after file_name[23]
    print(paste0('file_idx: ', b))
    file = file_names[[b]]
    
    # simulation scenario
    txt = substr(file, start=81, stop=1000)
    txt = substr(txt,1,nchar(txt)-6) # drop ".Rdata"
    case = read.table(text = as.character(txt), sep = "_")
    
    num_mod = case$V4
    num_sample = case$V9
    num_time = case$V12
    prop_miss = case$V14
    
    summary_matrix[b, 1:4] = c(num_mod, num_sample, num_time, prop_miss)

    for (i in 1:nsim){
        #print(paste0('nsim: ', i))
        result = paste0(dir, '/results_sim_NB/', i, '_', substr(file, start=81, stop=1000))
        error_catch = try(load(result), TRUE)
        if (class(error_catch) == 'try-error'){
            time_matrix[b, i] = NA
        }else{
            time_matrix[b, i] = round(as.numeric(diff_time, units = "mins"),3)
        }
        
    }
    
    summary_matrix[b, 5] = round(mean(time_matrix[b, ], na.rm=T),2)
    summary_matrix[b, 6] = round(sd(time_matrix[b, ], na.rm=T),2)

}    

save(summary_matrix, time_matrix, 
     file=paste0(dir, '/computation_time/time_est_all_3B_totalSims', nsim, '.Rdata'))

