args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1] # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations/sim_data_NB'

source("basis_setup_sparse.R")
source("sim_settings.R") # for simulation parameters
source("simulate_data_sparse.R") # for generate simulated data

# #------------------------------------
# # simulation setting
# # data modalities = 3
# # number of samples = 50, 100, 500
# # number of time pts = 5, 10, 20
# # missingness: 20%, 50%, 80%
# #------------------------------------
# N_range = c(50, 100, 500)
# N_T_max_range = c(5, 10, 20) # time pts
# missing_range = c(0.2, 0.5, 0.8) # % missing

# sim_idx = 0
# file_names = list()
# for (N in N_range){
#     for (N_T_max in N_T_max_range){
#         for (missing_prop in missing_range){
#             sim_idx = sim_idx + 1
#             print(paste0('sim_idx: ', sim_idx, ' #samples: ', N, ' # time pts: ', N_T_max, ' missing %: ', missing_prop))
            
#             N_T_mean = (1-missing_prop)*N_T_max
#             #print(paste0('N_T_mean: ', N_T_mean))      
#             sim_param = sim_setting_3B_T1D(N=N, P=3, N_T_max=N_T_max, N_T_mean=N_T_mean, R_corr=3)
            
#             # generate data
#             TIME_SPARSE=list()
#             Y_SPARSE=list()
#             MU_SPARSE=list()
#             F_SPARSE=list()
#             OMEGA_SPARSE=list()
#             ALPHA=list()
#             phi_t=list()
#             nsims = 1000  # number of replicates for each simulated scenario: needed for coverage probability calculation

#             for(ii in 1:nsims){
#                 print(paste0('nsims: ', ii))
#                 simdata=simulate_data_sparse(T_range=sim_param$T_range, N_T_max=sim_param$N_T_max, N_T_mean=sim_param$N_T_mean,
#                                              P=sim_param$P, N=sim_param$N, params=sim_param$params, orth=sim_param$orth,
#                                              nseeds=ii)
#                 TIME_SPARSE[[ii]]=simdata[[1]]
#                 Y_SPARSE[[ii]]=simdata[[2]]
#                 MU_SPARSE[[ii]]=simdata[[3]]
#                 F_SPARSE[[ii]]=simdata[[4]]
#                 OMEGA_SPARSE[[ii]]=simdata[[5]]
#                 ALPHA[[ii]]=simdata[[6]]
#                 phi_t[[ii]] = simdata[[9]]
#             }

#             file_names[[sim_idx]] = paste0(dir, '/', paste(paste('sim_dat_B', sim_param$P, 'R_corr', sim_param$R_corr,
#                                                    'N', sim_param$N, 'time_points', N_T_max, 
#                                                    'missing', missing_prop, sep='_'), '.Rdata', sep=''))
            
#             save(nsims, TIME_SPARSE, Y_SPARSE, MU_SPARSE, F_SPARSE, OMEGA_SPARSE, ALPHA, phi_t, sim_param, 
#                  file=file_names[[sim_idx]])
#         }
#     }    
# }
# save(file_names, file = paste0(dir, '/', 'sim_3B_file_names.Rdata', sep=''))


#-----------------------------------------------------
# simulate missing data for N = 500 selected cases
#-----------------------------------------------------
N_range = c(500)
N_T_max_range = c(10, 20) # time pts
missing_range = c(0.2, 0.5, 0.8) # % missing

sim_idx = 21 
#file_names = list() load file_names from existing file
load('/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations/sim_data_NB/sim_3B_file_names.Rdata') 
for (N in N_range){
    for (N_T_max in N_T_max_range){
        for (missing_prop in missing_range){
            sim_idx = sim_idx + 1
            print(paste0('sim_idx: ', sim_idx, ' #samples: ', N, ' # time pts: ', N_T_max, ' missing %: ', missing_prop))
            
            N_T_mean = (1-missing_prop)*N_T_max
            #print(paste0('N_T_mean: ', N_T_mean))      
            sim_param = sim_setting_3B_T1D(N=N, P=3, N_T_max=N_T_max, N_T_mean=N_T_mean, R_corr=3)
            
            # generate data
            TIME_SPARSE=list()
            Y_SPARSE=list()
            MU_SPARSE=list()
            F_SPARSE=list()
            OMEGA_SPARSE=list()
            ALPHA=list()
            phi_t=list()
            nsims = 1000  # number of replicates for each simulated scenario: needed for coverage probability calculation

            for(ii in 1:nsims){
                print(paste0('nsims: ', ii))
                simdata=simulate_data_sparse(T_range=sim_param$T_range, N_T_max=sim_param$N_T_max, N_T_mean=sim_param$N_T_mean,
                                             P=sim_param$P, N=sim_param$N, params=sim_param$params, orth=sim_param$orth,
                                             nseeds=ii)
                TIME_SPARSE[[ii]]=simdata[[1]]
                Y_SPARSE[[ii]]=simdata[[2]]
                MU_SPARSE[[ii]]=simdata[[3]]
                F_SPARSE[[ii]]=simdata[[4]]
                OMEGA_SPARSE[[ii]]=simdata[[5]]
                ALPHA[[ii]]=simdata[[6]]
                phi_t[[ii]] = simdata[[9]]
            }

            file_output = paste0(dir, '/', paste(paste('sim_dat_B', sim_param$P, 'R_corr', sim_param$R_corr,
                                                   'N', sim_param$N, 'time_points', N_T_max, 
                                                   'missing', missing_prop, sep='_'), '.Rdata', sep=''))
            
            save(nsims, TIME_SPARSE, Y_SPARSE, MU_SPARSE, F_SPARSE, OMEGA_SPARSE, ALPHA, phi_t, sim_param, 
                 file=file_output)
        }
    }    
}



