args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1] # dir = '/panfs/panfs1.ucsd.edu/panscratch/lij014/Oct_2020_mfpca/simulations'

source("basis_setup_sparse.R")
source("sim_settings_highMod.R") # for simulation parameters
source("simulate_data_sparse.R") # for generate simulated data

##############################
##############################
##  SIMULATION PARAMETERS   ##
##############################
##############################
# 100 subjects, 10 time points, 20% missing, 10 modalities
sim_param = sim_setting_10B(N=100, P=10, N_T_max=10, N_T_mean=8, R_corr=3) 


# #######################
# #######################
# ##   GENERATE DATA   ##
# #######################
# #######################
nsims=100

TIME_SPARSE=list()
Y_SPARSE=list()
MU_SPARSE=list()
F_SPARSE=list()
OMEGA_SPARSE=list()
ALPHA=list()
phi_t=list()

for(ii in 1:nsims){
    print(paste0('nsims:', ii))
    simdata=simulate_data_sparse(T_range=sim_param$T_range, N_T_max=sim_param$N_T_max, N_T_mean=sim_param$N_T_mean,
    	                         P=sim_param$P, N=sim_param$N, params=sim_param$params, orth=sim_param$orth, nseeds=ii)
	TIME_SPARSE[[ii]]=simdata[[1]]
	Y_SPARSE[[ii]]=simdata[[2]]
	MU_SPARSE[[ii]]=simdata[[3]]
	F_SPARSE[[ii]]=simdata[[4]]
	OMEGA_SPARSE[[ii]]=simdata[[5]]
	ALPHA[[ii]]=simdata[[6]]
	phi_t[[ii]] = simdata[[9]]
}

save(nsims, TIME_SPARSE, Y_SPARSE, MU_SPARSE, F_SPARSE, OMEGA_SPARSE, ALPHA, phi_t, sim_param, 
	 file=paste0(dir, '/sim_data_10B/', paste(paste('sim_dat_5B', 'N', sim_param$N, 'R_corr', sim_param$R_corr, 'Complete', 
	                  sim_param$complete, 'nsims', nsims, sep='_'), '.Rdata', sep=''))) 









