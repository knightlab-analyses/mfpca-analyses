#-------------------------------------------
# select optimal model + diagnostics
#-------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('../../mfpca_apply_fcns.R')
                        
#---------------------------------------------------------------------------------------------------------------------------
# model selection 
# ref: https://discourse.mc-stan.org/t/understanding-looic/13409/9
# Rather than looking at LOOIC, you should look at the elpd (expected log predictive density), as LOOIC is just -2 * elpd
# For elpd, the larger the better; for looic: lower is better
# elpd_loo: if elpd_diff > 6 * se_diff, then confident that simpler model is better (at least 4-5 times)
#-------------------------------------------------------------------------------------------------------------------------
num_mdl = as.numeric(Sys.getenv("PBS_ARRAYID"))
looic.list = list()
for (i in 1:num_mdl){
    print(paste0('i:', i))
    load(paste0(dir, '/results/mfpca_t2d_A1c_PDGFBB_IGHA2_ratio_v2_', i, '.RData'))
    looic.list[[i]] = loo(extract(fit, "log_lik")[[1]])
}
looic.compare = loo::loo_compare(looic.list)
looic.compare = as.data.frame(looic.compare)[, c('elpd_diff', 'se_diff')]
save(looic.list, looic.compare, file=paste0(dir, '/results/mfpca_looic_models_A1c_PDGFBB_IGHA2_ratio_compare_v2.RData'))
                       


