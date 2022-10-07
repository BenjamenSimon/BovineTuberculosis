
include("LoadPackages.jl")
include("Likelihood.jl")
include("DataUpdaters.jl")
include("Proposals.jl")
include("Posterior.jl")
include("MCMCfuncs.jl")

#################
### LOAD DATA ###
#################

combi_array = load("Data/Set 1/combi_array.jld2")["array"]

f_to_p_dict = load("Data/Set 1/f_to_p_dict.jld2")["dict"]

record_of_movements = load("Data/Set 1/record_of_movements_oi.jld2")["array"]

dict_of_movements = load("Data/Set 1/dict_of_movements.jld2")["dict"]

ids_to_pos_dict = load("Data/Set 1/ids_to_pos_dict.jld2")["dict"]

###########
### RUN ###
###########

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015

F_tr = 0.004
ϵ_tr = 0.05

ρ_tr = 0.75
ρ_E_tr = 0.2

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52

epi_params_true = [β_c_tr, β_b_tr, γ_tr, F_tr, ϵ_tr, ρ_tr, ρ_E_tr, θ_bb_tr, θ_bd_tr]

# d_β_c = Uniform(0, 0.04)
# d_β_b = Uniform(0, 0.08)
# d_γ = Uniform(0, 0.3)
# d_F = Uniform(0, 0.08)
# d_ϵ = Uniform(0, 0.1)
# d_ρ = Uniform(0, 1)
# d_ρ_E = Uniform(0, 1)

d_β_c = Gamma(11, 0.0002)
d_β_b = Gamma(21, 0.0002)
d_γ = Gamma(21, 0.00075)
d_F = Gamma(21, 0.0002)
d_ϵ = Gamma(11, 0.005)
d_ρ = Beta(2.5, 1.5)
d_ρ_E = Beta(2, 5)


epi_params_dists = [d_β_c, d_β_b, d_γ, d_F, d_ϵ, d_ρ, d_ρ_E]


############################
### Testing the Updaters ###
############################

r1, or1, ar1, tr1, ut1 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                          combi_array = combi_array, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                          dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                          ids_to_pos_dict = ids_to_pos_dict)

r1_cut = r1[1:5499, :]
or1_cut = or1[1:5499, :]
ar1_cut = ar1[1:5499, :]
tr1_cut = tr1[1:5499, :]
ut1_cut = ut1[1:5499, :]

r2, or2, ar2, tr2, ut2 = Blk_Adaptive_RWM_MCMC(;N_its = 55000, infer_block = [false, true], data_aug_infer = [false, false, false, false, false, false, false, false],
                          combi_array = combi_array, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 0.03, 0.2, 0.04],
                          dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                          ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [true, false, false, false, false, false, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [false, true, false, false, false, false, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [false, false, true, false, false, false, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [false, false, false, true, false, false, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [false, false, false, false, true, false, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, false], data_aug_infer = [false, false, false, false, false, true, false, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [false, false], data_aug_infer = [false, false, false, false, false, false, true, false],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)



Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [false, false], data_aug_infer = [false, false, false, false, false, false, false, true],
                        combi_array = combi_array, moves_record = record_of_movements,
                        params_init = epi_params_true, tuning = [0.01, 0.02, 0.02, 0.03],
                        dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                        ids_to_pos_dict = ids_to_pos_dict)

#################################
### Benchmarking the Updaters ###
#################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

@benchmark begin
  update_data_Move_SE(combi_array, 1, 26, 50, 1)
end



##########################
### Profiling function ###
##########################

using ProfileView

@profview Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                            combi_array = combi_array, moves_record = record_of_movements,
                            params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                            dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                            ids_to_pos_dict = ids_to_pos_dict)

warntype_last()

ProfileView.view(nothing)

##############################
### New profiling function ###
##############################

using Profile
using PProf

# collect a profile
@profile Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                          combi_array = combi_array, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                          dict_of_movements = dict_of_movements, f_to_p_dict = f_to_p_dict,
                          ids_to_pos_dict = ids_to_pos_dict)

# Export pprof profile and open interactive profiling web interface.
pprof()

# https://github.com/JuliaPerf/PProf.jl
# https://github.com/google/pprof/blob/main/doc/README.md

PProf.refresh(file="PProf_MCMC_inf.pb.gz")
