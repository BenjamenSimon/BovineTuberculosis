
include("LoadPackages.jl")
include("Likelihood.jl")
include("DataUpdaters.jl")
include("Proposals.jl")
include("Posterior.jl")
include("MCMCfuncs.jl")

#################
### LOAD DATA ###
#################

# combi_array = load("Data/Set 1/combi_array.jld2")["array"]

# combi_array = load("Data/Set 1/combi_array_unnamed.jld2")["array"]

DATA_res_and_track = load("Data/Set 1/DATA_res_and_track.jld2")["array"]
DATA_pers_and_parish = load("Data/Set 1/DATA_pers_and_parish.jld2")["array"]

# f_to_p_dict = load("Data/Set 1/f_to_p_dict.jld2")["dict"]

f_to_p_structs = load("Data/Set 1/f_to_p_structs.jld2")["struct"]

record_of_movements = Array(load("Data/Set 1/record_of_movements_oi.jld2")["array"])

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

r1, or1, ar1, tr1, ut1 = Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 1/res_infpars.csv", r1, header = true)
CSV.write("Inference/Test 1/other_res_infpars.csv", or1, header = true)
CSV.write("Inference/Test 1/aug_res_infpars.csv", ar1, header = true)
CSV.write("Inference/Test 1/tuning_res_infpars.csv", tr1, header = true)
CSV.write("Inference/Test 1/update_tracker_infpars.csv", ut1, header = true)

r2, or2, ar2, tr2, ut2 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [false, true], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 2/res_detpars.csv", r2, header = true)
CSV.write("Inference/Test 2/other_res_detpars.csv", or2, header = true)
CSV.write("Inference/Test 2/aug_res_detpars.csv", ar2, header = true)
CSV.write("Inference/Test 2/tuning_res_detpars.csv", tr2, header = true)
CSV.write("Inference/Test 2/update_tracker_detpars.csv", ut2, header = true)

r3, or3, ar3, tr3, ut3 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, true], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 3/res_allpars.csv", r3, header = true)
CSV.write("Inference/Test 3/other_res_allpars.csv", or3, header = true)
CSV.write("Inference/Test 3/aug_res_allpars.csv", ar3, header = true)
CSV.write("Inference/Test 3/tuning_res_allpars.csv", tr3, header = true)
CSV.write("Inference/Test 3/update_tracker_allpars.csv", ut3, header = true)

r4, or4, ar4, tr4, ut4 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [true, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 4/res_infpars_mSE.csv", r4, header = true)
CSV.write("Inference/Test 4/other_res_infpars_mSE.csv", or4, header = true)
CSV.write("Inference/Test 4/aug_res_infpars_mSE.csv", ar4, header = true)
CSV.write("Inference/Test 4/tuning_res_infpars_mSE.csv", tr4, header = true)
CSV.write("Inference/Test 4/update_tracker_infpars_mSE.csv", ut4, header = true)


r5, or5, ar5, tr5, ut5 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, true, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 5/res_infpars_mEI.csv", r5, header = true)
CSV.write("Inference/Test 5/other_res_infpars_mEI.csv", or5, header = true)
CSV.write("Inference/Test 5/aug_res_infpars_mEI.csv", ar5, header = true)
CSV.write("Inference/Test 5/tuning_res_infpars_mEI.csv", tr5, header = true)
CSV.write("Inference/Test 5/update_tracker_infpars_mEI.csv", ut5, header = true)


r6, or6, ar6, tr6, ut6 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, true, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 6/res_infpars_mEI.csv", r6, header = true)
CSV.write("Inference/Test 6/other_res_infpars_mEI.csv", or6, header = true)
CSV.write("Inference/Test 6/aug_res_infpars_mEI.csv", ar6, header = true)
CSV.write("Inference/Test 6/tuning_res_infpars_mEI.csv", tr6, header = true)
CSV.write("Inference/Test 6/update_tracker_infpars_mEI.csv", ut6, header = true)


r7, or7, ar7, tr7, ut7 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 7/res_infpars_mEI.csv", r7, header = true)
CSV.write("Inference/Test 7/other_res_infpars_mEI.csv", or7, header = true)
CSV.write("Inference/Test 7/aug_res_infpars_mEI.csv", ar7, header = true)
CSV.write("Inference/Test 7/tuning_res_infpars_mEI.csv", tr7, header = true)
CSV.write("Inference/Test 7/update_tracker_infpars_mEI.csv", ut7, header = true)



r12, or12, ar12, tr12, ut12 = Blk_Adaptive_RWM_MCMC(;N_its = 10000, infer_block = [false, true], data_aug_infer = [false, false, false, false, true, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 12/res_infpars_mEI.csv", r12, header = true)
CSV.write("Inference/Test 12/other_res_infpars_mEI.csv", or12, header = true)
CSV.write("Inference/Test 12/aug_res_infpars_mEI.csv", ar21, header = true)
CSV.write("Inference/Test 12/tuning_res_infpars_mEI.csv", tr12, header = true)
CSV.write("Inference/Test 12/update_tracker_infpars_mEI.csv", ut12, header = true)





r15, or15, ar15, tr15, ut15 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, true],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test 15/res_infpars_mEI.csv", r15, header = true)
CSV.write("Inference/Test 15/other_res_infpars_mEI.csv", or15, header = true)
CSV.write("Inference/Test 15/aug_res_infpars_mEI.csv", ar15, header = true)
CSV.write("Inference/Test 15/tuning_res_infpars_mEI.csv", tr15, header = true)
CSV.write("Inference/Test 15/update_tracker_infpars_mEI.csv", ut15, header = true)











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

Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                            DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                            params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                            dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                            ids_to_pos_dict = ids_to_pos_dict)

using ProfileView

@profview Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                            DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                            params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                            dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                            ids_to_pos_dict = ids_to_pos_dict)

warntype_last()

@code_warntype s(x)

# Articles on speeding up julia code
https://discourse.julialang.org/t/code-warntype-shows-lots-of-union-nothing-tuple-int64-int64/60748/4

ProfileView.view(nothing)

##############################
### New profiling function ###
##############################

using Profile
using PProf

# collect a profile
@profile Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                            DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                            params_init = epi_params_true, tuning = [0.05, 45, 0.1, 15],
                            dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                            ids_to_pos_dict = ids_to_pos_dict)

# Export pprof profile and open interactive profiling web interface.
pprof()

# https://github.com/JuliaPerf/PProf.jl
# https://github.com/google/pprof/blob/main/doc/README.md

PProf.refresh(file="PProf_MCMC_inf_fast.pb.gz")
