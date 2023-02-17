
include("LoadPackages.jl")
include("Likelihood.jl")
include("DataUpdaters.jl")
include("Proposals.jl")
include("Posterior.jl")
include("MCMCfuncs.jl")

#################
### LOAD DATA ###
#################

DATA_res_and_track = load("Data/Set 1/DATA_res_and_track.jld2")["array"]
DATA_pers_and_parish = load("Data/Set 1/DATA_pers_and_parish.jld2")["array"]

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


#################################################
#### 1. Infection parameters - Strong Priors ####
#################################################

r1, or1, ar1, tr1, ut1 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_1_[Inf]/res_1_[Inf].csv", r1, header = true)
CSV.write("Inference/Results_1_[Inf]/other_res_1_[Inf].csv", or1, header = true)
CSV.write("Inference/Results_1_[Inf]/aug_res_1_[Inf].csv", ar1, header = true)
CSV.write("Inference/Results_1_[Inf]/tuning_res_1_[Inf].csv", tr1, header = true)
CSV.write("Inference/Results_1_[Inf]/update_tracker_1_[Inf].csv", ut1, header = true)

#################################################
#### 2. Detection parameters - Strong Priors ####
#################################################

r2, or2, ar2, tr2, ut2 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [false, true], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_2_[Det]/res_2_[Det].csv", r2, header = true)
CSV.write("Inference/Results_2_[Det]/other_res_2_[Det].csv", or2, header = true)
CSV.write("Inference/Results_2_[Det]/aug_res_2_[Det].csv", ar2, header = true)
CSV.write("Inference/Results_2_[Det]/tuning_res_2_[Det].csv", tr2, header = true)
CSV.write("Inference/Results_2_[Det]/update_tracker_2_[Det].csv", ut2, header = true)

###########################################
#### 3. All parameters - Strong Priors ####
###########################################

r3, or3, ar3, tr3, ut3 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, true], data_aug_infer = [false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_3_[Inf][Det]/res_3_[Inf][Det].csv", r3, header = true)
CSV.write("Inference/Results_3_[Inf][Det]/other_res_3_[Inf][Det].csv", or3, header = true)
CSV.write("Inference/Results_3_[Inf][Det]/aug_res_3_[Inf][Det].csv", ar3, header = true)
CSV.write("Inference/Results_3_[Inf][Det]/tuning_res_3_[Inf][Det].csv", tr3, header = true)
CSV.write("Inference/Results_3_[Inf][Det]/update_tracker_3_[Inf][Det].csv", ut3, header = true)

############################################################################
#### 4. Infection parameters - Strong Priors - MoveSE Data Augmentation ####
############################################################################

r4, or4, ar4, tr4, ut4 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [true, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_4_[Inf][mSE]/res_4_[Inf][mSE].csv", r4, header = true)
CSV.write("Inference/Results_4_[Inf][mSE]/other_res_4_[Inf][mSE].csv", or4, header = true)
CSV.write("Inference/Results_4_[Inf][mSE]/aug_res_4_[Inf][mSE].csv", ar4, header = true)
CSV.write("Inference/Results_4_[Inf][mSE]/tuning_res_4_[Inf][mSE].csv", tr4, header = true)
CSV.write("Inference/Results_4_[Inf][mSE]/update_tracker_4_[Inf][mSE].csv", ut4, header = true)

############################################################################
#### 5. Infection parameters - Strong Priors - MoveEI Data Augmentation ####
############################################################################

r5, or5, ar5, tr5, ut5 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, true, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_5_[Inf][mEI]/res_5_[Inf][mEI].csv", r5, header = true)
CSV.write("Inference/Results_5_[Inf][mEI]/other_res_5_[Inf][mEI].csv", or5, header = true)
CSV.write("Inference/Results_5_[Inf][mEI]/aug_res_5_[Inf][mEI].csv", ar5, header = true)
CSV.write("Inference/Results_5_[Inf][mEI]/tuning_res_5_[Inf][mEI].csv", tr5, header = true)
CSV.write("Inference/Results_5_[Inf][mEI]/update_tracker_5_[Inf][mEI].csv", ut5, header = true)

###################################################################################
#### 6. Infection parameters - Strong Priors - Add/Remove SE Data Augmentation ####
###################################################################################

r6, or6, ar6, tr6, ut6 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, true, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 30, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_6_[Inf][arSE]/res_6_[Inf][arSE].csv", r6, header = true)
CSV.write("Inference/Results_6_[Inf][arSE]/other_res_6_[Inf][arSE].csv", or6, header = true)
CSV.write("Inference/Results_6_[Inf][arSE]/aug_res_6_[Inf][arSE].csv", ar6, header = true)
CSV.write("Inference/Results_6_[Inf][arSE]/tuning_res_6_[Inf][arSE].csv", tr6, header = true)
CSV.write("Inference/Results_6_[Inf][arSE]/update_tracker_6_[Inf][arSE].csv", ut6, header = true)

###################################################################################
#### 7. Infection parameters - Strong Priors - Add/Remove EI Data Augmentation ####
###################################################################################

r7, or7, ar7, tr7, ut7 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, false, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_7_[Inf][arEI]/res_7_[Inf][arEI].csv", r7, header = true)
CSV.write("Inference/Results_7_[Inf][arEI]/other_res_7_[Inf][arEI].csv", or7, header = true)
CSV.write("Inference/Results_7_[Inf][arEI]/aug_res_7_[Inf][arEI].csv", ar7, header = true)
CSV.write("Inference/Results_7_[Inf][arEI]/tuning_res_7_[Inf][arEI].csv", tr7, header = true)
CSV.write("Inference/Results_7_[Inf][arEI]/update_tracker_7_[Inf][arEI].csv", ut7, header = true)

#############################################################################
#### 8. Infection parameters - Strong Priors - mSE mEI Data Augmentation ####
#############################################################################

r8, or8, ar8, tr8, ut8 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [true, true, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_8_[Inf][mSE][mEI]/res_8_[Inf][mSE][mEI].csv", r8, header = true)
CSV.write("Inference/Results_8_[Inf][mSE][mEI]/other_res_8_[Inf][mSE][mEI].csv", or8, header = true)
CSV.write("Inference/Results_8_[Inf][mSE][mEI]/aug_res_8_[Inf][mSE][mEI].csv", ar8, header = true)
CSV.write("Inference/Results_8_[Inf][mSE][mEI]/tuning_res_8_[Inf][mSE][mEI].csv", tr8, header = true)
CSV.write("Inference/Results_8_[Inf][mSE][mEI]/update_tracker_8_[Inf][mSE][mEI].csv", ut8, header = true)

###############################################################################
#### 9. Infection parameters - Strong Priors - arSE arEI Data Augmentation ####
###############################################################################

r9, or9, ar9, tr9, ut9 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, true, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.02, 30, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_9_[Inf][arSE][arEI]/res_9_[Inf][arSE][arEI].csv", r9, header = true)
CSV.write("Inference/Results_9_[Inf][arSE][arEI]/other_res_9_[Inf][arSE][arEI].csv", or9, header = true)
CSV.write("Inference/Results_9_[Inf][arSE][arEI]/aug_res_9_[Inf][arSE][arEI].csv", ar9, header = true)
CSV.write("Inference/Results_9_[Inf][arSE][arEI]/tuning_res_9_[Inf][arSE][arEI].csv", tr9, header = true)
CSV.write("Inference/Results_9_[Inf][arSE][arEI]/update_tracker_9_[Inf][arSE][arEI].csv", ut9, header = true)

###############################################################################
#### 10. Infection parameters - Strong Priors - mSE arSE Data Augmentation ####
###############################################################################

r10, or10, ar10, tr10, ut10 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [true, false, true, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 30, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_10_[Inf][mSE][arSE]/res_10_[Inf][mSE][arSE].csv", r10, header = true)
CSV.write("Inference/Results_10_[Inf][mSE][arSE]/other_res_10_[Inf][mSE][arSE].csv", or10, header = true)
CSV.write("Inference/Results_10_[Inf][mSE][arSE]/aug_res_10_[Inf][mSE][arSE].csv", ar10, header = true)
CSV.write("Inference/Results_10_[Inf][mSE][arSE]/tuning_res_10_[Inf][mSE][arSE].csv", tr10, header = true)
CSV.write("Inference/Results_10_[Inf][mSE][arSE]/update_tracker_10_[Inf][mSE][arSE].csv", ut10, header = true)

###############################################################################
#### 11. Infection parameters - Strong Priors - mEI arEI Data Augmentation ####
###############################################################################

r11, or11, ar11, tr11, ut11 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, true, false, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_11_[Inf][mEI][arEI]/res_11_[Inf][mEI][arEI].csv", r11, header = true)
CSV.write("Inference/Results_11_[Inf][mEI][arEI]/other_res_11_[Inf][mEI][arEI].csv", or11, header = true)
CSV.write("Inference/Results_11_[Inf][mEI][arEI]/aug_res_11_[Inf][mEI][arEI].csv", ar11, header = true)
CSV.write("Inference/Results_11_[Inf][mEI][arEI]/tuning_res_11_[Inf][mEI][arEI].csv", tr11, header = true)
CSV.write("Inference/Results_11_[Inf][mEI][arEI]/update_tracker_11_[Inf][mEI][arEI].csv", ut11, header = true)

############################################################################
#### 12. Infection parameters - Strong Priors - arDet Data Augmentation ####
############################################################################

r12, or12, ar12, tr12, ut12 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [false, true], data_aug_infer = [false, false, false, false, true, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.5, 4],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_12_[Det][arDet]/res_12_[Det][arDet].csv", r12, header = true)
CSV.write("Inference/Results_12_[Det][arDet]/other_res_12_[Det][arDet].csv", or12, header = true)
CSV.write("Inference/Results_12_[Det][arDet]/aug_res_12_[Det][arDet].csv", ar12, header = true)
CSV.write("Inference/Results_12_[Det][arDet]/tuning_res_12_[Det][arDet].csv", tr12, header = true)
CSV.write("Inference/Results_12_[Det][arDet]/update_tracker_12_[Det][arDet].csv", ut12, header = true)

##############################################################################
#### 13. Infection parameters - Strong Priors - arDeath Data Augmentation ####
##############################################################################

r13, or13, ar13, tr13, ut13 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, true, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_13_[Inf][arDeath]/res_13_[Inf][arDeath].csv", r13, header = true)
CSV.write("Inference/Results_13_[Inf][arDeath]/other_res_13_[Inf][arDeath].csv", or13, header = true)
CSV.write("Inference/Results_13_[Inf][arDeath]/aug_res_13_[Inf][arDeath].csv", ar13, header = true)
CSV.write("Inference/Results_13_[Inf][arDeath]/tuning_res_13_[Inf][arDeath].csv", tr13, header = true)
CSV.write("Inference/Results_13_[Inf][arDeath]/update_tracker_13_[Inf][arDeath].csv", ut13, header = true)

#############################################################################
#### 14. Infection parameters - Strong Priors - arPEnv Data Augmentation ####
#############################################################################

r14, or14, ar14, tr14, ut14 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, true, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_14_[Inf][arPEnv]/res_14_[Inf][arPEnv].csv", r14, header = true)
CSV.write("Inference/Results_14_[Inf][arPEnv]/other_res_14_[Inf][arPEnv].csv", or14, header = true)
CSV.write("Inference/Results_14_[Inf][arPEnv]/aug_res_14_[Inf][arPEnv].csv", ar14, header = true)
CSV.write("Inference/Results_14_[Inf][arPEnv]/tuning_res_14_[Inf][arPEnv].csv", tr14, header = true)
CSV.write("Inference/Results_14_[Inf][arPEnv]/update_tracker_14_[Inf][arPEnv].csv", ut14, header = true)

##############################################################################
#### 15. Infection parameters - Strong Priors - arMoves Data Augmentation ####
##############################################################################

r15, or15, ar15, tr15, ut15 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, false, false, false, false, false, true],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.03, 50, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Results_15_[Inf][arMoves]/res_15_[Inf][arMoves].csv", r15, header = true)
CSV.write("Inference/Results_15_[Inf][arMoves]/other_res_15_[Inf][arMoves].csv", or15, header = true)
CSV.write("Inference/Results_15_[Inf][arMoves]/aug_res_15_[Inf][arMoves].csv", ar15, header = true)
CSV.write("Inference/Results_15_[Inf][arMoves]/tuning_res_15_[Inf][arMoves].csv", tr15, header = true)
CSV.write("Inference/Results_15_[Inf][arMoves]/update_tracker_15_[Inf][arMoves].csv", ut15, header = true)







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





### Testing for bug ###

r_bug, or_bug, ar_bug, tr_bug, ut_bug = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, true, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test Bug/res_infpars_arSEarEI.csv", r_bug, header = true)
CSV.write("Inference/Test Bug/other_res_infpars_arSEarEI.csv", or_bug, header = true)
CSV.write("Inference/Test Bug/aug_res_infpars_arSEarEI.csv", ar_bug, header = true)
CSV.write("Inference/Test Bug/tuning_res_infpars_arSEarEI.csv", tr_bug, header = true)
CSV.write("Inference/Test Bug/update_tracker_infpars_arSEarEI.csv", ut_bug, header = true)


r_bug2, or_bug2, ar_bug2, tr_bug2, ut_bug2 = Blk_Adaptive_RWM_MCMC(;N_its = 100000, infer_block = [true, false], data_aug_infer = [false, false, true, true, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish, moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 5],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Test Bug2/res_infpars_arSEarEI.csv", r_bug2, header = true)
CSV.write("Inference/Test Bug2/other_res_infpars_arSEarEI.csv", or_bug2, header = true)
CSV.write("Inference/Test Bug2/aug_res_infpars_arSEarEI.csv", ar_bug2, header = true)
CSV.write("Inference/Test Bug2/tuning_res_infpars_arSEarEI.csv", tr_bug2, header = true)
CSV.write("Inference/Test Bug2/update_tracker_infpars_arSEarEI.csv", ut_bug2, header = true)
