
include("LoadPackages.jl")
include("Likelihood.jl")
include("DataUpdaters.jl")
include("Proposals.jl")
include("Posterior.jl")
include("MCMCfuncs_badgers.jl")

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
d_θ_bb = Gamma(2, 0.25/104)
d_θ_bd = Gamma(2, 0.25/104)


epi_params_dists = [d_β_c, d_β_b, d_γ, d_F, d_ϵ, d_ρ, d_ρ_E, d_θ_bb, d_θ_bd]


##############################
#### 1. Badger parameters ####
##############################

r1, or1, ar1, tr1, ut1 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [false, false, true], data_aug_infer = [false, false, false, false, false, false, false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Badger/Results_1_[Bad]/res_1_[Bad].csv", r1, header = true)
CSV.write("Inference/Badger/Results_1_[Bad]/other_res_1_[Bad].csv", or1, header = true)
CSV.write("Inference/Badger/Results_1_[Bad]/aug_res_1_[Bad].csv", ar1, header = true)
CSV.write("Inference/Badger/Results_1_[Bad]/tuning_res_1_[Bad].csv", tr1, header = true)
CSV.write("Inference/Badger/Results_1_[Bad]/update_tracker_1_[Bad].csv", ut1, header = true)


###########################
#### 2. All parameters ####
###########################

r2, or2, ar2, tr2, ut2 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, true, true], data_aug_infer = [false, false, false, false, false, false, false, false, false, false, false, false, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Badger/Results_2_[Inf][Det][Bad]/res_2_[Inf][Det][Bad].csv", r2, header = true)
CSV.write("Inference/Badger/Results_2_[Inf][Det][Bad]/other_res_2_[Inf][Det][Bad].csv", or2, header = true)
CSV.write("Inference/Badger/Results_2_[Inf][Det][Bad]/aug_res_2_[Inf][Det][Bad].csv", ar2, header = true)
CSV.write("Inference/Badger/Results_2_[Inf][Det][Bad]/tuning_res_2_[Inf][Det][Bad].csv", tr2, header = true)
CSV.write("Inference/Badger/Results_2_[Inf][Det][Bad]/update_tracker_2_[Inf][Det][Bad].csv", ut2, header = true)


#########################################
#### 3. All parameters and SE and EI ####
#########################################

r3, or3, ar3, tr3, ut3 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, true, true], data_aug_infer = [false, false, false, false, false, false, false, false, true, true, true, true, false, false],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", r3, header = true)
CSV.write("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/other_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", or3, header = true)
CSV.write("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/aug_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", ar3, header = true)
CSV.write("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/tuning_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", tr3, header = true)
CSV.write("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/update_tracker_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", ut3, header = true)


######################################
#### 4. All parameters All Badger ####
######################################

r4, or4, ar4, tr4, ut4 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, true, true], data_aug_infer = [false, false, false, false, false, false, false, false, true, true, true, true, true, true],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", r4, header = true)
CSV.write("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/other_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", or4, header = true)
CSV.write("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/aug_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", ar4, header = true)
CSV.write("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/tuning_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", tr4, header = true)
CSV.write("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/update_tracker_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", ut4, header = true)


#################
#### 5. All  ####
#################

r5, or5, ar5, tr5, ut5 = Blk_Adaptive_RWM_MCMC(;N_its = 1000, infer_block = [true, true, true], data_aug_infer = [true, true, true, true, true, true, true, true, true, true, true, true, true, true],
                          DATA_res_and_track = DATA_res_and_track, DATA_pers_and_parish = DATA_pers_and_parish,
                          moves_record = record_of_movements,
                          params_init = epi_params_true, tuning = [0.05, 45, 0.01, 15, 0.01, 15],
                          dict_of_movements = dict_of_movements, f_to_p_structs = f_to_p_structs,
                          ids_to_pos_dict = ids_to_pos_dict)

CSV.write("Inference/Badger/Results_5_[All]/res_5_[All].csv", r5, header = true)
CSV.write("Inference/Badger/Results_5_[All]/other_res_5_[All].csv", or5, header = true)
CSV.write("Inference/Badger/Results_5_[All]/aug_res_5_[All].csv", ar5, header = true)
CSV.write("Inference/Badger/Results_5_[All]/tuning_res_5_[All].csv", tr5, header = true)
CSV.write("Inference/Badger/Results_5_[All]/update_tracker_5_[All].csv", ut5, header = true)
