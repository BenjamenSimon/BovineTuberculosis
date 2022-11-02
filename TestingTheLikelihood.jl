
include("LoadPackages.jl")
include("Likelihood.jl")

#################
### LOAD DATA ###
#################

# combi_array = load("Data/Set 1/combi_array.jld2")["array"]

# combi_array = load("Data/Set 1/combi_array_unnamed.jld2")["array"]

DATA_res_and_track = load("Data/Set 1/DATA_res_and_track.jld2")["array"]
DATA_pers_and_parish = load("Data/Set 1/DATA_pers_and_parish.jld2")["array"]

# f_to_p_dict = load("Data/Set 1/f_to_p_dict.jld2")["dict"]

f_to_p_structs = load("Data/Set 1/f_to_p_structs.jld2")["struct"]

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

#########################################
### Constructing the likelihood array ###
#########################################

llh_array = zeros(size(combi_array[1], 1), 360, 13)

# 1. m_off_llh, 2. m_on_out_llh, 3. c_exp_llh, 4. c_inf_llh, 5. exp_det_llh,
# 6. inf_det_llh, 7. c_dth_llh, 8. b_exp_llh, 9. b_inf_llh, 10. b_bths_llh,
# 11. bS_dths_llh, 12. bE_dths_llh, 13. bI_dths_llh

p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)


##############################
### Testing the likelihood ###
##############################

# simple_movements_llh_i_t(;position = 1, t = 1, combi_array)
#
# individual_movements_off_llh_i_t(;position = 1, t = 1, movement_records = record_of_movements)

movements_llh_i_t(;position = 1, t = 1, combi_array, movement_records = record_of_movements, movement_dict = dict_of_movements)


c_epidemic_llh_i_t(;position = 1, t = 353, combi_array)

b_epidemic_llh_i_t(;position = 1, t = 353, combi_array)


exposures_llh_i_t(;position = 1, t = 353, combi_array)

infections_llh_i_t(;position = 1, t = 353, combi_array)


detection_llh_i_t(;position = 1, t = 353, combi_array, epi_params = epi_params_true)


c_birth_death_llh_i_t(;position = 1, t = 353, combi_array)

b_birth_death_llh_i_t(;position = 1, t = 353, combi_array, epi_params = epi_params_true)


p_env_llh_k_t(;p_position = 1, t = 353, combi_array, epi_params = epi_params_true)


scope = [1, 360, Vector(1:size(combi_array[1]), 1), 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements,
                                                          epi_params_true, dict_of_movements, f_to_p_dict)

calc_llh_h(scope, llh_array_cur)

calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)


#######################
### Checking values ###
#######################



###################################
### Benchmarking the likelihood ###
###################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

@benchmark begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true, dict_of_movements)
end

# Range (min … max):  47.284 s … 65.815 s  ┊ GC (min … max): 4.83% … 4.82%
# Time  (median):     51.466 s             ┊ GC (median):    4.72%
# Time  (mean ± σ):   52.266 s ±  5.060 s  ┊ GC (mean ± σ):  4.88% ± 0.44%
# Memory estimate: 122.74 GiB, allocs estimate: 41577218.

# Range (min … max):  6.607 s …    7.388 s  ┊ GC (min … max): 6.32% … 11.96%
# Time  (median):     6.879 s               ┊ GC (median):    6.26%
# Time  (mean ± σ):   6.882 s ± 119.459 ms  ┊ GC (mean ± σ):  6.33% ±  0.67%
# Memory estimate: 2.30 GiB, allocs estimate: 42884363.

# Range (min … max):   8.522 s … 17.698 s  ┊ GC (min … max): 6.69% … 5.38%
# Time  (median):     10.762 s             ┊ GC (median):    6.51%
# Time  (mean ± σ):   11.711 s ±  2.366 s  ┊ GC (mean ± σ):  6.21% ± 0.87%
# Memory estimate: 2.41 GiB, allocs estimate: 44729867.

using ProfileView

@profview begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true, dict_of_movements)
end


#########################################################
### Benchmarking the likelihood with not named arrays ###
#########################################################

scope = Scope(1, 360, Vector(1:size(combi_array[1], 1)), [3,4,8,9])


llh_array = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_structs)

update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, DATA_res_and_track, DATA_pers_and_parish, epi_params_true, f_to_p_structs)


using ProfileView

@profview update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, DATA_res_and_track, DATA_pers_and_parish, epi_params_true, f_to_p_structs)

warntype_last()

ProfileView.view(nothing)


using Profile
using PProf

# collect a profile
@profile update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)

# Export pprof profile and open interactive profiling web interface.
pprof()
