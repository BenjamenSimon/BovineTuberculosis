
include("Likelihood.jl")

#########################################
### Constructing the likelihood array ###
#########################################

# llh_array = Array{Array{BigFloat, 2}, 1}(undef,100)
# for k in 1:100
#   llh_array[k] = Array{BigFloat, 2}(undef, 360, 14)
# end

llh_array = zeros(100, 360, 13)

# 1. m_off_llh, 2. m_on_out_llh, 3. c_exp_llh, 4. c_inf_llh, 5. exp_det_llh,
# 6. inf_det_llh, 7. c_dth_llh, 8. b_exp_llh, 9. b_inf_llh, 10. b_bths_llh,
# 11. bS_dths_llh, 12. bE_dths_llh, 13. bI_dths_llh, 14. indv_moves_off_llh

# p_env_llh_array = Array{Array{BigFloat, 2}, 1}(undef,77)
# for k in 1:77
#   p_env_llh_array[k] = Array{BigFloat, 2}(undef, 360, 2)
# end

p_env_llh_array = zeros(77, 360, 2)

##############################
### Testing the likelihood ###
##############################

# simple_movements_llh_i_t(;pos_id = 1, t = 1, combi_array)
#
# individual_movements_off_llh_i_t(;pos_id = 1, t = 1, movement_records = record_of_movements)

movements_llh_i_t(;position = 1, t = 1, combi_array, movement_records = record_of_movements, movement_dict = dict_of_movements)


c_epidemic_llh_i_t(;pos_id = 1, t = 353, combi_array)

b_epidemic_llh_i_t(;pos_id = 1, t = 353, combi_array)


exposures_llh_i_t(;pos_id = 1, t = 353, combi_array)

infections_llh_i_t(;pos_id = 1, t = 353, combi_array)


detection_llh_i_t(;pos_id = 1, t = 353, combi_array, epi_params = epi_params_true)


c_birth_death_llh_i_t(;pos_id = 1, t = 353, combi_array)

b_birth_death_llh_i_t(;pos_id = 1, t = 353, combi_array, epi_params = epi_params_true)


p_env_llh_k_t(;p_pos_id = 1, t = 353, combi_array, epi_params = epi_params_true)


scope = [1, 360, 1:100, 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array,
                                                            record_of_movements, epi_params_true)

calc_llh_h(scope, llh_array_cur)

calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)



###################################
### Benchmarking the likelihood ###
###################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

@benchmark begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true)
end

# Range (min … max):  47.284 s … 65.815 s  ┊ GC (min … max): 4.83% … 4.82%
# Time  (median):     51.466 s             ┊ GC (median):    4.72%
# Time  (mean ± σ):   52.266 s ±  5.060 s  ┊ GC (mean ± σ):  4.88% ± 0.44%
# Memory estimate: 122.74 GiB, allocs estimate: 41577218.


using ProfileView

@profview begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true)
end



record_of_movements[findall((record_of_movements[:,1] .== 1) .& (record_of_movements[:,2] .== 1)), :]
