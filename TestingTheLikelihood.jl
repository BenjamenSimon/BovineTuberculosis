
include("Likelihood.jl")

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


scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements,
                                                          epi_params_true, dict_of_movements, f_to_p_dict)

calc_llh_h(scope, llh_array_cur)

calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)



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

array1 = Array{Int64, 3}(combi_array[1])
array2 = Array{Int64, 3}(combi_array[2])
array3 = Array{Float64, 3}(combi_array[3])
array4 = Array{Float64, 3}(combi_array[4])

combi_array = Union{Array{Int},Array{Float64}}[array1, array2, array3, array4]

f_to_p_dict

record_of_movements = Array(record_of_movements)

dict_of_movements



scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true, dict_of_movements)

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

@benchmark begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true, dict_of_movements)
end

# Range (min … max):  2.909 s …   3.609 s  ┊ GC (min … max): 4.13% … 19.77%
# Time  (median):     3.038 s              ┊ GC (median):    5.74%
# Time  (mean ± σ):   3.046 s ± 82.212 ms  ┊ GC (mean ± σ):  5.91% ±  1.48%
# Memory estimate: 491.97 MiB, allocs estimate: 12525523.

# Range (min … max):  2.587 s …    3.215 s  ┊ GC (min … max): 4.22% … 20.02%
# Time  (median):     2.661 s               ┊ GC (median):    6.01%
# Time  (mean ± σ):   2.689 s ± 104.071 ms  ┊ GC (mean ± σ):  6.16% ±  1.49%
# Memory estimate: 486.84 MiB, allocs estimate: 12368563.

# Range (min … max):  1.940 s …    3.897 s  ┊ GC (min … max): 7.34% … 8.95%
# Time  (median):     2.058 s               ┊ GC (median):    7.64%
# Time  (mean ± σ):   2.116 s ± 203.493 ms  ┊ GC (mean ± σ):  8.15% ± 1.52%
# Memory estimate: 486.84 MiB, allocs estimate: 12368563.

# Range (min … max):  1.966 s …   2.382 s  ┊ GC (min … max): 6.00% … 18.87%
# Time  (median):     2.055 s              ┊ GC (median):    7.80%
# Time  (mean ± σ):   2.068 s ± 59.504 ms  ┊ GC (mean ± σ):  8.10% ±  1.45%
# Memory estimate: 486.84 MiB, allocs estimate: 12368563.

# Range (min … max):  2.630 s …   3.270 s  ┊ GC (min … max): 5.06% … 19.17%
# Time  (median):     2.743 s              ┊ GC (median):    5.99%
# Time  (mean ± σ):   2.748 s ± 59.499 ms  ┊ GC (mean ± σ):  6.04% ±  1.31%
# Memory estimate: 492.79 MiB, allocs estimate: 12552367.


using ProfileView

@profview begin
  update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements, epi_params_true, dict_of_movements)
end
