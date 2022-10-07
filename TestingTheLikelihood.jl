
include("LoadPackages.jl")
include("Likelihood.jl")

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


#######################
### Checking values ###
#######################

scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

### ~~~~~~~ 1. ~~~~~~~ ###

ep1 = [0.002, 0.004, 0.015, 0.004, 0.05, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_1 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_1 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 1. ~~~~~~~")
println("the sum of LLH null pre = ", sum(llh_array_1))
println("the sum of penv LLH null pre = ", sum(p_env_llh_array_1))

llh_array_prime_1, p_env_llh_array_prime_1 = update_llh_array_ALL(scope, llh_array_1, p_env_llh_array_1, combi_array, record_of_movements, ep1, dict_of_movements, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_1))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_1))

println("the sum of LLH prime = ", sum(llh_array_prime_1))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_1))

### ~~~~~~~ 3. ~~~~~~~ ###

include("DataUpdaters.jl")
ep3 = [0.002, 0.004, 0.015, 0.004, 0.05, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_3 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_3 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 3. ~~~~~~~")
println("the sum of LLH null pre = ", sum(llh_array_3))
println("the sum of penv LLH null pre = ", sum(p_env_llh_array_3))

combi_array_prime_3 = update_pers_EPIDEMIC(combi_array, log.(ep3), f_to_p_dict, scope)

llh_array_prime_3, p_env_llh_array_prime_3 = update_llh_array_ALL(scope, llh_array_3, p_env_llh_array_3, combi_array_prime_3, record_of_movements, ep3, dict_of_movements, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_3))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_3))

println("the sum of LLH prime = ", sum(llh_array_prime_3))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_3))

println("the sum of the p_exp and p_inf before is = ", sum(combi_array[3][:, :, [4,5,6,7]]))
println("the sum of the p_exp and p_inf after is = ", sum(combi_array_prime_3[3][:, :, [4,5,6,7]]))
println("which before breaks down into c_exp_prob = ", sum(combi_array[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array[3][:, :, [7]]))
println("which after  breaks down into c_exp_prob = ", sum(combi_array_prime_3[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array_prime_3[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array_prime_3[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array_prime_3[3][:, :, [7]]))

combi_array_prime_3 = update_pers_EPIDEMIC(combi_array_prime_3, log.(ep3), f_to_p_dict, scope)

llh_array_prime_3, p_env_llh_array_prime_3 = update_llh_array_ALL(scope, llh_array_3, p_env_llh_array_3, combi_array_prime_3, record_of_movements, ep3, dict_of_movements, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_3))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_3))

println("the sum of LLH prime = ", sum(llh_array_prime_3))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_3))

println("the sum of the p_exp and p_inf before is = ", sum(combi_array[3][:, :, [4,5,6,7]]))
println("the sum of the p_exp and p_inf after is = ", sum(combi_array_prime_3[3][:, :, [4,5,6,7]]))
println("which before breaks down into c_exp_prob = ", sum(combi_array[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array[3][:, :, [7]]))
println("which after  breaks down into c_exp_prob = ", sum(combi_array_prime_3[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array_prime_3[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array_prime_3[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array_prime_3[3][:, :, [7]]))

#######################
### Checking values ###
#######################

scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

ep1 = [0.002, 0.004, 0.015, 0.004, 0.05, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_1 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_1 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 1. ~~~~~~~")
println("the sum of LLH null pre = ", sum(llh_array_1))
println("the sum of penv LLH null pre = ", sum(p_env_llh_array_1))

llh_array_prime_1, p_env_llh_array_prime_1 = update_llh_array_ALL(scope, llh_array_1, p_env_llh_array_1, combi_array, record_of_movements, ep1, dict_of_movements, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_1))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_1))

println("the sum of LLH prime = ", sum(llh_array_prime_1))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_1))



ep2 = [0.002041307650203201, 0.003953796422198999, 0.014898734975778802, 0.00399454576769233, 0.04925481068027904, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_2 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_2 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 2. ~~~~~~~")
  println("the sum of LLH null pre = ", sum(llh_array_2))
  println("the sum of penv LLH null pre = ", sum(p_env_llh_array_2))

  llh_array_prime_2, p_env_llh_array_prime_2 = update_llh_array_ALL(scope, llh_array_2, p_env_llh_array_2, combi_array, record_of_movements, ep2, dict_of_movements, f_to_p_dict)

  println("the sum of LLH null post = ", sum(llh_array_2))
  println("the sum of penv LLH null post = ", sum(p_env_llh_array_2))

  println("the sum of LLH prime = ", sum(llh_array_prime_2))
  println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_2))



include("DataUpdaters.jl")
ep3 = [0.002041307650203201, 0.003953796422198999, 0.014898734975778802, 0.00399454576769233, 0.04925481068027904, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_3 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_3 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 3. ~~~~~~~")
  println("the sum of LLH null pre = ", sum(llh_array_3))
  println("the sum of penv LLH null pre = ", sum(p_env_llh_array_3))

  combi_array_prime_3 = update_pers_EPIDEMIC(combi_array, log.(ep3), f_to_p_dict, scope)

  llh_array_prime_3, p_env_llh_array_prime_3 = update_llh_array_ALL(scope, llh_array_3, p_env_llh_array_3, combi_array_prime_3, record_of_movements, ep3, dict_of_movements, f_to_p_dict)

  println("the sum of LLH null post = ", sum(llh_array_3))
  println("the sum of penv LLH null post = ", sum(p_env_llh_array_3))

  println("the sum of LLH prime = ", sum(llh_array_prime_3))
  println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_3))

  println("the sum of the p_exp and p_inf before is = ", sum(combi_array[3][:, :, [4,5,6,7]]))
  println("the sum of the p_exp and p_inf after is = ", sum(combi_array_prime_3[3][:, :, [4,5,6,7]]))
  println("which before breaks down into c_exp_prob = ", sum(combi_array[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array[3][:, :, [7]]))
  println("which after  breaks down into c_exp_prob = ", sum(combi_array_prime_3[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array_prime_3[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array_prime_3[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array_prime_3[3][:, :, [7]]))

  combi_array_prime_3 = update_pers_EPIDEMIC(combi_array_prime_3, log.(ep3), f_to_p_dict, scope)

  llh_array_prime_3, p_env_llh_array_prime_3 = update_llh_array_ALL(scope, llh_array_3, p_env_llh_array_3, combi_array_prime_3, record_of_movements, ep3, dict_of_movements, f_to_p_dict)

  println("the sum of LLH null post = ", sum(llh_array_3))
  println("the sum of penv LLH null post = ", sum(p_env_llh_array_3))

  println("the sum of LLH prime = ", sum(llh_array_prime_3))
  println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_3))

  println("the sum of the p_exp and p_inf before is = ", sum(combi_array[3][:, :, [4,5,6,7]]))
  println("the sum of the p_exp and p_inf after is = ", sum(combi_array_prime_3[3][:, :, [4,5,6,7]]))
  println("which before breaks down into c_exp_prob = ", sum(combi_array[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array[3][:, :, [7]]))
  println("which after  breaks down into c_exp_prob = ", sum(combi_array_prime_3[3][:, :, [4]]), " b_exp_prob = ", sum(combi_array_prime_3[3][:, :, [5]]), " c_inf_prob = ", sum(combi_array_prime_3[3][:, :, [6]]), " b_inf_prob = ", sum(combi_array_prime_3[3][:, :, [7]]))




ep4 = [0.002, 0.004, 0.015, 0.004, 0.05, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_4 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_4 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 4. ~~~~~~~")
println("the sum of LLH null pre = ", sum(llh_array_4))
println("the sum of penv LLH null pre = ", sum(p_env_llh_array_4))

combi_array_prime_4 = update_pers_EPIDEMIC(combi_array, log.(ep4), f_to_p_dict, scope)

llh_array_prime_4, p_env_llh_array_prime_4 = update_llh_array_ALL(scope, llh_array_4, p_env_llh_array_4, combi_array_prime_4, record_of_movements, ep4, dict_of_movements, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_4))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_4))

println("the sum of LLH prime = ", sum(llh_array_prime_4))
println("the sum of LLH prime [3,4,8,9] = ", sum(llh_array_prime_4[:,:,[3,4,8,9]]))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_4))



ep5 = [0.002, 0.004, 0.015, 0.004, 0.05, 0.75, 0.2, 0.0048076923, 0.0048076923]

llh_array_5 = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array_5 = zeros(size(combi_array[4], 1), 360, 2)

println("~~~~~~~ 5. ~~~~~~~")
println("the sum of LLH null pre = ", sum(llh_array_5))
println("the sum of penv LLH null pre = ", sum(p_env_llh_array_5))

combi_array_prime_5 = update_pers_EPIDEMIC(combi_array, log.(ep5), f_to_p_dict, scope)

llh_array_prime_5, p_env_llh_array_prime_5 = update_llh_array_EPIDEMIC(scope, llh_array_5, p_env_llh_array_5, combi_array_prime_5, ep5, f_to_p_dict)

println("the sum of LLH null post = ", sum(llh_array_5))
println("the sum of penv LLH null post = ", sum(p_env_llh_array_5))

println("the sum of LLH prime [3,4,8,9] = ", sum(llh_array_prime_5))
println("the sum of penv LLH prime = ", sum(p_env_llh_array_prime_5))


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




############################################################
### Profiling the infection process parameter components ###
############################################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

scope = [1, 360, 1:size(combi_array[1], 1), [3,4,8,9]]

### With simd macro

llh_array = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

@benchmark begin
  update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)
end

# Range (min … max):  568.829 ms …   1.136 s  ┊ GC (min … max): 2.80% … 2.78%
# Time  (median):     600.062 ms              ┊ GC (median):    2.78%
# Time  (mean ± σ):   625.868 ms ± 73.789 ms  ┊ GC (mean ± σ):  3.38% ± 1.12%
# Memory estimate: 143.45 MiB, allocs estimate: 3453352.

# Range (min … max):  555.873 ms … 813.649 ms  ┊ GC (min … max): 2.82% … 3.11%
# Time  (median):     575.547 ms               ┊ GC (median):    2.79%
# Time  (mean ± σ):   587.611 ms ±  37.675 ms  ┊ GC (mean ± σ):  3.42% ± 1.14%
# Memory estimate: 143.45 MiB, allocs estimate: 3453352.

# Range (min … max):  554.967 ms … 609.033 ms  ┊ GC (min … max): 2.80% … 2.83%
# Time  (median):     568.087 ms               ┊ GC (median):    2.82%
# Time  (mean ± σ):   570.999 ms ±   8.830 ms  ┊ GC (mean ± σ):  3.43% ± 1.12%
# Memory estimate: 143.45 MiB, allocs estimate: 3453352.

# Range (min … max):  555.009 ms … 816.497 ms  ┊ GC (min … max): 2.84% … 4.69%
# Time  (median):     578.012 ms               ┊ GC (median):    2.78%
# Time  (mean ± σ):   582.758 ms ±  21.902 ms  ┊ GC (mean ± σ):  3.39% ± 1.12%
# Memory estimate: 143.45 MiB, allocs estimate: 3453352.

### Without simd macro

llh_array = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

@benchmark begin
  update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)
end

# Range (min … max):  559.294 ms … 667.230 ms  ┊ GC (min … max): 2.78% … 5.50%
# Time  (median):     578.001 ms               ┊ GC (median):    2.72%
# Time  (mean ± σ):   580.698 ms ±  12.198 ms  ┊ GC (mean ± σ):  3.34% ± 1.11%
# Memory estimate: 143.45 MiB, allocs estimate: 3453460.

# Range (min … max):  609.629 ms … 997.158 ms  ┊ GC (min … max): 2.68% … 4.80%
# Time  (median):     780.922 ms               ┊ GC (median):    2.71%
# Time  (mean ± σ):   765.492 ms ±  68.692 ms  ┊ GC (mean ± σ):  3.28% ± 1.12%
# Memory estimate: 143.45 MiB, allocs estimate: 3453460.

# Range (min … max):  639.210 ms …   1.018 s  ┊ GC (min … max): 2.61% … 3.51%
# Time  (median):     799.836 ms              ┊ GC (median):    2.73%
# Time  (mean ± σ):   804.944 ms ± 35.219 ms  ┊ GC (mean ± σ):  3.27% ± 1.11%
# Memory estimate: 143.45 MiB, allocs estimate: 3453460.

# Range (min … max):  730.403 ms …    1.431 s  ┊ GC (min … max): 2.64% … 2.92%
# Time  (median):     854.116 ms               ┊ GC (median):    2.78%
# Time  (mean ± σ):   903.513 ms ± 122.074 ms  ┊ GC (mean ± σ):  3.17% ± 1.09%
# Memory estimate: 143.45 MiB, allocs estimate: 3453460.

# Range (min … max):  665.320 ms … 763.391 ms  ┊ GC (min … max): 2.30% … 4.60%
# Time  (median):     701.870 ms               ┊ GC (median):    2.26%
# Time  (mean ± σ):   699.618 ms ±  16.495 ms  ┊ GC (mean ± σ):  2.81% ± 0.94%
# Memory estimate: 143.45 MiB, allocs estimate: 3453460.

using ProfileView

array1 = Array{Int64, 3}(combi_array[1])
array2 = Array{Int64, 3}(combi_array[2])
array3 = Array{Float64, 3}(combi_array[3])
array4 = Array{Float64, 3}(combi_array[4])

combi_array = Union{Array{Int},Array{Float64}}[array1, array2, array3, array4]

scope = [1, 360, 1:size(combi_array[1], 1), [3,4,8,9]]

llh_array = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

@profview begin
  update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)
end

@profview begin
  for i in 1:10000
    exposures_llh_i_t(;position = 1, t = 100, combi_array)
  end
end

@profview begin
  for i in 1:10000
    infections_llh_i_t(;position = 1, t = 100, combi_array)
  end
end

@profview begin
  for i in 1:10000
    p_env_llh_k_t(;p_position = 1, t = 100, combi_array, epi_params = epi_params_true)
  end
end

@profview begin
  position = 1
  t = 100
  for i in 1:10000
    combi_array[1][position, t, [10,11,12]]
    combi_array[2][position, t, [13,14]]
    combi_array[3][position, t, 4]

    exposures(States_postM = combi_array[1][position, t, [10,11,12]],
                             new_EandI = combi_array[2][position, t, [13,14]],
                              exp_prob = combi_array[3][position, t, 4])
  end
end

@profview begin
  position = 1
  t = 100
  for i in 1:10000
    exposures(States_postM = combi_array[1][position, t, [10,11,12]],
                             new_EandI = combi_array[2][position, t, [13,14]],
                              exp_prob = combi_array[3][position, t, 4])
  end
end

@profview begin
  position = 1
  t = 100
  for i in 1:10000
    exposures(States_postM = combi_array[1][position, t, [10,11,12]],
                             new_EandI = combi_array[2][position, t, [13,14]],
                              exp_prob = combi_array[3][position, t, 4])
  end
end




###########################
### Profiling Functions ###
###########################

array1 = Array{Int64, 3}(combi_array[1])
array2 = Array{Int64, 3}(combi_array[2])
array3 = Array{Float64, 3}(combi_array[3])
array4 = Array{Float64, 3}(combi_array[4])

combi_array = Union{Array{Int},Array{Float64}}[array1, array2, array3, array4]

scope = [1, 360, 1:size(combi_array[1], 1), [3,4,8,9]]

struct Scope
  t_start::Int64
  t_end::Int64
  range_::UnitRange{Int64}
  llh_inx::Vector{Int64}
end

scope = Scope(1, 360, 1:size(combi_array[1], 1), [3,4,8,9])


llh_array = zeros(size(combi_array[1], 1), 360, 13)
p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)


using ProfileView

@profview update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)

warntype_last()

ProfileView.view(nothing)


using Profile
using PProf

# collect a profile
@profile update_llh_array_EPIDEMIC(scope, llh_array, p_env_llh_array, combi_array, epi_params_true, f_to_p_dict)

# Export pprof profile and open interactive profiling web interface.
pprof()
