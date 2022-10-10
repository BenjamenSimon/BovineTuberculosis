
#######################
### Data Structures ###
#######################

struct Scope
  t_start::Int64
  t_end::Int64
  h_positions::UnitRange{Int64}
  h_llh_indices::Vector{Int64}
end

# scope = Scope(1, 360, 1:size(DATA_res_and_track[1], 1), [3,4,8,9])
# scope.t_start

################################################
### Multivariate Hypergeometric pdf function ###
################################################

function log_pdf_mvhyper(K, k)
  # K is the number of objects of each type i
  # k is a vector of the number of succeses of each type i
  # N is the total number of objects
  # n is the number of draws total

  N = sum(K)
  n = sum(k)

  m = length(K) # number of groups

  x = fill(0.0, m) # results

  for i in 1:m
    Kck_i = binomial(BigInt(K[i]), BigInt(k[i]))

    x[i] = log(Kck_i)
  end

  top = sum(x)
  bot = binomial(BigInt(N), BigInt(n))

  prob = top-log(bot)

  return(prob)
end


###################################################
### Functions to calculate movement likelihoods ###
###################################################

function moves_MHG(;States_Moves::Vector{Int64} = [90000, 5000, 5000], Moves::Vector{Int64})

  prob::Float64 = log_pdf_mvhyper(States_Moves, Moves)

  return(prob)
end


#################################################################
### Functions to calculate exposure and infection likelihoods ###
#################################################################

function exposures(;States_postM::Vector{Int64}, new_EandI::Vector{Int64}, exp_prob::Float64)

  d = Binomial(States_postM[1], exp_prob)

  prob = logpdf(d, new_EandI[1])

  return(prob::Float64)
end

function infections(;States_postM, new_EandI, inf_prob)

  d = Binomial(States_postM[2], inf_prob)

  prob = logpdf(d, new_EandI[2])

  return(prob)
end


####################################################
### Functions to calculate Detection likelihoods ###
####################################################

function exposed_detections(;States_postEI, dets, exp_det_prob)

  d = Binomial(States_postEI[2], exp_det_prob)

  prob = logpdf(d, dets[1])

  return(prob)
end

function infectious_detections(;States_postEI, dets, inf_det_prob)

  d = Binomial(States_postEI[3], inf_det_prob)

  prob = logpdf(d, dets[2])

  return(prob)
end


########################################################
### Functions to calculate Birth & Death likelihoods ###
########################################################

function cattle_death(;States_postDet, Deaths)

  prob = log_pdf_mvhyper(States_postDet, Deaths)

  return(prob)
end

function badger_births(;States_postEI, b_births, θ_bb)

  bN = sum(States_postEI)

  rate = θ_bb * bN

  d = Poisson(rate)

  prob = logpdf(d, b_births)

  return(prob)
end


function badger_S_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(States_postEI[1], θ_bd)

  prob = logpdf(d, Deaths[1])

  return(prob)
end

function badger_E_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(States_postEI[2], θ_bd)

  prob = logpdf(d, Deaths[2])

  return(prob)
end

function badger_I_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(States_postEI[3], θ_bd)

  prob = logpdf(d, Deaths[3])

  return(prob)
end


##################################################################
### Functions to calculate Environmental Reservoir likelihoods ###
##################################################################

function new_pressure(;pI_init, new_pres)

  d = Poisson( sum(pI_init) )

  prob = logpdf(d, new_pres)

  return(prob)
end


function remaining_pressure(;scaled_p_env_prev, ϵ, remaining_pres)

  d = Binomial( round(Int, scaled_p_env_prev), (1-ϵ) )

  prob = logpdf(d, remaining_pres)

  return(prob)
end


#########################################
### Constructing the likelihood array ###
#########################################

# llh_array = Array{Array{BigFloat, 2}, 1}(undef,100)
# for k in 1:100
#   llh_array[k] = Array{BigFloat, 2}(undef, 360, 14)
# end

# 1. m_off_llh, 2. indv_moves_off_llh, 3. c_exp_llh, 4. c_inf_llh, 5. exp_det_llh,
# 6. inf_det_llh, 7. c_dth_llh, 8. b_exp_llh, 9. b_inf_llh, 10. b_bths_llh,
# 11. bS_dths_llh, 12. bE_dths_llh, 13. bI_dths_llh, 14. m_on_out_llh


# p_env_llh_array = Array{Array{BigFloat, 2}, 1}(undef,77)
# for k in 1:77
#   p_env_llh_array[k] = Array{BigFloat, 2}(undef, 360, 2)
# end



# cStates_Moves_i_t = DATA_res_and_track[1][position, t, [7,8,9]]   # :cS_Moves, :cE_Moves, :cI_Moves
# cStates_postM_i_t = DATA_res_and_track[1][position, t, [10,11,12]]  # :cS_postM, :cE_postM, :cI_postM
# cStates_postEI_i_t = DATA_res_and_track[1][position, t, [13,14,15]]  # :cS_postEI, :cE_postEI, :cI_postEI
# cStates_postDet_i_t = DATA_res_and_track[1][position, t, [16,17,18]]  #  :cS_postDet, :cE_postDet, :cI_postDet
#
# bStates_init_i_t = DATA_res_and_track[1][position, t, [22,23,24]]  #  :bS_init, :bE_init, :bI_init
# bStates_postEI_i_t = DATA_res_and_track[1][position, t, [25,26,27]]  #  :bS_postEI, :bE_postEI, :bI_postEI
#
# Moves_off_i_t = DATA_res_and_track[2][position, t, [7,8,9]]  #  :sus_off, :exp_off, :inf_off
# Moves_on_out_i_t = DATA_res_and_track[2][position, t, [10,11,12]]  #  :S_on_out, :E_on_out, :I_on_out
#
# cEandI_i_t = DATA_res_and_track[2][position, t, [13,14]]  #  :c_new_exp, :c_new_inf
# cEandI_probs_i_t = DATA_pers_and_parish[1][position, t, [4,6]]  #  :c_exp_prob, :c_inf_prob
#
# bEandI_i_t = DATA_res_and_track[2][position, t, [15,16]]  #  :b_new_exp, :b_new_inf
# bEandI_probs_i_t = DATA_pers_and_parish[1][position, t, [5,7]]  #  :b_exp_prob, :b_inf_prob
#
# dets_i_t = DATA_res_and_track[2][position, t, [19,20]]  #  :E_detected, :I_detected
# test_occur_i_t = DATA_res_and_track[2][position, t, [17]]  #  :test_occur
#
# cDeaths_i_t = DATA_res_and_track[2][position, t, [22,23,24]]  #  :c_S_death, :c_E_death, :c_I_death
#
# b_birth_i_t = DATA_res_and_track[2][position, t, [25]]  #  :b_birth
# bDeaths_i_t = DATA_res_and_track[2][position, t, [26,27,28]]  #  :b_S_death, :b_E_death, :b_I_death
#
# inf_det_prob = epi_params[6]
# exp_det_prob = inf_det_prob * epi_params[7]
#
# θ_bb = epi_params[8]
# θ_bd = epi_params[9]
#
#
#
# pcStates_init_k_t = DATA_pers_and_parish[2][position, t, [4,5,6]]
# pbStates_init_k_t = DATA_pers_and_parish[2][position, t, [10,11,12]]
#
# pI_init_k_t = sum(DATA_pers_and_parish[2][position, t, [6,12]])
#
# remaining_pressure_k_t = DATA_pers_and_parish[2][position, t, [16]]
# new_pressure_k_t = DATA_pers_and_parish[2][position, t, [17]]
# scaling_k = DATA_pers_and_parish[2][position, t, [18]]
#
# p_env_prev_k_t = DATA_pers_and_parish[2][position, t, [19]]
#
# ϵ = epi_params[5]



########################################
### Caclulate LLH Elements Functions ###
########################################

function movements_llh_i_t(;position, t, DATA_res_and_track, movement_records, movement_dict)

  # MOVEMENT PROCESS

  movement_index = movement_dict[(position, t)]

  movement_records_i_t = movement_records[movement_index, :]

  # Probability of moving animal states off farm of all animals on farm
  m_off_llh::Float64 = moves_MHG(States_Moves = DATA_res_and_track[1][position, t, [7,8,9]],
                               Moves = DATA_res_and_track[2][position, t, [7,8,9]])

  # Probability of moving animal states to each farm of animals moved off
  indv_moves_off_llh::Float64 = 0.0

  @inbounds for j in 1:size(movement_records_i_t, 1)

    indv_moves_off_llh = indv_moves_off_llh + Float64(moves_MHG(States_Moves = movement_records_i_t[j, 4:6],
                                                                        Moves = movement_records_i_t[j, 7:9]))

  end

  # Probability of animal states moving onto the farm from outside farms of interest
  # m_on_out_llh = moves_MHG(States_Moves = [10800, 600, 600],
  #                                 Moves = DATA_res_and_track[2][position, t, [10,11,12]])

  # [1,2]
  # not returning m_on_out_llh at the moment

  llh_res = [m_off_llh, indv_moves_off_llh]

  return(llh_res)
end # ALL MOVEMENTS

function movements_total_llh_i_t(;position, t, DATA_res_and_track, movement_records, movement_dict)

  # MOVEMENT PROCESS

  # Probability of moving animal states off farm of all animals on farm
  m_off_llh = moves_MHG(States_Moves = DATA_res_and_track[1][position, t, [7,8,9]],
                               Moves = DATA_res_and_track[2][position, t, [7,8,9]])

  # [1]
  return(m_off_llh)
end # JUST TOTAL MOVEMENTS OFF


function c_epidemic_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # CATTLE EPIDEMIC PROCESS

  c_exp_llh = exposures(States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                           new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                            exp_prob = DATA_pers_and_parish[1][position, t, 4])

  c_inf_llh = infections(States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                            new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                             inf_prob = DATA_pers_and_parish[1][position, t, 6])

  # [3,4]
  return([c_exp_llh, c_inf_llh])
end


function b_epidemic_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # BADGER EPIDEMIC PROCESS

  b_exp_llh = exposures(;States_postM = DATA_res_and_track[1][position, t, [22,23,24]],
                            new_EandI = DATA_res_and_track[2][position, t, [15,16]],
                             exp_prob = DATA_pers_and_parish[1][position, t, 5])

  b_inf_llh = infections(;States_postM = DATA_res_and_track[1][position, t, [22,23,24]],
                             new_EandI = DATA_res_and_track[2][position, t, [15,16]],
                              inf_prob = DATA_pers_and_parish[1][position, t, 7])

  # [8,9]
  return([b_exp_llh, b_inf_llh])
end


function exposures_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # EXPOSURE PROCESS

  c_exp_llh::Float64 = exposures(States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                           new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                            exp_prob = DATA_pers_and_parish[1][position, t, 4])

  b_exp_llh::Float64 = exposures(;States_postM = DATA_res_and_track[1][position, t, [22,23,24]],
                            new_EandI = DATA_res_and_track[2][position, t, [15,16]],
                             exp_prob = DATA_pers_and_parish[1][position, t, 5])

  # [3,8]
  return([c_exp_llh::Float64, b_exp_llh::Float64])
end


function infections_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # INFECTION PROCESS

  c_inf_llh = infections(;States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                             new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                              inf_prob = DATA_pers_and_parish[1][position, t, 6])

  b_inf_llh = infections(;States_postM = DATA_res_and_track[1][position, t, [22,23,24]],
                             new_EandI = DATA_res_and_track[2][position, t, [15,16]],
                              inf_prob = DATA_pers_and_parish[1][position, t, 7])

  # [4,9]
  return([c_inf_llh, b_inf_llh])
end


function detection_llh_i_t(;position, t, DATA_res_and_track, epi_params)

  test_occur_i_t = DATA_res_and_track[2][position, t, [17]][1]  #  :test_occur

  # DETECTION PROCESS
    # where test_occur_i_t = 1 if its a test week on farm i, 0 otherwise

  exp_det_llh = test_occur_i_t * exposed_detections(States_postEI = DATA_res_and_track[1][position, t, [13,14,15]],
                                                           dets = DATA_res_and_track[2][position, t, [19,20]],
                                                   exp_det_prob = (epi_params[6] * epi_params[7]))

  inf_det_llh = test_occur_i_t * infectious_detections(States_postEI = DATA_res_and_track[1][position, t, [13,14,15]],
                                                              dets = DATA_res_and_track[2][position, t, [19,20]],
                                                      inf_det_prob = epi_params[6])

  # [5,6]
  return([exp_det_llh, inf_det_llh])
end


function c_birth_death_llh_i_t(;position, t, DATA_res_and_track)

  # CATTLE BIRTH AND DEATH PROCESS

  c_dth_llh = cattle_death(States_postDet = DATA_res_and_track[1][position, t, [16,17,18]],
                                   Deaths = DATA_res_and_track[2][position, t, [22,23,24]])

  # [7]
  return(c_dth_llh)
end


function b_birth_death_llh_i_t(;position, t, DATA_res_and_track, epi_params)

  # BADGER BIRTH AND DEATH PROCESS

  b_bths_llh = badger_births(;States_postEI = DATA_res_and_track[1][position, t, [25,26,27]],
                                   b_births = DATA_res_and_track[2][position, t, 25],
                                       θ_bb = epi_params[8])

  bS_dths_llh = badger_S_deaths(;States_postEI = DATA_res_and_track[1][position, t, [25,26,27]],
                                        Deaths = DATA_res_and_track[2][position, t, [26,27,28]],
                                          θ_bd = epi_params[9])

  bE_dths_llh = badger_E_deaths(;States_postEI = DATA_res_and_track[1][position, t, [25,26,27]],
                                        Deaths = DATA_res_and_track[2][position, t, [26,27,28]],
                                          θ_bd = epi_params[9])

  bI_dths_llh = badger_I_deaths(;States_postEI = DATA_res_and_track[1][position, t, [25,26,27]],
                                        Deaths = DATA_res_and_track[2][position, t, [26,27,28]],
                                          θ_bd = epi_params[9])

  # [10,11,12,13]
  return([b_bths_llh, bS_dths_llh, bE_dths_llh, bI_dths_llh])
end

function p_env_llh_k_t(;p_position, t, DATA_pers_and_parish, epi_params)

  # PARISH ENVIRONMENTAL RESERVOIR

  new_pressure_llh = new_pressure(;pI_init = DATA_pers_and_parish[2][p_position, t, [6,12]],
                                   new_pres = DATA_pers_and_parish[2][p_position, t, 17])


  remaining_pressure_llh = remaining_pressure(;scaled_p_env_prev = prod(DATA_pers_and_parish[2][p_position, t, [18, 19]]),
                                               ϵ = epi_params[5],
                                               remaining_pres = DATA_pers_and_parish[2][p_position, t, 16])

  return([new_pressure_llh, remaining_pressure_llh])
end



##################################
### Update LLH Array Functions ###
##################################

@views function update_llh_array_ALL(scope::Scope, llh_array_cur, p_env_llh_array_cur, DATA_res_and_track, DATA_pers_and_parish, movement_records, epi_params, movement_dict, f_to_p_structs::Vector{Farm_Parish_info})

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions

  llh_array_new = deepcopy(llh_array_cur)
  p_env_llh_array_new = deepcopy(p_env_llh_array_cur)

  parishes_to_update = Int64[]

  @inbounds for pos in positions
    for t in t_start:t_end

      llh_array_new[pos, t, [1,2]] = movements_llh_i_t(;position = pos, t = t, DATA_res_and_track, movement_records = movement_records, movement_dict = movement_dict)

      llh_array_new[pos, t, [3,8]] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [4,9]] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [5,6]] = detection_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

      llh_array_new[pos, t, 7] = c_birth_death_llh_i_t(;position = pos, t = t, DATA_res_and_track)

      llh_array_new[pos, t, [10,11,12,13]] = b_birth_death_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

    end

    push!(parishes_to_update, f_to_p_structs[pos].parish_position)
  end

  @inbounds for p_pos in parishes_to_update
    for t in t_start:t_end

      p_env_llh_array_new[p_pos, t, [1,2]] = p_env_llh_k_t(;p_position = p_pos, t = t, DATA_pers_and_parish, epi_params = epi_params)

    end
  end

  return(llh_array_new, p_env_llh_array_new)
end

@views function update_llh_array_ALL_excindvmoves(scope::Scope, llh_array_cur, p_env_llh_array_cur, DATA_res_and_track, DATA_pers_and_parish, movement_records, epi_params, movement_dict, f_to_p_structs::Vector{Farm_Parish_info})

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions

  llh_array_new = deepcopy(llh_array_cur)
  p_env_llh_array_new = deepcopy(p_env_llh_array_cur)

  parishes_to_update = Int64[]

  @inbounds for pos in positions
    for t in t_start:t_end

      llh_array_new[pos, t, 1] = movements_total_llh_i_t(;position = pos, t = t, DATA_res_and_track, movement_records = movement_records, movement_dict = movement_dict)

      llh_array_new[pos, t, [3,8]] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [4,9]] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [5,6]] = detection_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

      llh_array_new[pos, t, 7] = c_birth_death_llh_i_t(;position = pos, t = t, DATA_res_and_track)

      llh_array_new[pos, t, [10,11,12,13]] = b_birth_death_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

    end

    push!(parishes_to_update, f_to_p_structs[pos].parish_position)
  end

  @inbounds for p_pos in parishes_to_update
    for t in t_start:t_end

      p_env_llh_array_new[p_pos, t, [1,2]] = p_env_llh_k_t(;p_position = p_pos, t = t, DATA_pers_and_parish, epi_params = epi_params)

    end
  end

  return(llh_array_new, p_env_llh_array_new)
end

@views function update_llh_array_EPIDEMIC(scope::Scope, llh_array_cur, p_env_llh_array_cur, DATA_res_and_track, DATA_pers_and_parish, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions

  llh_array_new = deepcopy(llh_array_cur)
  p_env_llh_array_new = deepcopy(p_env_llh_array_cur)

  parishes_to_update = Int64[]

  @inbounds for pos in positions
    for t in t_start:t_end

      llh_array_new[pos, t, [3,8]] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [4,9]] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

    end

    push!(parishes_to_update, f_to_p_structs[pos].parish_position)
  end

  @inbounds for p_pos in parishes_to_update
    for t in t_start:t_end

      p_env_llh_array_new[p_pos, t, [1,2]] = p_env_llh_k_t(;p_position = p_pos, t = t, DATA_pers_and_parish, epi_params = epi_params)

    end
  end

  return(llh_array_new, p_env_llh_array_new)
end

@views function update_llh_array_DETECTION(scope::Scope, llh_array_cur, p_env_llh_array_cur, DATA_res_and_track, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions

  llh_array_new = deepcopy(llh_array_cur)

  parishes_to_update = Int64[]

  @inbounds for pos in positions
    for t in t_start:t_end

      llh_array_new[pos, t, [5,6]] = detection_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

    end
  end

  return(llh_array_new, p_env_llh_array_cur)
end

@views function update_llh_array_BBD(scope::Scope, llh_array_cur, DATA_res_and_track, epi_params)

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions

  llh_array_new = deepcopy(llh_array_cur)
  p_env_llh_array_new = deepcopy(p_env_llh_array_cur)

  parishes_to_update = Int64[]

  @inbounds for pos in positions
    for t in t_start:t_end

      llh_array_new[pos, t, [10,11,12,13]] = b_birth_death_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

    end
  end

  return(llh_array_new)
end

###############################
### Sum LLH Array Functions ###
###############################

function calc_llh_h(scope::Scope, llh_array_cur)

  t_start = scope.t_start
  t_end = scope.t_end
  positions = scope.h_positions
  h_llh_indices = scope.h_llh_indices

  llh = 0.0

  for pos in positions
    llh += sum(llh_array_cur[pos, t_start:t_end, h_llh_indices])
  end

  return(llh)
end

function calc_llh_h_and_p(scope::Scope, llh_array_cur, p_env_llh_array_cur)

  llh = calc_llh_h(scope, llh_array_cur)

  t_start = scope.t_start
  t_end = scope.t_end

  p_env_llh = 0.0

  for p_pos in 1:size(p_env_llh_array_cur, 1) #p_positions
    p_env_llh += sum(p_env_llh_array_cur[p_pos, t_start:t_end, :])
  end

  return(llh + p_env_llh)
end
