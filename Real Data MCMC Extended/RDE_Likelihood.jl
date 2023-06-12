
#######################
### Data Structures ###
#######################

struct Scope
  t_start::Int64
  t_end::Int64
  h_positions::Vector{Int64}
  h_llh_indices::Vector{Int64}
end

# scope = Scope(1, 360, 1:size(DATA_res_and_track[1], 1), [3,4,8,9])
# scope.t_start

###################
### The indicies###
###################

###~~~ Likelihood ~~~###

# 1. m_off_llh
# 2. indv_moves_off_llh
# 3. c_exp_llh
# 4. c_inf_llh
# 5. exp_det_llh
# 6. inf_det_llh
# 7. c_dth_llh
# 8. m_on_out_llh

###~~~ Results ~~~###

# 1. t
# 2. farm_uid
# 3. pos
# 4. cS_init      5. cE_init      6. cI_init
# 7. cS_Moves     8. cE_Moves     9. cI_Moves
# 10. cS_postM    11. cE_postM    12. cI_postM
# 13. cS_postEI   14. cE_postEI   15. cI_postEI
# 16. cS_postDet  17. cE_postDet  18. cI_postDet
# 19. cS_final    20. cE_final    21. cI_final
# 22. county      23. parish

###~~~ Track ~~~###

# 1. t
# 2. row_id
# 3. pos
# 4. sus_on         5. exp_on        6. inf_on
# 7. sus_off        8. exp_off       9. inf_off
# 10. S_on_out      11. E_on_out     12. I_on_out
# 13. c_new_exp     14. c_new_inf
# 15. number_tested
# 16. number_slaughtered
# 17. E_detected    18. I_detected
# 19. c_birth
# 20. c_S_death     21. c_E_death    22. c_I_death
# 23. number_movements_off
# 24. number_deaths


###~~~ Persistents ~~~###

# 1. t
# 2. row_id
# 3. pos
# 4. c_exp_prob
# 5. c_inf_prob


###~~~ Parish ~~~###

# 1. t
# 2. parish_uid
# 3. parish_pos
# 4. pcS_init        5. pcE_init        6. pcI_init
# 7. pcS_final       8. pcE_final       9. pcI_final
# 10. remaining_pressure
# 11. new_pressure
# 12. scaling
# 13. p_env_prev     14. p_env_cur
# 15. county
# 16. parish
# 17. first_h_of_p_uid

###~~~ Parameters ~~~###

# 1. β
# 2. γ
# 3. F
# 4. ϵ
# 5. ρ
# 6. ρE

###~~~ Record of Movements ~~~###

# 1. t
# 2. off_row_id     3. on_row_id
# 4. cSgen          5. cEgen          6. cIgen
# 7. cS_off_ij      8. cE_off_ij      9. cI_off_ij
# 10. num_moves_ij


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


################################################
### Functions to calculate Deaths likelihood ###
################################################

function cattle_death(;States_postDet, Deaths)

  prob = log_pdf_mvhyper(States_postDet, Deaths)

  return(prob)
end

##################################################################
### Functions to calculate Environmental Reservoir likelihoods ###
##################################################################

function new_pressure(;pI_init, new_pres)

  d = Poisson( pI_init )

  prob = logpdf(d, new_pres)

  return(prob)
end


function remaining_pressure(;scaled_p_env_prev, ϵ, remaining_pres)

  d = Binomial( round(Int, scaled_p_env_prev), (1-ϵ) )

  prob = logpdf(d, remaining_pres)

  return(prob)
end




########################################
### Caclulate LLH Elements Functions ###
########################################

function movements_llh_i_t(;position, t, DATA_res_and_track, movement_records, movement_dict)

  # MOVEMENT PROCESS

  ############
  ### Movmements from outside of farm of interest
  ############

  m_on_out_llh = moves_MHG(States_Moves = [10800, 600, 600],
                                  Moves = DATA_res_and_track[2][position, t, [10,11,12]])

  ############
  #### States of total number of animals moved off farm i at t
  ############

  m_off_llh::Float64 = moves_MHG(States_Moves = DATA_res_and_track[1][position, t, [7,8,9]],
                               Moves = DATA_res_and_track[2][position, t, [7,8,9]])

  ############
  ### States of animals moved to each farm off farm i at t
  ############

  movement_index = movement_dict[(position, t)]

  movement_records_i_t = movement_records[movement_index, :]

  indv_moves_off_llh::Float64 = 0.0

  @inbounds for j in 1:size(movement_records_i_t, 1)

    indv_moves_off_llh = indv_moves_off_llh + Float64(moves_MHG(States_Moves = movement_records_i_t[j, 4:6],
                                                                        Moves = movement_records_i_t[j, 7:9]))

  end

  # [1,2,8]
  llh_res = [m_off_llh, indv_moves_off_llh, m_on_out_llh]

  return(llh_res)
end # ALL MOVEMENTS

function movements_total_llh_i_t(;position, t, DATA_res_and_track, movement_records, movement_dict)

  # MOVEMENT PROCESS

  ############
  #### States of total number of animals moved off farm i at t
  ############

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
                             inf_prob = DATA_pers_and_parish[1][position, t, 5])

  # [3,4]
  return([c_exp_llh, c_inf_llh])
end


function exposures_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # EXPOSURE PROCESS

  c_exp_llh::Float64 = exposures(States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                           new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                            exp_prob = DATA_pers_and_parish[1][position, t, 4])

  return(c_exp_llh::Float64)
end


function infections_llh_i_t(;position, t, DATA_res_and_track, DATA_pers_and_parish)

  # INFECTION PROCESS

  c_inf_llh = infections(;States_postM = DATA_res_and_track[1][position, t, [10,11,12]],
                             new_EandI = DATA_res_and_track[2][position, t, [13,14]],
                              inf_prob = DATA_pers_and_parish[1][position, t, 5])

  return(c_inf_llh)
end


function detection_llh_i_t(;position, t, DATA_res_and_track, epi_params)

  test_occur_i_t = DATA_res_and_track[2][position, t, [17]][1]  #  :test_occur

  # DETECTION PROCESS
    # where test_occur_i_t = 1 if its a test week on farm i, 0 otherwise

  exp_det_llh = test_occur_i_t * exposed_detections(States_postEI = DATA_res_and_track[1][position, t, [13,14,15]],
                                                           dets = DATA_res_and_track[2][position, t, [17, 18]],
                                                   exp_det_prob = (epi_params[5] * epi_params[6]))

  inf_det_llh = test_occur_i_t * infectious_detections(States_postEI = DATA_res_and_track[1][position, t, [13,14,15]],
                                                              dets = DATA_res_and_track[2][position, t, [17, 18]],
                                                      inf_det_prob = epi_params[5])

  # [5,6]
  return([exp_det_llh, inf_det_llh])
end


function c_death_llh_i_t(;position, t, DATA_res_and_track)

  # CATTLE DEATH PROCESS

  c_dth_llh = cattle_death(States_postDet = DATA_res_and_track[1][position, t, [16,17,18]],
                                   Deaths = DATA_res_and_track[2][position, t, [20,21,22]])

  # [7]
  return(c_dth_llh)
end


function p_env_llh_k_t(;p_position, t, DATA_pers_and_parish, epi_params)

  # PARISH ENVIRONMENTAL RESERVOIR

  new_pressure_llh = new_pressure(;pI_init = DATA_pers_and_parish[2][p_position, t, 6],
                                   new_pres = DATA_pers_and_parish[2][p_position, t, 11])


  remaining_pressure_llh = remaining_pressure(;scaled_p_env_prev = prod(DATA_pers_and_parish[2][p_position, t, [12,13]]),
                                               ϵ = epi_params[4],
                                               remaining_pres = DATA_pers_and_parish[2][p_position, t, 10])

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

      llh_array_new[pos, t, [1,2,8]] = movements_llh_i_t(;position = pos, t = t, DATA_res_and_track, movement_records = movement_records, movement_dict = movement_dict)

      llh_array_new[pos, t, 3] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, 4] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [5,6]] = detection_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

      llh_array_new[pos, t, 7] = c_death_llh_i_t(;position = pos, t = t, DATA_res_and_track)

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

      llh_array_new[pos, t, 3] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, 4] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, [5,6]] = detection_llh_i_t(;position = pos, t = t, DATA_res_and_track, epi_params = epi_params)

      llh_array_new[pos, t, 7] = c_death_llh_i_t(;position = pos, t = t, DATA_res_and_track)

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

      llh_array_new[pos, t, 3] = exposures_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

      llh_array_new[pos, t, 4] = infections_llh_i_t(;position = pos, t = t, DATA_res_and_track, DATA_pers_and_parish)

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

@views function update_llh_array_DETECTION(scope::Scope, llh_array_cur, p_env_llh_array_cur, DATA_res_and_track, DATA_pers_and_parish, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

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


### TEST ###

llh_array = zeros(size(DATA_res_and_track[1], 1), 360, 8)

p_env_llh_array = zeros(size(DATA_pers_and_parish[2], 1), 360, 2)

global_scope = Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), Vector(1:8))

llh_array, p_env_llh_array = update_llh_array_ALL(global_scope,
                                                          llh_array, p_env_llh_array,
                                                          DATA_res_and_track, DATA_pers_and_parish, record_of_movements,
                                                          [0.002, 0.015, 0.004, 0.05, 0.75, 0.2], dict_of_movements, f_to_p_structs)

calc_llh_h(global_scope, llh_array)

calc_llh_h_and_p(global_scope, llh_array, p_env_llh_array)



DATA_res_and_track[1][:, 1, :]


DATA_res_and_track[2][:, 1, :]


calc_exp_prop(;States_init, p_env_prev, β, F)



for pos in 1:307
  for t in 1:360
    prob = calc_exp_prop(;States_init = DATA_res_and_track[1][pos, t, 4:6],
                                                     p_env_prev = DATA_pers_and_parish[2][f_to_p_structs[pos].parish_position, t, 13],
                                                     β = epi_params_true[1],
                                                     F = epi_params_true[3])

   if prob != DATA_pers_and_parish[1][pos, t, 4]
     println("calculated probabilty = ", prob)
     println("existing probabilty   = ", DATA_pers_and_parish[1][pos, t, 4])
   end
  end
end

#
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [1]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [2]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [3]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [4]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [5]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [6]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [7]), llh_array)
# calc_llh_h(Scope(1, 360, Vector(1:size(DATA_res_and_track[1], 1)), [8]), llh_array)
#
# which(llh_array[:, :, 3] .== -Inf)
#
# findall(x -> x == -Inf, llh_array[:, :, 3])
#
# llh_array[170, 60:70, :]
