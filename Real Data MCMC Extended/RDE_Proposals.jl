
function rng_mvhyper(n, k)
  # n is a vector of size of pop for each group
  # k is the number of trials (total number moving)

  N = sum(n) # total pop size
  m = length(n) # number of groups
  n_otr = N - n[1] # number not in first group

  x = fill(0, m) # results

  x[1] = rand(Hypergeometric(n[1], n_otr, k))

  for i in 2:(m-1)
    n_otr = n_otr - n[i]
    k = k - x[i-1]
    x[i] = rand(Hypergeometric(n[i], n_otr, k))
  end

  x[m] = k - x[m-1]
  return(x)
end

###############
### Move SE ###
###############

function propose_Move_SE(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_SE_events_at_t = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[2][:, :, 13] .> 0), [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_SE_events_at_t[rand(1:size(farms_with_SE_events_at_t, 1)), :]

  # println("  ", "t = ", t)
  # println("  ", "position = ", position)

  Δ = 0
  while Δ == 0
    # Δ = rand(Poisson(3))
    Δ = 1
    sgn_Δ = rand([-1, 1])
    Δ = sgn_Δ * Δ
  end

  # println("   ", "t + Δ = ", t + Δ)

  if 360 < t + Δ || t + Δ <= 0
    # println("   ", "Out of bounds! ")
    log_q_ratio = -Inf
    Move_SE_track = [position, t, 0., 4., Δ, 1.]
                   # :mSE_position, :mSE_t, :mSE_is_accepted, :mSE_reason, :mSE_Δ_time, :mSE_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Move_SE_track)
  end

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_Move_SE(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, 1, f_to_p_structs)
  # num_SE_moved set to 1

  if valid != 1
    log_q_ratio = -Inf
    Move_SE_track = [position, t, 0., 2., Δ, 1.]
                   # :mSE_position, :mSE_t, :mSE_is_accepted, :mSE_reason, :mSE_Δ_time, :mSE_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Move_SE_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_SE_events_at_t = size(farms_with_SE_events_at_t, 1)

  num_farms_with_SE_events_at_t_prime = size(DATA_res_and_track_prime[2][(DATA_res_and_track_prime[2][:, :, 13] .> 0), [1,3]], 1)

  # println("   ", "num_farms_with_SE_events_at_t = ", num_farms_with_SE_events_at_t)
  # println("   ", "num_farms_with_SE_events_at_t_prime = ", num_farms_with_SE_events_at_t_prime)

  log_q_ratio = log(num_farms_with_SE_events_at_t / num_farms_with_SE_events_at_t_prime)


  Move_SE_track = [position, t, 0., 0., Δ, 1.]
                 # :mSE_position, :mSE_t, :mSE_is_accepted, :mSE_reason, :mSE_Δ_time, :mSE_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Move_SE_track)
end


###############
### Move EI ###
###############

function propose_Move_EI(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_EI_events_at_t = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[2][:, :, 14] .> 0), [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_EI_events_at_t[rand(1:size(farms_with_EI_events_at_t, 1)), :]

  Δ = 0
  while Δ == 0
    Δ = rand(Poisson(3))
    sgn_Δ = rand([-1, 1])
    Δ = sgn_Δ * Δ
  end

  if 360 < t + Δ || t + Δ <= 0
    log_q_ratio = -Inf
    Move_EI_track = [position, t, 0., 4., Δ, 1.]
                   # :mEI_position, :mEI_t, :mEI_is_accepted, :mEI_reason, :mEI_Δ_time, :mEI_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Move_EI_track)
  end

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_Move_EI(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, 1, epi_params, f_to_p_structs)
  # num_EI_moved set to 1

  if valid != 1
    log_q_ratio = -Inf
    Move_EI_track = [position, t, 0., 2., Δ, 1.]
                  # :mEI_position, :mEI_t, :mEI_is_accepted, :mEI_reason, :mEI_Δ_time, :mEI_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Move_EI_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_EI_events_at_t = size(farms_with_EI_events_at_t, 1)

  num_farms_with_EI_events_at_t_prime = size(DATA_res_and_track_prime[2][(DATA_res_and_track_prime[2][:, :, 14] .> 0), [1,3]], 1)

  log_q_ratio = log(num_farms_with_EI_events_at_t / num_farms_with_EI_events_at_t_prime)


  Move_EI_track = [position, t, 0., 0., Δ, 1.]
                 # :mEI_position, :mEI_t, :mEI_is_accepted, :mEI_reason, :mEI_Δ_time, :mEI_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Move_EI_track)
end


#################
### AddRem SE ###
#################

function propose_AddRem_SE(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_S_and_exp_prob_at_t = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, :, 10] .> 0  .&& DATA_pers_and_parish_cur[1][:, :, 4] .> 0), [1,3]]
  # Extract all farms that have cS_postM_t > 0 AND prob_exp_t > 0

  ### Generate the update parameters ###

  t, position = farms_with_S_and_exp_prob_at_t[rand(1:size(farms_with_S_and_exp_prob_at_t, 1)), :]

  Δ = rand([-1, 1])

  AddRem_SE_track = [position, t, 0., 0., 1., -999., -999., -999., -999.]
                   # :arSE_position, :arSE_t, :arSE_is_accepted, :arSE_reason, :arSE_Δ, :arSE_SE_before, :arSE_SE_after, :arSE_cS, :arSE_prob

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_SE_track = update_data_AddRem_SE(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, f_to_p_structs, AddRem_SE_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_SE_track[4] = 2
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_SE_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_S_and_exp_prob_at_t = size(farms_with_S_and_exp_prob_at_t, 1)

  num_farms_with_S_and_exp_prob_at_t_prime = size(DATA_res_and_track_prime[2][(DATA_res_and_track_prime[1][:, :, 10] .> 0  .&& DATA_pers_and_parish_prime[1][:, :, 4] .> 0), [1,3]], 1)

  log_q_ratio = log(num_farms_with_S_and_exp_prob_at_t / num_farms_with_S_and_exp_prob_at_t_prime)


  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_SE_track)
end


#################
### AddRem EI ###
#################

function propose_AddRem_EI(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_E_at_t = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, :, 11] .> 0), [1,3]]
  # Extract all farms that have cE_postM_t > 0

  ### Generate the update parameters ###

  t, position = farms_with_E_at_t[rand(1:size(farms_with_E_at_t, 1)), :]

  Δ = rand([-1, 1])

  AddRem_EI_track = [position, t, 0., 0., 1., -999., -999., -999., -999.]
                   # :arEI_position, :arEI_t, :arEI_is_accepted, :arEI_reason, :arEI_Δ, :arEI_EI_before, :arEI_EI_after, :arEI_cE, :arEI_prob

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_EI_track = update_data_AddRem_EI(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs, AddRem_EI_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_EI_track[4] = 2

    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_EI_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_E_at_t = size(farms_with_E_at_t, 1)

  num_farms_with_E_at_t_prime = size(DATA_res_and_track_prime[2][(DATA_res_and_track_prime[1][:, :, 11] .> 0), [1,3]], 1)

  log_q_ratio = log(num_farms_with_E_at_t / num_farms_with_E_at_t_prime)


  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_EI_track)
end


#########################
### AddRem Detections ###
#########################

function detection_permutations(cStates_postEI, cDet_cur, epi_params)

  total_det_t = sum(cDet_cur)

  permutations = fill(-99., Int64(factorial(big(total_det_t+2-1))/(factorial(big(total_det_t))*factorial(2-1))), 3)

  iter = 1

  for nE in 0:total_det_t

    nI = (total_det_t - nE)

    prob_Edet = logpdf(Binomial(convert(Int64, cStates_postEI[2]), prod(epi_params[[5,6]])), nE)
    prob_Idet = logpdf(Binomial(convert(Int64, cStates_postEI[3]), epi_params[5]), nI)

    permutations[iter, :] = [nE, nI, (exp(prob_Edet) + exp(prob_Idet))]
    iter += 1
  end

  chosen_perm = [permutations[wsample(1:size(permutations, 1), permutations[:, 3]), :] ; permutations[convert(Int64, (cDet_cur[1]+1)), 3]]

  return(chosen_perm)
end

function propose_AddRem_Det(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_permutable_detections_at_t = DATA_res_and_track_cur[2][((DATA_res_and_track_cur[1][:, :, 14] .> 0 .&& DATA_res_and_track_cur[1][:, :, 15] .> 0) .&&
                                                   (DATA_res_and_track_cur[2][:, :, 17] .> 0) .|| (DATA_res_and_track_cur[2][:, :, 18] .> 0)),
                                                  [1,3]]
  # Extract all farms that have cE_postEI_t > 0 AND cI_postEI_t > 0 AND (Edet_t >0 OR I_det_t > 0)


  ### Generate the update parameters ###

  t, position = farms_with_permutable_detections_at_t[rand(1:size(farms_with_permutable_detections_at_t, 1)), :]

  cStates_postEI_t = Array(DATA_res_and_track_cur[1][position, t, 13:15])
  detections_t = Array(DATA_res_and_track_cur[2][position, t, 17:18])

  # detections_t_new = rng_mvhyper(cStates_postEI_t[2:3], sum(detections_t))
  # NOTE: The probabilities of an E detection event is not that same as an I detection event
  # so we may need an alternate method that accounts for this:

  detections_t_new_and_probs = detection_permutations(cStates_postEI_t, detections_t, epi_params)

  Δs = detections_t_new_and_probs[1:2] - detections_t

  AddRem_Det_track = [position, t, 0., 0., Δs[1], Δs[2], -999., -999., -999., -999., -999., -999.]
                  # :arDet_position, :arDet_t, :arDet_is_accepted, :arDet_reason, :arDet_ΔE, :arDet_ΔI,
                                   # :arDet_Edet_before, :arDet_Idet_before, :arDet_Edet_after, :arDet_Idet_after,
                                   # :arDet_cE, :arDet_cI
  if Δs[1] == 0
    log_q_ratio = -Inf
    AddRem_Det_track[4] = 3
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Det_track)
  end


  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_Det_track = update_data_AddRem_Det(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δs, epi_params, f_to_p_structs, AddRem_Det_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_Det_track[4] = 2
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Det_track)
  end


  ### Calculate the log q ratio (proposal density ratio) ###
  log_q_ratio = detections_t_new_and_probs[4] - detections_t_new_and_probs[3] # q_cur_given_prime - q_prime_given_cur

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_Det_track)
end


#####################
### AddRem Deaths ###
#####################

# function deaths_permutations(cStates_postDet, cDeaths, epi_params)
#
#   total_deaths = sum(cDeaths)
#
#   permutations = fill(-99., Int64(factorial(total_deaths+3-1)/(factorial(total_deaths)*factorial(3-1))), 4)
#
#   iter = 1
#
#   for nS in 0:total_deaths
#     for nE in 0:(total_deaths - nS)
#
#       nI = (total_deaths - nS - nE)
#
#       perm_prob = log_pdf_mvhyper(cStates_postDet, [nS, nE, nI])
#
#       permutations[iter, :] = [nS, nE, nI, exp(perm_prob)]
#       iter += 1
#     end
#   end
#
#   chosen_perm = [permutations[wsample(1:size(permutations, 1), permutations[:, 4]), :] ; log_pdf_mvhyper(cStates_postDet, cDeaths)]
#
#   return(chosen_perm)
# end

function propose_AddRem_Deaths(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_permutable_deaths_at_t = DATA_res_and_track_cur[2][( ((DATA_res_and_track_cur[1][:, :, 16] .> 0 .&& DATA_res_and_track_cur[1][:, :, 17] .> 0) .||
                                                            (DATA_res_and_track_cur[1][:, :, 16] .> 0 .&& DATA_res_and_track_cur[1][:, :, 18] .> 0) .||
                                                            (DATA_res_and_track_cur[1][:, :, 17] .> 0 .&& DATA_res_and_track_cur[1][:, :, 18] .> 0)) .&&
                                                            (DATA_res_and_track_cur[2][:, :, 20] .> 0 .|| DATA_res_and_track_cur[2][:, :, 21] .> 0 .|| DATA_res_and_track_cur[2][:, :, 22] .> 0) ),
                                                          [1,3]]
  # Extract all farms that have atleast 2 cStates and at least 1 death.


  ### Generate the update parameters ###

  t, position = farms_with_permutable_deaths_at_t[rand(1:size(farms_with_permutable_deaths_at_t, 1)), :]

  cStates_postDet_t = Array(DATA_res_and_track_cur[1][position, t, 16:18])
  deaths_t = Array(DATA_res_and_track_cur[2][position, t, 20:22])

  deaths_t_new = rng_mvhyper(cStates_postDet_t, sum(deaths_t))

  Δs = deaths_t_new - deaths_t

  AddRem_Deaths_track = [position, t, 0., 0., Δs[1], Δs[2], Δs[3], -999., -999., -999., -999., -999., -999., -999., -999., -999.]
                      # :arDeaths_position, :arDeaths_t, :arDeaths_is_accepted, :arDeaths_reason,
                                          # :arDeaths_ΔS, :arDeaths_ΔE, :arDeaths_ΔI,
                                          # :arDeaths_Sdths_before, :arDeaths_Edths_before, :arDeaths_Idths_before,
                                          # :arDeaths_Sdths_after, :arDeaths_Edths_after, :arDeaths_Idths_after,
                                          # :arDeaths_cS, :arDeaths_cE, :arDeaths_cI

  if sum(Δs) == 0
    log_q_ratio = -Inf
    AddRem_Deaths_track[4] = 3
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Deaths_track)
  end


  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_Deaths_track = update_data_AddRem_Deaths(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δs, epi_params, f_to_p_structs, AddRem_Deaths_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_Deaths_track[4] = 2
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Deaths_track)
  end

  ### Calculate the log q ratio (proposal density ratio) ###
  log_q_ratio = log_pdf_mvhyper(cStates_postDet, deaths_t) - log_pdf_mvhyper(cStates_postDet, deaths_t_new) # q_cur_given_prime - q_prime_given_cur

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_Deaths_track)
end


########################
### AddRem Movements ###
########################

function log_q_ratio_movements(log_q_ratio_data)

  q_prime_given_cur = 0.0
  q_cur_given_prime = 0.0

  idx_size = Int(size(log_q_ratio_data, 1)/2)

  for i in 1:idx_size
    q_prime_given_cur += log_pdf_mvhyper(log_q_ratio_data[i, 1], log_q_ratio_data[i, 2])
    q_cur_given_prime += log_pdf_mvhyper(log_q_ratio_data[(i+idx_size), 1], log_q_ratio_data[(i+idx_size), 1])
  end

  log_q_ratio = q_cur_given_prime - q_prime_given_cur

  return(log_q_ratio)
end

function generate_new_movement(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, epi_params, movement_record, movement_dict, f_to_p_structs::Vector{Farm_Parish_info}, ids_to_pos_dict, tracker)

  DATA_res_and_track_prime = deepcopy(DATA_res_and_track_cur)
  DATA_pers_and_parish_prime = deepcopy(DATA_pers_and_parish_cur)


  movement_record_prime = deepcopy(movement_record)

  ### Scope ###

  lower_t = t
  upper_t = size(DATA_res_and_track_cur[1], 2) #T
  changed_h_positions = [position] # extended dependent on movement choice
  h_llh_indices = 1:8

  ####################
  ### Extract Data ###
  ####################

  # :cS_Moves, :cE_Moves, :cI_Moves
  states = DATA_res_and_track_prime[1][position, t, [7,8,9]]
  # :sus_off, :exp_off, :inf_off
  moves_off = DATA_res_and_track_prime[2][position, t, [7,8,9]]

  total_moves = sum(moves_off)


  ############################################
  ### Extract the movement data for time t ###
  ############################################

  move_record_rows = movement_dict[(position, t)]

  # This extracts the movement data at time t for the off farm
  moves_data_cur = movement_record[move_record_rows, :]
  moves_data_prime = deepcopy(moves_data_cur)

  # Extract on_row_ids that recieved animals from this farm
  on_cph_row_ids = moves_data_prime[:, 3]

  # Extract number of movements to each farm
  m_farm = moves_data_prime[:, 10]


  ##################################
  ### Generate new movement data ###
  ##################################

  new_move_off = rng_mvhyper(states, total_moves)

  running_move_off = deepcopy(new_move_off)

  # FOR each (j) of the farms that recieved animals
  for j in 1:size(on_cph_row_ids, 1)

    # IF there is more than one farm that recieved animals from farm i
    if (size(on_cph_row_ids, 1) > 1)
      # Extract the total movements onto this farm
      total_moved = m_farm[j]
      # Generate the states
      states_on = rng_mvhyper(running_move_off, total_moved)
    else # ELSE there is only one farm
      states_on = running_move_off
    end

    moves_data_prime[j, 4:6] = running_move_off
    moves_data_prime[j, 7:9] = states_on

    # Update the running total of the moved states
    # Reduce the moved states by those that just got assigned
    running_move_off = running_move_off - states_on

  end #end of for farms that recieved animals


  ########################
  ### Update moves off ###
  ########################

  # :sus_off, :exp_off, :inf_off
  DATA_res_and_track_prime[2][position, t, [7,8,9]] = new_move_off

  movement_record_prime[move_record_rows, :] .= moves_data_prime


  ############################################
  ### Calculate the farm level differences ###
  ############################################

  # :Δ_S_on, :Δ_E_on, :Δ_I_on, :on_row_id, :on_position, :on_parish_pos
  differences = fill(-99, size(on_cph_row_ids, 1), 6)

  for j in 1:size(on_cph_row_ids, 1)

      if on_cph_row_ids[j] > 0
        differences[j, 1:6] = [ Array(moves_data_cur[j, 7:9]) - Array(moves_data_prime[j, 7:9]) ;
                                on_cph_row_ids[j] ;
                                ids_to_pos_dict[on_cph_row_ids[j]] ;
                                f_to_p_structs[ ids_to_pos_dict[ on_cph_row_ids[j] ] ].parish_position ]


        changed_h_positions = [changed_h_positions ; differences[j, 5]]
      end
  end

  changed_h_positions_oi = unique(changed_h_positions[changed_h_positions .> 0])

  differences_oi = differences[(differences[:,4] .> 0), :]

  ##############################################
  ### Calculate the parish level differences ###
  ##############################################

  affected_parishes = unique([differences_oi[:, 6]; f_to_p_structs[position].parish_position ])

  # :Δ_pS, :Δ_pE, :Δ_pI, :parish_pos
  parish_differences = fill(0, size(affected_parishes, 1), 4)

  for (idx, par) in enumerate(affected_parishes)

    parish_differences[idx, 4] = par

    for i in 1:size(differences_oi, 1)
      if differences_oi[i, 5] == par
        parish_differences[idx, 1:3] += differences_oi[i, 1:3]
      end
    end

    if par == f_to_p_structs[position].parish_position #off_p_id
      Δ_off = DATA_res_and_track_prime[2][position, t, [7,8,9]] - DATA_res_and_track_cur[2][position, t, [7,8,9]]
      parish_differences[idx, 1:3] -= Δ_off
    end

  end

  ###############################
  ### Prepare for log q ratio ###
  ###############################

  log_q_ratio_data = fill([-99,-99,-99], 2*(size(on_cph_row_ids, 1)+1), 2)

  temp_moves_cur = Array(moves_data_cur)[:, 4:9]
  temp_moves_prime = Array(moves_data_prime)[:, 4:9]

  size_moves_cur = size(moves_data_cur, 1)

  log_q_ratio_data[1, :] = [states, moves_off]
  log_q_ratio_data[(2+size_moves_cur), :] = [states, new_move_off]

  tracker[5:13] = [moves_off ; new_move_off ; states]

  for k in 1:size_moves_cur
    log_q_ratio_data[(1+k), :] = [temp_moves_cur[k, 1:3], temp_moves_cur[k, 4:6]]
    log_q_ratio_data[(2+size_moves_cur+k), :] = [temp_moves_prime[k, 1:3], temp_moves_prime[k, 4:6]]
  end

  scope = Scope(lower_t, upper_t, changed_h_positions_oi, h_llh_indices)

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, movement_record_prime, scope, differences_oi, parish_differences, log_q_ratio_data, tracker)
end

function propose_AddRem_Movements(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, movement_record, movement_dict, f_to_p_structs::Vector{Farm_Parish_info}, ids_to_pos_dict)

  ### Choose a farm and a timestep ###

  farms_with_permutable_movements_at_t = DATA_res_and_track_cur[2][( ((DATA_res_and_track_cur[1][:, :, 7] .> 0 .&& DATA_res_and_track_cur[1][:, :, 8] .> 0) .||
                                                              (DATA_res_and_track_cur[1][:, :, 7] .> 0 .&& DATA_res_and_track_cur[1][:, :, 9] .> 0) .||
                                                              (DATA_res_and_track_cur[1][:, :, 8] .> 0 .&& DATA_res_and_track_cur[1][:, :, 9] .> 0)) .&&
                                                              (DATA_res_and_track_cur[2][:, :, 7] .> 0 .|| DATA_res_and_track_cur[2][:, :, 8] .> 0 .|| DATA_res_and_track_cur[2][:, :, 9] .> 0) ),
                                                            [1,3]]
  # Extract all farms that have atleast 2 cStates and at least 1 movement.


  ### Generate the update parameters ###

  t, position = farms_with_permutable_movements_at_t[rand(1:size(farms_with_permutable_movements_at_t, 1)), :]


  AddRem_Movements_track = [position, t, 0., 0., -999., -999., -999., -999., -999., -999., -999., -999., -999.]
                      # :arMoves_position, :arMoves_t, :arMoves_is_accepted, :arMoves_reason,
                                         # :arMoves_Soff_before, :arMoves_Eoff_before, :arMoves_Ioff_before,
                                         # :arMoves_Soff_after, :arMoves_Eoff_after, :arMoves_Ioff_after,
                                         # :arMoves_cS, :arMoves_cE, :arMoves_cI

  ### Generate a new set of movements ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, movement_record_prime, scope, differences_oi, parish_differences, log_q_ratio_data, AddRem_Movements_track = generate_new_movement(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, epi_params, movement_record,
                                                                                                                            movement_dict, f_to_p_structs, ids_to_pos_dict, AddRem_Movements_track)

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, valid = update_data_AddRem_Movement(DATA_res_and_track_cur, DATA_pers_and_parish_cur, scope, position, DATA_res_and_track_prime, DATA_pers_and_parish_prime, epi_params, differences_oi, parish_differences)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_Movements_track[4] = 2
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, scope, movement_record, AddRem_Movements_track)
  end

  ### Calculate the log q ratio (proposal density ratio) ###
  log_q_ratio = log_q_ratio_movements(log_q_ratio_data)

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, movement_record_prime, AddRem_Movements_track)
end

function propose_AddRem_outside_move_states(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ##################
  ### Choose a farm and a timestep
  ##################

  farms_with_permutable_movements_at_t = DATA_res_and_track_cur[2][ (DATA_res_and_track_cur[2][:, :, 10] .> 0 .|| DATA_res_and_track_cur[2][:, :, 11] .> 0 .|| DATA_res_and_track_cur[2][:, :, 12] .> 0), [1,3]]

  # Extract all farms that have at least 1 movement on from outside


  ##################
  ### Generate the update parameters
  ##################

  t, position = farms_with_permutable_movements_at_t[rand(1:size(farms_with_permutable_movements_at_t, 1)), :]


  ##################
  ### Generate a new set of movements
  ##################

  states_on_out_t = DATA_res_and_track_cur[2][position, t, 10:12]

  states_on_out_t_new = rng_mvhyper([10800,600,600], sum(states_on_out_t))

  Δs = states_on_out_t_new - states_on_out_t



  AddRem_Movements_outside_track = [position, t, 0., 0., Δs[1], Δs[2], Δs[3], -999., -999., -999., -999., -999., -999.]
                     # :arMoves_out_position, :arMoves_out_t, :arMoves_out_is_accepted, :arMoves_out_reason,
                                         # :arMoves_out_ΔS, :arMoves_out_ΔE, :arMoves_out_ΔI,
                                         # :arMoves_out_S_on_out_before, :arMoves_out_E_on_out_before, :arMoves_out_I_on_out_before,
                                         # :arMoves_out_S_on_out_after, :arMoves_out_E_on_out_after, :arMoves_out_I_on_out_after



   if sum(Δs) == 0
      log_q_ratio = -Inf
      AddRem_Movements_outside_track[4] = 3
      return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Movements_outside_track)
   end

  ##################
  ### Calculate the update
  ##################

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_Movements_outside_track = update_data_AddRem_outside_move_states(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δs, epi_params, f_to_p_structs, AddRem_Movements_outside_track)

  if valid != 1
     log_q_ratio = -Inf
     AddRem_Movements_outside_track[4] = 2
     return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_Movements_outside_track)
  end

  ##################
  ### Calculate the log q ratio (proposal density ratio)
  ##################

  log_q_ratio = log_pdf_mvhyper([10800,600,600], states_on_out_t) - log_pdf_mvhyper([10800,600,600], states_on_out_t_new) # q_cur_given_prime - q_prime_given_cur

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_Movements_outside_track)
end

#####################################
### AddRem Environmental Pressure ###
#####################################

function propose_AddRem_penv(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_env_pressure_at_t = DATA_pers_and_parish_cur[2][( (DATA_pers_and_parish_cur[2][:, :, 6] .> 0  .||
                                                     DATA_pers_and_parish_cur[2][:, :, 13] .> 0) .&&
                                                     DATA_pers_and_parish_cur[2][:, :, 1] .< 360 ),
                                                   [1,3]]
  # Extract all farms that have p_env_prev_t > 0 or pcI_t > 0

  ### Generate the update parameters ###

  t, p_position = Int.(farms_with_env_pressure_at_t[rand(1:size(farms_with_env_pressure_at_t, 1)), :])

  ### Generate the update ###

  penv_data = DATA_pers_and_parish_cur[2][p_position, t, [6,10,11,12,13]]
  # pcI_t,   remaining_pressure_t, new_pressure_t,    scaling, p_env_prev_t

  ϵ = epi_params[4]

  propose_r_pressure_t = rand(Binomial(round(Int, penv_data[4] * penv_data[5]), 1-ϵ))
  propose_n_pressure_t = rand(Poisson(sum(penv_data[1])))

  Δs = [propose_r_pressure_t - penv_data[2],  propose_n_pressure_t - penv_data[3]]

  AddRem_penv_track = [p_position, t, 0., 0., Δs[1], Δs[2], -999., -999., -999., -999., -999., -999.]
                      # :arpenv_position, :arpenv_t, :arpenv_is_accepted, :arpenv_reason, :arpenv_Δr, :arpenv_Δn,
                                        # :arpenv_r_pres_before, :arpenv_n_pres_before, :arpenv_r_pres_after, :arpenv_n_pres_after,
                                        # :arpenv_p_env_prev, :arpenv_pI

  if sum(Δs) == 0
    log_q_ratio = -Inf
    AddRem_penv_track[4] = 3

    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_penv_track)
  end


  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_penv_track = update_data_AddRem_penv(DATA_res_and_track_cur, DATA_pers_and_parish_cur, p_position, t, Δs, epi_params, f_to_p_structs, ids_to_pos_dict, AddRem_penv_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_penv_track[4] = 2

    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_penv_track)
  end


  ### Calculate log q ratio ###

  q_prime_r_press_given_cur = logpdf(Binomial(round(Int, penv_data[4] * penv_data[5]), 1-ϵ), propose_r_pressure_t)
  q_prime_n_press_given_cur = logpdf(Poisson(sum(penv_data[1])), propose_n_pressure_t)

  q_cur_r_press_given_prime = logpdf(Binomial(round(Int, penv_data[4] * penv_data[5]), 1-ϵ), penv_data[2])
  q_cur_n_press_given_prime = logpdf(Poisson(sum(penv_data[1])), penv_data[3])

  log_q_ratio = (q_cur_r_press_given_prime + q_cur_n_press_given_prime) - (q_prime_r_press_given_cur + q_prime_n_press_given_cur)

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_penv_track)
end


####################################
### PROPOSE PARAMETERS (Generic) ###
####################################

function accepted_prop(n_tune, other_res)

  ## Calculate acceptance rate of last batch ##
  batch_start = ((n_tune - 1) * 25) + 1
  batch_end = (n_tune*25) - 1

  # Accepted Proportion
  return( sum(other_res[batch_start:batch_end]) / 24 )
end

function calculate_delta_n(n_tune, acc_prop)

  if acc_prop < 0.33
    delta_n = -1 * min(0.05, 1/sqrt(n_tune))
  else
    delta_n = min(0.05, 1/sqrt(n_tune))
  end

  return(delta_n)
end

function tune_λ(λ_cur, n_tune, other_res)

    acc_prop = accepted_prop(n_tune, other_res)

    ## Update tuning parameter ##
    delta_n = calculate_delta_n(n_tune, acc_prop)

    # log(λ_new) = log(λ_old) + Δn
    # λ_new = exp( log(λ_old) + Δn )
    # λ_new = exp(log(λ_old)) * exp(Δn)
  return(λ_cur * exp(delta_n))
end

function mixture_less(d, λ, params_cur)

  params_draw = rand(MvNormal(params_cur, ( (1/d) * (λ)^2 * I ) ))

  return(params_draw)
end

function mixture_more(covarM, m, params_cur)

  params_draw = rand(MvNormal(params_cur, ((m)^2 * cov(covarM)) ))

  return(params_draw)
end

function propose_params(N_its, results, other_res,
                        it, log_params_cur, oi, accepted_idx,
                        n_tune, m, λ, d, covarM)

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

    # log_params_draw = deepcopy(log_params_cur)

  ######################
  ### Initial tuning ###
  ######################

    if it < min(5000, (N_its/10))

      if it == (n_tune*25)
        λ = tune_λ(λ, n_tune, other_res[:, accepted_idx])
      end

      draw = mixture_less(3, λ, log_params_cur[oi])

      mixture = 0.0
    end

  ##########################
  ### Draw the paramters ###
  ##########################

    if it >= min(5000, (N_its/10))

      mixture = rand()

      ### Non-adaptive draw ###

      if mixture < 0.05
        draw = mixture_less(d, λ, log_params_cur[oi])
      end

      ### Adaptive Draw ###

      if mixture >= 0.05
        draw = mixture_more(covarM, m, log_params_cur[oi])
      end

    end

    log_params_draw = deepcopy(log_params_cur)
    log_params_draw[oi] = draw

  #log_params_draw, log_q_ratio, mixture, λ
  return(log_params_draw, 0.0, mixture, λ)
end

###################################
### PROPOSE EPIDEMIC PARAMETERS ###
###################################

function propose_epidemic_params(N_its, results, other_res,
                                  it, log_params_cur, DATA_res_and_track_cur, DATA_pers_and_parish_cur, f_to_p_structs::Vector{Farm_Parish_info},
                                  n_tune, m, λ, d, covarM)

  ##########################
  ### Propose parameters ###
  ##########################

    log_params_draw, log_q_ratio, mixture, λ = propose_params(N_its, results, other_res,
                                                          it, log_params_cur, [1,2,3,4], 1,
                                                          n_tune, m, λ, d, covarM)

  #######################
  ### Update the data ###
  #######################

    scope = Scope(1, 360, Vector(1:size(DATA_res_and_track_cur[1], 1)), [3,4])
    # lower_t, upper_t, h_pos_ids, h_llh_indicies

    DATA_res_and_track_prime, DATA_pers_and_parish_prime = update_pers_EPIDEMIC(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_params_draw, f_to_p_structs, scope)

  return(log_params_draw, 0.0, mixture, λ, DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope)
end

####################################
### PROPOSE DETECTION PARAMETERS ###
####################################

function propose_detection_params(N_its, results, other_res,
                                  it, log_params_cur, DATA_res_and_track_cur, DATA_pers_and_parish_cur, f_to_p_structs,
                                  n_tune, m, λ, d, covarM)

  ##########################
  ### Propose parameters ###
  ##########################

    log_params_draw, log_q_ratio, mixture, λ = propose_params(N_its, results, other_res,
                                                              it, log_params_cur, [5,6], 6,
                                                              n_tune, m, λ, d, covarM)

    scope = Scope(1, 360, Vector(1:size(DATA_res_and_track_cur[1], 1)), [5,6])
    # lower_t, upper_t, h_pos_ids, h_llh_indices

  return(log_params_draw, 0.0, mixture, λ, DATA_res_and_track_cur, DATA_pers_and_parish_cur, scope)
end



#################################################
### CHANGE INITIAL CONDITIONS (CHANGE S TO E) ###
#################################################

function propose_change_init_S_to_E(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_S_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 4] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_S_init[rand(1:size(farms_with_S_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_S_to_E(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_S_to_E_track = [position, t, 0., 2., Δ]
                   # :ciSE_position, :ciSE_t, :ciSE_is_accepted, :ciSE_reason, :ciSE_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_S_to_E_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_S_init = size(farms_with_S_init, 1)

  num_farms_with_E_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 5] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_S_init / num_farms_with_E_init_prime)


  Init_S_to_E_track = [position, t, 0., 0., Δ]
                 # :ciSE_position, :ciSE_t, :ciSE_is_accepted, :ciSE_reason, :ciSE_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_S_to_E_track)
end

#################################################
### CHANGE INITIAL CONDITIONS (CHANGE S TO I) ###
#################################################

function propose_change_init_S_to_I(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_S_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 4] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_S_init[rand(1:size(farms_with_S_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_S_to_I(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_S_to_I_track = [position, t, 0., 2., Δ]
                   # :ciSI_position, :ciSI_t, :ciSI_is_accepted, :ciSI_reason, :ciSI_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_S_to_I_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_S_init = size(farms_with_S_init, 1)

  num_farms_with_I_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 6] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_S_init / num_farms_with_I_init_prime)


  Init_S_to_I_track = [position, t, 0., 0., Δ]
                 # :ciSI_position, :ciSI_t, :ciSI_is_accepted, :ciSI_reason, :ciSI_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_S_to_I_track)
end

#################################################
### CHANGE INITIAL CONDITIONS (CHANGE E TO S) ###
#################################################

function propose_change_init_E_to_S(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_E_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 5] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_E_init[rand(1:size(farms_with_E_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_E_to_S(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_E_to_S_track = [position, t, 0., 2., Δ]
                   # :ciES_position, :ciES_t, :ciES_is_accepted, :ciES_reason, :ciES_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_E_to_S_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_E_init = size(farms_with_E_init, 1)

  num_farms_with_S_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 4] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_E_init / num_farms_with_S_init_prime)


  Init_E_to_S_track = [position, t, 0., 0., Δ]
                 # :ciES_position, :ciES_t, :ciES_is_accepted, :ciES_reason, :ciES_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_E_to_S_track)
end

#################################################
### CHANGE INITIAL CONDITIONS (CHANGE E TO I) ###
#################################################

function propose_change_init_E_to_I(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_E_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 5] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_E_init[rand(1:size(farms_with_E_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_E_to_I(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_E_to_I_track = [position, t, 0., 2., Δ]
                   # :ciEI_position, :ciEI_t, :ciEI_is_accepted, :ciEI_reason, :ciEI_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_E_to_I_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_E_init = size(farms_with_E_init, 1)

  num_farms_with_I_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 6] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_E_init / num_farms_with_I_init_prime)


  Init_E_to_I_track = [position, t, 0., 0., Δ]
                 # :ciEI_position, :ciEI_t, :ciEI_is_accepted, :ciEI_reason, :ciEI_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_E_to_I_track)
end

#################################################
### CHANGE INITIAL CONDITIONS (CHANGE I TO S) ###
#################################################

function propose_change_init_I_to_S(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_I_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 6] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_I_init[rand(1:size(farms_with_I_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_I_to_S(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_I_to_S_track = [position, t, 0., 2., Δ]
                   # :ciIS_position, :ciIS_t, :ciIS_is_accepted, :ciIS_reason, :ciIS_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_I_to_S_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_I_init = size(farms_with_I_init, 1)

  num_farms_with_S_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 4] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_I_init / num_farms_with_S_init_prime)


  Init_I_to_S_track = [position, t, 0., 0., Δ]
                 # :ciIS_position, :ciIS_t, :ciIS_is_accepted, :ciIS_reason, :ciIS_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_I_to_S_track)
end

#################################################
### CHANGE INITIAL CONDITIONS (CHANGE I TO E) ###
#################################################

function propose_change_init_I_to_E(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  farms_with_I_init = DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 6] .> 0), 1, [1,3]]

  ### Generate the update parameters ###

  t, position = farms_with_I_init[rand(1:size(farms_with_I_init, 1)), :]

  Δ = 1

  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid = update_data_change_init_I_to_E(DATA_res_and_track_cur, DATA_pers_and_parish_cur, position, t, Δ, epi_params, f_to_p_structs)

  if valid != 1
    log_q_ratio = -Inf
    Init_I_to_E_track = [position, t, 0., 2., Δ]
                   # :ciIE_position, :ciIE_t, :ciIE_is_accepted, :ciIE_reason, :ciIE_num_moved
    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], Init_I_to_E_track)
  end


  ### Log q ratio calculation ###

  num_farms_with_I_init = size(farms_with_I_init, 1)

  num_farms_with_E_init_prime = size(DATA_res_and_track_cur[2][(DATA_res_and_track_cur[1][:, 1, 5] .> 0), 1, [1,3]], 1)

  log_q_ratio = log(num_farms_with_I_init / num_farms_with_E_init_prime)


  Init_I_to_E_track = [position, t, 0., 0., Δ]
                 # :ciSE_position, :ciSE_t, :ciSE_is_accepted, :ciSE_reason, :ciSE_num_moved

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, Init_I_to_E_track)
end

####################################################
### CHANGE INITIAL CONDITIONS (CHANGE penv INIT) ###
####################################################

function propose_change_init_penv(DATA_res_and_track_cur, DATA_pers_and_parish_cur, epi_params, f_to_p_structs::Vector{Farm_Parish_info})

  ### Choose a farm and a timestep ###

  t, p_position = [1, rand(1:size(DATA_pers_and_parish_cur[2], 1))]

  ### Generate the update ###

  Δ = rand(1)[1] * rand([-1,1]) #change in p_env_init

  AddRem_penv_track = [p_position, t, 0., 0., Δ, -999., -999.]
                      # :arpenv_position, :arpenv_t, :arpenv_is_accepted, :arpenv_reason, :arpenv_Δ,
                                        # :arpenv_p_env_prev_before, :aarpenv_p_env_prev_after

  if Δ < 0 && -Δ > DATA_pers_and_parish_cur[2][p_position, t, 13]
    log_q_ratio = -Inf
    AddRem_penv_track[4] = 3

    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_penv_track)
  end


  ### Calculate the update ###

  DATA_res_and_track_prime, DATA_pers_and_parish_prime, scope, valid, AddRem_penv_track = update_data_change_init_penv(DATA_res_and_track_cur, DATA_pers_and_parish_cur, p_position, t, Δ, epi_params, f_to_p_structs, ids_to_pos_dict, AddRem_penv_track)

  if valid != 1
    log_q_ratio = -Inf
    AddRem_penv_track[4] = 2

    return(DATA_res_and_track_cur, DATA_pers_and_parish_cur, log_q_ratio, [0.,0.,0.,0.], AddRem_penv_track)
  end


  ### Calculate log q ratio ###

  log_q_ratio = 0.

  return(DATA_res_and_track_prime, DATA_pers_and_parish_prime, log_q_ratio, scope, AddRem_penv_track)
end
