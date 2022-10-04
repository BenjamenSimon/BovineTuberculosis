
###################################################################
### Functions to calculate exposure and infection probabilities ###
###################################################################

function calc_exp_prop(;States_init, p_env_prev, β, F)

  I = States_init[3]
  N = sum(States_init)

  if N > 0
    exp_prob = 1 - exp(- (I/N) * β - F * p_env_prev)
  else
    exp_prob = 0
  end

  return(exp_prob)
end

function calc_inf_prob(γ)

  inf_prob = 1 - exp(-γ)

  return(inf_prob)
end


###################################################
### Update the persistents after parameter draw ###
###################################################

function update_pers_EPIDEMIC(combi_array_cur, log_epi_params_draw, f_to_p_dict, scope)

  combi_array_prime = deepcopy(combi_array_cur)
  epi_params_draw = exp.(log_epi_params_draw)

  t_start = 1
  t_end = size(combi_array_cur[1], 2)
  positions = 1:size(combi_array_cur[1], 1)


  @inbounds for pos in positions
    @inbounds for t in t_start:t_end

      # Update Cattle Exposure Probability

      combi_array_prime[3][pos, t, 4] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 4:6],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2], t, 19],
                                                       β = epi_params_draw[1],
                                                       F = epi_params_draw[4])

      # Update Badger Exposure Probability

      combi_array_prime[3][pos, t, 5] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 22:24],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2], t, 19],
                                                       β = epi_params_draw[2],
                                                       F = epi_params_draw[4])


      # Update Cattle Infection Probability

      combi_array_prime[3][pos, t, 6] = calc_inf_prob(epi_params_draw[3])

      # Update Badger Infection Probability

      combi_array_prime[3][pos, t, 7] = calc_inf_prob(epi_params_draw[3])

    end
  end

  return(combi_array_prime)
end

function update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, scope)

  t_start = scope[1]
  t_end = scope[2]
  positions = scope[3]


  @inbounds for pos in positions
    @inbounds for t in t_start:t_end

      # Update Cattle Exposure Probability

      combi_array_prime[3][pos, t, 4] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 4:6],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2], t, 19],
                                                       β = epi_params[1],
                                                       F = epi_params[4])


      # Update Cattle Infection Probability

      combi_array_prime[3][pos, t, 6] = calc_inf_prob(epi_params[3])

    end
  end

  return(combi_array_prime)
end


####################################################
### Update the data after data augmentation draw ###
####################################################

### MOVE SE ###

function update_data_Move_SE(combi_array_cur, position, t, Δ, num_SE_moved, f_to_p_dict)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Functional objects ###
  sgnΔ = sign(Δ)
  A = convert(Int, (1 - sgnΔ)/2)
  B = 1-A

  ### Scope ###

  lower_t = (t+(A*Δ))
  upper_t = (t+(B*Δ))
  h_positions = [position]
  h_element_range = 1:13

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, 13] -= num_SE_moved
  combi_array_prime[2][position, (t+Δ), 13] += num_SE_moved


  ############
  ### Update the states
  ############

  # :cS_init, :cS_Moves, :cS_postM
  combi_array_prime[1][position, (lower_t+1):(upper_t), [4,7,10]] .+= sgnΔ * num_SE_moved
  # :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, (lower_t):(upper_t-1), [13,16,19]] .+= sgnΔ * num_SE_moved

  # :cE_init, :cE_Moves, :cE_postM
  combi_array_prime[1][position, (lower_t+1):(upper_t), [5,8,11]] .-= sgnΔ * num_SE_moved
  # :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t):(upper_t-1), [14,17,20]] .-= sgnΔ * num_SE_moved


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcS_init
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 4] .+= sgnΔ * num_SE_moved
  # :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 7] .+= sgnΔ * num_SE_moved
  # :pcE_init
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):(upper_t), 5] .-= sgnΔ * num_SE_moved
  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 8] .-= sgnΔ * num_SE_moved


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end

### MOVE EI ###

function update_data_Move_EI(combi_array_cur, position, t, Δ, num_EI_moved, epi_params, f_to_p_dict)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Functional objects ###
  sgnΔ = sign(Δ)
  A = convert(Int, (1 - sgnΔ)/2)
  B = 1-A

  ### Scope ###

  lower_t = (t+(A*Δ))
  upper_t = (t+(B*Δ))
  h_positions = [position]
  h_element_range = 1:13

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, 14] -= num_EI_moved
  combi_array_prime[2][position, (t+Δ), 14] += num_EI_moved


  ############
  ### Update the states
  ############

  # :cE_init, :cE_Moves, :cE_postM
  combi_array_prime[1][position, (lower_t+1):(upper_t), [5,8,11]] .+= sgnΔ * num_EI_moved
  # :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t):(upper_t-1), [14,17,20]] .+= sgnΔ * num_EI_moved

  # :cI_init, :cI_Moves, :cI_postM
  combi_array_prime[1][position, (lower_t+1):(upper_t), [6,9,12]] .-= sgnΔ * num_EI_moved
  # :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, (lower_t):(upper_t-1), [15,18,21]] .-= sgnΔ * num_EI_moved


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), [5,6,8,9,11,12,14,15,17,18,20,21]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcE_init
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 5] .+= sgnΔ * num_EI_moved
  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 8] .+= sgnΔ * num_EI_moved
  # :pcI_init
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):(upper_t), 6] .-= sgnΔ * num_EI_moved
  # :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t):(upper_t-1), 9] .-= sgnΔ * num_EI_moved


  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])



  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end


### ADD/REM SE ###

function update_data_AddRem_SE(combi_array_cur, position, t, Δ, f_to_p_dict, tracker)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = size(combi_array_cur[1], 2) #T
  h_positions = [position]
  h_element_range = 1:13

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, 13] += Δ


  ############
  ### Update the states
  ############

  # :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, lower_t, [13,16,19]] .-= Δ
  # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [4,7,10,13,16,19]] .-= Δ

  # :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t+1):(upper_t), [5,8,11]] .+= Δ
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, lower_t, [5,8,11,14,17,20]] .+= Δ

  tracker[6:9] = [combi_array_cur[2][position, t, 13], combi_array_prime[2][position, t, 13],  combi_array_cur[1][position, t, 4], combi_array_cur[3][position, t, 4]]
               #  :arSE_SE_before, :arSE_SE_after, :arSE_cS, :arSE_prob

  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0, tracker)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 7] -= Δ
  # :pcS_init, :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [4,7]] .-= Δ

  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 8] += Δ
  # :pcE_init, :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [5,8]] .+= Δ


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1, tracker)
                                                                          # valid
end

### ADD/REM EI ###

function update_data_AddRem_EI(combi_array_cur, position, t, Δ, epi_params, f_to_p_dict, tracker)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = size(combi_array_cur[1], 2) #T
  h_positions = [position]
  h_element_range = 1:13

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, 14] += Δ


  ############
  ### Update the states
  ############

  # :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, lower_t, [14,17,20]] .-= Δ
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [5,8,11,14,17,20]] .-= Δ

  # :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, (lower_t+1):(upper_t), [6,9,12]] .+= Δ
  # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, lower_t, [6,9,12,15,18,21]] .+= Δ


  tracker[6:9] = [combi_array_cur[2][position, t, 14], combi_array_prime[2][position, t, 14],  combi_array_cur[1][position, t, 5], combi_array_cur[3][position, t, 5]]
                # :arEI_EI_before, :arEI_EI_after, :arEI_cE, :arEI_prob

  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), [5,6,8,9,11,12,14,15,17,18,20,21]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0, tracker)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 8] -= Δ
  # :pcE_init, :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [5,8]] .-= Δ

  # :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 9] += Δ
  # :pcI_init, :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [6,9]] .+= Δ


  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1, tracker)
                                                                      # valid
end


### ADD/REM Detections ###

function update_data_AddRem_Det(combi_array_cur, position, t, Δs, epi_params, f_to_p_dict, tracker)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = size(combi_array_cur[1], 2) #T
  h_positions = [position]
  h_element_range = 1:13

  # Δs = [ΔE, ΔI]
  # Check ΔE = 0 outside of func

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, [19,20]] += Δs


  ############
  ### Update the states
  ############


  # cE_postDet, :cE_final
  combi_array_prime[1][position, lower_t, [17, 20]] .-= Δs[1]
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [5,8,11,14,17,20]] .-= Δs[1]


  # :cI_postDet, :cI_final
  combi_array_prime[1][position, lower_t, [18, 21]] .-= Δs[2]
  # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [6,9,12,15,18,21]] .-= Δs[2]


  tracker[7:12] = [combi_array_cur[2][position, t, 19:20] ; combi_array_prime[2][position, t, 19:20] ;  combi_array_cur[1][position, t, 14:15]]
                # :arDet_Edet_before, :arDet_Idet_before, :arDet_Edet_after, :arDet_Idet_after, :arDet_cE, :arDet_cI


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), [5,6,8,9,11,12,14,15,17,18,20,21]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0, tracker)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 8] -= Δs[1]
  # :pcE_init, :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [5,8]] .-= Δs[1]

  # :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 9] -= Δs[2]
  # :pcI_init, :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [6,9]] .-= Δs[2]


  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1, tracker)
                                                                      # valid
end


### ADD/REM Deaths ###

function update_data_AddRem_Deaths(combi_array_cur, position, t, Δs, epi_params, f_to_p_dict, tracker)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = size(combi_array_cur[1], 2) #T
  h_positions = [position]
  h_element_range = 1:13

  # Δs = [ΔS, ΔE, ΔI]
  # Check sum(Δs) != 0 outside of func

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][position, t, [22,23,24]] += Δs


  ############
  ### Update the states
  ############

  # :cS_final
  combi_array_prime[1][position, lower_t, [19]] .-= Δs[1]
  # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [4,7,10,13,16,19]] .-= Δs[1]


  # :cE_final
  combi_array_prime[1][position, lower_t, [20]] .-= Δs[2]
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [5,8,11,14,17,20]] .-= Δs[2]


  # :cI_final
  combi_array_prime[1][position, lower_t, [21]] .-= Δs[3]
  # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, (lower_t+1):upper_t, [6,9,12,15,18,21]] .-= Δs[3]


  tracker[8:16] = [combi_array_cur[2][position, t, 22:24] ; combi_array_prime[2][position, t, 22:24] ;  combi_array_cur[1][position, t, 16:18]]
                # :arDeaths_Sdths_before, :arDeaths_Edths_before, :arDeaths_Idths_before,
                # :arDeaths_Sdths_after, :arDeaths_Edths_after, :arDeaths_Idths_after,
                # :arDeaths_cS, :arDeaths_cE, :arDeaths_cI


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), 4:21] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0, tracker)
                                                                      # invalid
  end


  ##############
  ### Update the parish states
  ##############

  # :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 7] -= Δs[1]
  # :pcS_init, :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [4,7]] .-= Δs[1]

  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 8] -= Δs[2]
  # :pcE_init, :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [5,8]] .-= Δs[2]

  # :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 9] -= Δs[3]
  # :pcI_init, :pcI_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [6,9]] .-= Δs[3]


  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1, tracker)
                                                                      # valid
end


### ADD/REM Movements ###

function update_data_AddRem_Movement(combi_array_cur, scope, combi_array_prime, epi_params, differences_oi, parish_differences)

  t, T, h_positions, position = scope[[1,2,3,5]]

  ##########################################
  ### Update the farm that moved animals ###
  ##########################################

  Δ_off = combi_array_prime[2][position, t, [7,8,9]] - combi_array_cur[2][position, t, [7,8,9]]

  # :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, t, [10,13,16,19]] .-= Δ_off[1]
  # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][position, (t+1):T, [4,7,10,13,16,19]] .-= Δ_off[1]


  # :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, t, [11,14,17,20]] .-= Δ_off[2]
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][position, (t+1):T, [5,8,11,14,17,20]] .-= Δ_off[2]


  # :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, t, [12,15,18,21]] .-= Δ_off[3]
  # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_array_prime[1][position, (t+1):T, [6,9,12,15,18,21]] .-= Δ_off[3]


  ##############################################
  ### Update the farms that recieved animals ###
  ##############################################

  for j in 1:size(differences_oi, 1)

      pos = differences_oi[j, 5]

      Δ_on_j = differences_oi[j, 1:3]

      # :cS_postM, :cS_postEI, :cS_postDet, :cS_final
      combi_array_prime[1][pos, t, [10,13,16,19]] .+= Δ_on_j[1]
      # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
      combi_array_prime[1][pos, (t+1):T, [4,7,10,13,16,19]] .+= Δ_on_j[1]


      # :cE_postM, :cE_postEI, :cE_postDet, :cE_final
      combi_array_prime[1][pos, t, [11,14,17,20]] .+= Δ_on_j[2]
      # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
      combi_array_prime[1][pos, (t+1):T, [5,8,11,14,17,20]] .+= Δ_on_j[2]


      # :cI_postM, :cI_postEI, :cI_postDet, :cI_final
      combi_array_prime[1][pos, t, [12,15,18,21]] .+= Δ_on_j[3]
      # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
      combi_array_prime[1][pos, (t+1):T, [6,9,12,15,18,21]] .+= Δ_on_j[3]


      # :sus_on, :exp_on, :inf_on
      combi_array_prime[2][pos, t, [4,5,6]] .+= Δ_on_j

      ##############
      ### Update the probabilities
      ##############

      combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [t, T, pos])

  end # end of for each move


  ################################
  ### Update the parish totals ###
  ################################

  for j in 1:size(parish_differences, 1)

    p_pos = parish_differences[j, 4]

    Δ_pj = parish_differences[j, 1:3]


    if p_pos > 0

      # :pcS_final
      combi_array_prime[4][p_pos, t, 7] += Δ_pj[1]
      # :pcS_init, :pcS_final
      combi_array_prime[4][p_pos, (t+1):T, [4,7]] .+= Δ_pj[1]

      # :pcE_final
      combi_array_prime[4][p_pos, t, 8] += Δ_pj[2]
      # :pcE_init, :pcE_final
      combi_array_prime[4][p_pos, (t+1):T, [5,8]] .+= Δ_pj[2]

      # :pcI_final
      combi_array_prime[4][p_pos, t, 9] += Δ_pj[3]
      # :pcI_init, :pcI_final
      combi_array_prime[4][p_pos, (t+1):T, [6,9]] .+= Δ_pj[3]

    end # end if p_pos > 0

  end # end of for each parish


  ####################################
  ### Early return: Invalid Update ###
  ####################################

  for chpos in h_positions
    posi_check = combi_array_prime[1][h_positions, t:T, 4:21] .>= 0

    if sum(posi_check) != prod(size(posi_check))
      # returns changed position as a Int instead of a Vector
      return(combi_array_prime, 2)
    end
  end

  return(combi_array_prime, 1)
                          # valid
end


### ADD/REM Environmental Pressure ###

function update_data_AddRem_penv(combi_array_cur, p_position, t, Δs, epi_params, f_to_p_dict, ids_to_pos_dict, tracker)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = t+1
  h_positions = f_to_p_dict[ids_to_pos_dict[combi_array_prime[4][p_position, t, 23]]][4]
  h_element_range = [3,8]

  # Δs = [Δr, Δn]
  # Check sum(Δs) != 0 outside of func

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  # :remaining_pressure, :new_pressure
  combi_array_prime[4][p_position, t, [16, 17]] += Δs

  Δpenv = sum(Δs)/combi_array_prime[4][p_position, t, 18] # /scaling

  ############
  ### Update the states
  ############

  # :p_env_prev
  combi_array_prime[4][p_position, (t+1), 19] += Δpenv
  # :p_env_cur
  combi_array_prime[4][p_position, t, 20] += Δpenv


  tracker[7:12] = [combi_array_cur[4][p_position, t, 16:17] ; combi_array_prime[4][p_position, t, 16:17] ;  combi_array_cur[4][p_position, (t+1), [19,6]]]
                # :arpenv_r_pres_before, :arpenv_n_pres_before, :arpenv_r_pres_after, :arpenv_n_pres_after, :arpenv_p_env_prev, :arpenv_pI

  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_pers_EPIDEMIC(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1, tracker)
                                                                      # valid
end
