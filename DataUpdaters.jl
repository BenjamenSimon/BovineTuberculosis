
include("Likelihood.jl")

###################################################################
### Functions to calculate exposure and infection probabilities ###
###################################################################

function calc_exp_prop(;States_init, p_env_prev, β, F)

  I = States_init[3]
  N = sum(States_init)

  if N > 0
    exp_prob = 1 - exp(- I/N * β - F * p_env_prev)
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

function update_pers_EPIDEMIC(combi_array_cur, epi_params_draw, f_to_p_dict, scope)

  combi_array_prime = deepcopy(combi_array_cur)

  t_start = 1
  t_end = size(combi_array_cur[1], 2)
  positions = 1:size(combi_array_cur[1], 1)


  @inbounds for pos in positions
    @inbounds for t in t_start:t_end

      # Update Cattle Exposure Probability

      combi_array_prime[3][pos, t, 4] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 4:6],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2]],
                                                       β = epi_params_draw[1],
                                                       F = epi_params_draw[4])

      # Update Badger Exposure Probability

      combi_array_prime[3][pos, t, 5] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 22:24],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2]],
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

# combi_array[4][f_to_p_dict[combi_array[1][1,1,2]][2]]
#
# combi_array[4][Int64(f_to_p_dict[1][2])]

function update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, scope)

  t_start = scope[1]
  t_end = scope[2]
  positions = scope[3]


  @inbounds for pos in positions
    @inbounds for t in t_start:t_end

      # Update Cattle Exposure Probability

      combi_array_prime[3][pos, t, 4] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 4:6],
                                                       p_env_prev = combi_array_prime[4][f_to_p_dict[pos][2]],
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

function update_data_AddRem_SE(combi_array_cur, position, t, Δ, f_to_p_dict)

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

  # :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 7] -= Δ
  # :pcS_init, :pcS_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [4,7]] .-= Δ

  # :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], lower_t, 8] += Δ
  # :pcE_init, :pcE_final
  combi_array_prime[4][f_to_p_dict[position][2], (lower_t+1):upper_t, [5,8]] .+= Δ


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end

### ADD/REM EI ###

function update_data_AddRem_EI(combi_array_cur, position, t, Δ, epi_params, f_to_p_dict)

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


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end


### ADD/REM Detections ###

function detection_permutations(cStates_postEI, cDet_cur, epi_params)

  total_det_t = sum(cDet_cur)

  permutations = fill(-99., Int64(factorial(total_det_t+2-1)/(factorial(total_det_t)*factorial(2-1))), 3)

  iter = 1

  for nE in 0:total_det_t

    nI = (total_det_t - nE)

    prob_Edet = logpdf(Binomial(convert(Int64, cStates_postEI[2]), prod(epi_params[[6,7]])), nE)
    prob_Idet = logpdf(Binomial(convert(Int64, cStates_postEI[3]), epi_params[6]), nI)

    permutations[iter, :] = [nE, nI, (exp(prob_Edet) + exp(prob_Idet))]
    iter += 1
  end

  println(permutations)

  chosen_perm = [permutations[wsample(1:size(permutations, 1), permutations[:, 3]), :] ; permutations[convert(Int64, (cDet_cur[1]+1)), 3]]

  return(chosen_perm)
end #to test


function update_data_AddRem_Det(combi_array_cur, position, t, Δs, epi_params, f_to_p_dict)

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


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end


### ADD/REM Deaths ###

function deaths_permutations(cStates_postDet, cDeaths, epi_params)

  total_deaths = sum(cDeaths)

  permutations = fill(-99., Int64(factorial(total_deaths+3-1)/(factorial(total_deaths)*factorial(3-1))), 4)

  iter = 1

  for nS in 0:total_deaths
    for nE in 0:(total_deaths - nS)

      nI = (total_deaths - nS - nE)

      perm_prob = log_pdf_mvhyper(cStates_postDet, [nS, nE, nI])

      permutations[iter, :] = [nS, nE, nI, exp(perm_prob)]
      iter += 1
    end
  end

  chosen_perm = [permutations[wsample(1:size(permutations, 1), permutations[:, 4]), :] ; log_pdf_mvhyper(cStates_postDet, cDeaths)]

  return(chosen_perm)
end #to test


function update_data_AddRem_Deaths(combi_array_cur, position, t, Δs, epi_params, f_to_p_dict)

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


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][position, (lower_t):(upper_t), 4:21] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_positions, h_element_range], 0)
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


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end


### ADD/REM Movements ###

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

function generate_new_movement(combi_array_cur, position, t, epi_params, movement_record, dict_of_movements, f_to_p_dict, ids_to_pos_dict)

  combi_array_prime = deepcopy(combi_array_cur)
  movement_record_prime = deepcopy(movement_record)

  ### Scope ###

  lower_t = t
  upper_t = size(combi_array_cur[1], 2) #T
  h_positions = [position] # extended dependent on movement choice
  h_element_range = 1:13


  ####################
  ### Extract Data ###
  ####################

  # :cS_Moves, :cE_Moves, :cI_Moves
  states = combi_array_prime[1][position, t, [7,8,9]]
  # :sus_off, :exp_off, :inf_off
  moves_off = combi_array_prime[2][position, t, [7,8,9]]

  total_moves = sum(moves_off)


  ############################################
  ### Extract the movement data for time t ###
  ############################################

  move_record_rows = dict_of_movements[(position, t)]

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

  new_move_states = rng_mvhyper(states, total_moves)

    ###################################
    ### Early return: No difference ###
    ###################################

    if all((new_move_states - moves_off) .== 0)
      # returns changed position as an Int instead of a Vector
      return(combi_array, movement_record, 2)
    end

  running_move_states = deepcopy(new_move_states)

  # FOR each (j) of the farms that recieved animals
  for j in 1:size(on_cph_row_ids, 1)

    # IF there is more than one farm that recieved animals from farm i
    if (size(on_cph_row_ids, 1) > 1)
      # Extract the total movements onto this farm
      total_moved = m_farm[j]
      # Generate the states
      states_on = rng_mvhyper(running_move_states, total_moved)
    else # ELSE there is only one farm
      states_on = running_move_states
    end

    moves_data_prime[j, 4:6] = running_move_states
    moves_data_prime[j, 7:9] = states_on

    # Update the running total of the moved states
    # Reduce the moved states by those that just got assigned
    running_move_states = running_move_states - states_on

  end #end of for farms that recieved animals


  ########################
  ### Update moves off ###
  ########################

  # :sus_off, :exp_off, :inf_off
  combi_array_prime[2][position, t, [7,8,9]] = new_move_states

  movement_record_prime[move_record_rows, :] .= moves_data_prime


  ############################################
  ### Calculate the farm level differences ###
  ############################################

  # :Δ_S_on, :Δ_E_on, :Δ_I_on, :on_row_id, :on_position, :on_parish_pos
  differences = fill(-99., size(on_cph_row_ids, 1), 6)

  for j in 1:size(on_cph_row_ids, 1)

      if on_cph_row_ids[j] > 0
        differences[j, 1:6] = [Array(moves_data_cur[j, 7:9]) - Array(moves_data_prime[j, 7:9]) ;
                              on_cph_row_ids[j] ;
                              ids_to_pos_dict[on_cph_row_ids[j]] ;
                              f_to_p_dict[ids_to_pos_dict[on_cph_row_ids[j]]][2] ]


        h_positions = [h_positions ; differences[j, 5]]
      end
  end

  h_positions_oi = unique(h_positions[h_positions .> 0])

  differences_oi = differences[(differences[:,4] .> 0), :]

  ##############################################
  ### Calculate the parish level differences ###
  ##############################################

  affected_parishes = unique([differences_oi[:, 6]; f_to_p_dict[position][2]])

  # :Δ_pS, :Δ_pE, :Δ_pI, :parish_pos
  parish_differences = fill(0, size(affected_parishes, 1), 4)

  for (idx, par) in enumerate(affected_parishes)

    parish_differences[idx, 4] = par

    for i in 1:size(differences_oi, 1)
      if differences_oi[i, 5] == par
        parish_differences[idx, 1:3] += differences_oi[i, 1:3]
      end
    end

    if par == f_to_p_dict[position][2] #off_p_id
      Δ_off = combi_array_prime[2][position, t, [7,8,9]] - combi_array_cur[2][position, t, [7,8,9]]
      parish_differences[idx, 1:3] -= Δ_off
    end

  end


  ############################################
  ### Update the data with the differences ###
  ############################################

  combi_array_prime = update_data_after_move(combi_array_cur, h_positions, t, upper_t, combi_array_prime, epi_params, differences_oi, parish_differences)


  ####################################
  ### Early return: Invalid Update ###
  ####################################

  for chpos in h_positions
    posi_check = combi_array_prime[1][h_positions, lower_t:upper_t, 4:21] .>= 0

    if sum(posi_check) != prod(size(posi_check))
      # returns changed position as a Int instead of a Vector
      return(combi_array, movement_record, 3)
    end
  end

  return(combi_array_prime, movement_record_prime, 1)
                                                  #valid
end


function update_data_after_move(combi_array_cur, position, t, T, combi_array_prime, epi_params, differences_oi, parish_differences)

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

      combi_array_prime = update_cattle_pers_general(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, pos])

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

  return(combi_array_prime)
end


### ADD/REM Environmental Pressure ###

function update_data_AddRem_penv(combi_array_cur, p_position, t, Δs, epi_params, f_to_p_dict, ids_to_pos_dict)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Scope ###

  lower_t = t
  upper_t = t+1
  h_positions = f_to_p_dict[ids_to_pos_dict[combi_array_prime[4][p_position, t, 23]]][4]
  h_element_range = 1:13

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
  combi_array_prime[4][p_position, t, 19] += Δpenv
  # :p_env_cur
  combi_array_prime[4][p_position, t, 20] += Δpenv


  ##############
  ### Update the probabilities
  ##############

  combi_array_prime = update_pers_EPIDEMIC(combi_array_prime, epi_params, f_to_p_dict, [lower_t, upper_t, h_positions])


  return(combi_array_prime, [lower_t, upper_t, h_positions, h_element_range], 1)
                                                                      # valid
end
