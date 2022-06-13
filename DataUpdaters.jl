
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

function update_pers_EPIDEMIC(combi_array_cur, epi_params_draw, scope)

  combi_array_prime = deepcopy(combi_array_cur)

  t_start = 1
  t_end = size(combi_array_cur[1], 2)
  pos_ids = 1:size(combi_array_cur[1], 1)


  @inbounds for pos in pos_ids
    @inbounds for t in t_start:t_end

      # Update Cattle Exposure Probability

      combi_array_prime[3][pos, t, 4] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 4:6],
                                                       p_env_prev = combi_array_prime[4][Int64(f_to_p_dict[pos][2])],
                                                       β = epi_params_draw[1],
                                                       F = epi_params_draw[4])

      # Update Badger Exposure Probability

      combi_array_prime[3][pos, t, 5] = calc_exp_prop(;States_init = combi_array_prime[1][pos, t, 22:24],
                                                       p_env_prev = combi_array_prime[4][Int64(f_to_p_dict[pos][2])],
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


####################################################
### Update the data after data augmentation draw ###
####################################################

function update_data_Move_SE(combi_array_cur, pos_id, t, Δ, num_SE_moved)

  combi_array_prime = deepcopy(combi_array_cur)

  ### Functional objects ###
  sgnΔ = sign(Δ)
  A = convert(Int, (1 - sgnΔ)/2)
  B = 1-A

  ### Scope ###

  lower_t = (t+(A*Δ))
  upper_t = (t+(B*Δ))
  h_pos_ids = [pos_id]
  h_element_range = 1:13

  ###########################
  ### Generate new states ###
  ###########################

  ##########
  ### Update the events
  ##########

  combi_array_prime[2][pos_id, t, 13] -= num_SE_moved
  combi_array_prime[2][pos_id, (t+Δ), 13] += num_SE_moved


  ############
  ### Update the states
  ############

  # :cS_init, :cS_Moves, :cS_postM
  combi_array_prime[1][pos_id, (lower_t+1):(upper_t), [4,7,10]] .+= sgnΔ * num_SE_moved
  # :cS_postEI, :cS_postDet, :cS_final
  combi_array_prime[1][pos_id, (lower_t):(upper_t-1), [13,16,19]] .+= sgnΔ * num_SE_moved

  # :cE_init, :cE_Moves, :cE_postM
  combi_array_prime[1][pos_id, (lower_t+1):(upper_t), [5,8,11]] .-= sgnΔ * num_SE_moved
  # :cE_postEI, :cE_postDet, :cE_final
  combi_array_prime[1][pos_id, (lower_t):(upper_t-1), [14,17,20]] .-= sgnΔ * num_SE_moved


  ##############
  ### Update the parish states
  ##############

  # :pcS_init
  combi_array_prime[4][f_to_p_dict[pos_id][2], (lower_t):(upper_t-1), 4] .+= sgnΔ * num_SE_moved
  # :pcS_final
  combi_array_prime[4][f_to_p_dict[pos_id][2], (lower_t):(upper_t-1), 7] .+= sgnΔ * num_SE_moved
  # :pcE_init
  combi_array_prime[4][f_to_p_dict[pos_id][2], (lower_t+1):(upper_t), 5] .-= sgnΔ * num_SE_moved
  # :pcE_final
  combi_array_prime[4][f_to_p_dict[pos_id][2], (lower_t):(upper_t-1), 8] .-= sgnΔ * num_SE_moved


  ###############
  ### Quick check for validity
  ###############

  posi_check = (combi_array_prime[1][pos_id, (lower_t):(upper_t), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0)

  if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
    return(combi_array, [lower_t, upper_t, h_pos_ids, h_element_range], 0)
                                                                      # invalid
  end

  return(combi_array_prime, [lower_t, upper_t, h_pos_ids, h_element_range], 1)
                                                                      # valid
end
