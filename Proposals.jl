
include("Likelihood.jl")
include("DataUpdaters.jl")

###############
### Move SE ###
###############

function propose_Move_SE(combi_array_cur, f_to_p_dict)

  ### Choose a farm and a timestep ###

  farms_with_SE_events_at_t = combi_array_cur[2][(combi_array_cur[2][:, :, 13] .> 0), [1,3]]

  ### Generate the update parameters ###

  position, t = farms_with_SE_events_at_t[rand(1:size(farms_with_SE_events_at_t, 1)), :]

  Δ = 0
  while Δ == 0
    Δ = rand(Poisson(3))
    sgn_Δ = rand([-1, 1])
    Δ = sgn_Δ * Δ
  end

  if 360 < t + Δ || t + Δ <= 0
    log_q_ratio = -Inf
    Move_SE_track = [position, t, 0, 4, Δ, NaN]
    return(combi_array_cur, log_q_ratio, [0,0,0,0], Move_SE_track)
  end

  ### Calculate the update ###

  combi_array_prime, scope, valid = update_data_Move_SE(combi_array_cur, position, t, Δ, 1, f_to_p_dict)
  # num_SE_moved set to 1

  if valid != 1
    log_q_ratio = -Inf
    Move_SE_track = [position, t, 0, 2, Δ, NaN]
    return(combi_array_cur, log_q_ratio, [0,0,0,0], Move_SE_track)
  end

  Move_SE_track = [position, t, 0, 0, Δ, 1]
  # num_SE_moved set to 1

  return(combi_array_prime, 0, scope, Move_SE_track)
end


###############
### Move EI ###
###############

function propose_Move_EI(combi_array_cur, epi_params, f_to_p_dict)

  ### Choose a farm and a timestep ###

  farms_with_EI_events_at_t = combi_array_cur[2][(combi_array_cur[2][:, :, 14] .> 0), [1,3]]

  ### Generate the update parameters ###

  position, t = farms_with_EI_events_at_t[rand(1:size(farms_with_EI_events_at_t, 1)), :]

  Δ = 0
  while Δ == 0
    Δ = rand(Poisson(3))
    sgn_Δ = rand([-1, 1])
    Δ = sgn_Δ * Δ
  end

  if 360 < t + Δ || t + Δ <= 0
    log_q_ratio = -Inf
    Move_EI_track = [position, t, 0, 4, Δ, NaN]
    return(combi_array_cur, log_q_ratio, [0,0,0,0], Move_EI_track)
  end

  ### Calculate the update ###

  combi_array_prime, scope, valid = update_data_Move_EI(combi_array_cur, position, t, Δ, 1, epi_params, f_to_p_dict)
  # num_EI_moved set to 1

  if valid != 1
    log_q_ratio = -Inf
    Move_EI_track = [position, t, 0, 2, Δ, NaN]
    return(combi_array_cur, log_q_ratio, [0,0,0,0], Move_EI_track)
  end

  Move_EI_track = [position, t, 0, 0, Δ, 1]
  # num_EI_moved set to 1

  return(combi_array_prime, 0, scope, Move_EI_track)
end


#################
### AddRem SE ###
#################

function propose_AddRem_SE(combi_array_cur, f_to_p_dict)

  ### Choose a farm and a timestep ###

  farms_with_S_and_exp_prob_at_t = combi_array_cur[2][(combi_array_cur[1][:, :, 10] .> 0  .&& combi_array_cur[3][:, :, 4] .> 0), [1,3]]
  # Extract all farms that have cS_postM_t > 0 AND prob_exp_t > 0

  ### Generate the update parameters ###

  position, t = farms_with_S_and_exp_prob_at_t[rand(1:size(farms_with_S_and_exp_prob_at_t, 1)), :]

  ### Calculate the update ###

  combi_array_prime, scope, valid = update_data_AddRem_SE(combi_array_cur, position, t, 1, f_to_p_dict)
  # Δ = 1

  if valid != 1
    log_q_ratio = -Inf
    AddRem_SE_track = [position, t, 0, 2, 1, NaN]
    # Δ set to 1
    return(combi_array_cur, log_q_ratio, [0,0,0,0], AddRem_SE_track)
  end

  AddRem_SE_track = [position, t, 0, 0, 1, 1]
  # Δ set to 1

  return(combi_array_prime, 0, scope, AddRem_SE_track)
end


#################
### AddRem EI ###
#################

function propose_AddRem_EI(combi_array_cur, epi_params, f_to_p_dict)

  ### Choose a farm and a timestep ###

  farms_with_E_at_t = combi_array_cur[2][(combi_array_cur[1][:, :, 11] .> 0), [1,3]]
  # Extract all farms that have cE_postM_t > 0

  ### Generate the update parameters ###

  position, t = farms_with_E_at_t[rand(1:size(farms_with_E_at_t, 1)), :]

  ### Calculate the update ###

  combi_array_prime, scope, valid = update_data_AddRem_EI(combi_array_cur, position, t, 1, epi_params, f_to_p_dict)
  # Δ = 1

  if valid != 1
    log_q_ratio = -Inf
    AddRem_EI_track = [position, t, 0, 2, 1, NaN]
    # Δ set to 1
    return(combi_array_cur, log_q_ratio, [0,0,0,0], AddRem_EI_track)
  end

  AddRem_EI_track = [position, t, 0, 0, 1, 1]
  # Δ set to 1

  return(combi_array_prime, 0, scope, AddRem_EI_track)
end
