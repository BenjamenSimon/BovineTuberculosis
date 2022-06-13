using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames
using DataFramesMeta
using LinearAlgebra
using RData

##################################
##### Moving an S -> E event #####
##################################

function propose_Move_SE_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    choosing_a_farm = fill(0, 500, 2)
    tr_nm = 1

    for timestep in 1:360
      for position_id in 1:100
        if combi_dfs[position_id, 2][timestep , 12] > 0
          choosing_a_farm[tr_nm, :] = [position_id, timestep]
          tr_nm += 1
        end
      end
    end

    choosing_a_farm = choosing_a_farm[findall(choosing_a_farm[:,2] .> 0), :]

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time, 1]
    t = choosing_a_farm[choose_farm_time, 2]
    Δ = 0
    while Δ == 0
      Δ = rand(Poisson(3))
      sgn_Δ = rand([-1, 1])
      Δ = sgn_Δ * Δ
    end

    if 360 < t + Δ || t + Δ <= 0
      log_q_ratio = -Inf
      Move_SE_track = [pos_id, t, 0, 4, Δ, NaN]
      return(combi_dfs, log_q_ratio, [0,0,0,0], Move_SE_track)
    end

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of S to E events at time t
    num_SE_t =  track[t, 12]

    # Generate a number to be moved
    # prob = rand()
    # num_SE_move_t = rand(Binomial(convert(Int64, num_SE_t), prob))
    num_SE_move_t = num_SE_t

  #### Generate new states####

    # Update the events
    track_prime[t , 12] -= num_SE_move_t
    track_prime[(t+Δ) , 12] += num_SE_move_t

    sgnΔ = sign(Δ)
    A = convert(Int, (1 - sgnΔ)/2)
    B = 1-A

    # Update the states

    # :cS_init, :cS_Moves, :cS_postM, :pcS_init
    res_prime[( (t+1+(A*Δ)):(t+(B*Δ)) ) , [3, 6, 9, 32]] .+= sgnΔ * num_SE_move_t
    # :cS_postEI, :cS_postDet, :cS_final, :pcS_final
    res_prime[( (t+(A*Δ)):(t-1+(B*Δ)) ) , [12, 15, 18, 35]] .+= sgnΔ * num_SE_move_t

    # :cE_init, :cE_Moves, :cE_postM, :pcE_init
    res_prime[( (t+1+(A*Δ)):(t+(B*Δ)) ) , [4, 7, 10, 33]] .-= sgnΔ * num_SE_move_t
    # :cE_postEI, :cE_postDet, :cE_final, :pcE_final
    res_prime[( (t+(A*Δ)):(t-1+(B*Δ)) ) , [13, 16, 19, 36]] .-= sgnΔ * num_SE_move_t

    #### Quick check for validity ####

    posi_check = res_prime[(t+(A*Δ)):(t+(B*Δ)), [3,4,6,7,9,10,12,13,15,16,18,19]] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      Move_SE_track = [pos_id, t, 0, 3, Δ, num_SE_move_t]
      return(combi_dfs, log_q_ratio, [0,0,0,0], Move_SE_track)
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcS_init = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 32]
    pcE_init = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 33]

    pcS_final = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 35]
    pcE_final = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 36]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in each state in the parish
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 32] = pcS_init
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 33] = pcE_init

      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 35] = pcS_final
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 36] = pcE_final
    end

    # Calculate the log q ratio (proposal density ratio)
    log_q_ratio = 0

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = (t+(A*Δ))
    upper_t = (t+(B*Δ))
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    Move_SE_track = [pos_id, t, 0, 0, Δ, num_SE_move_t]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], Move_SE_track)
end

# propose_Move_SE_update(combi_dfs=combi_dfs, epi_params = epi_params_true)

##################################
##### Moving an E -> I event #####
##################################

function propose_Move_EI_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        if combi_dfs[position_id, 2][timestep , 13] > 0
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    Δ = 0
    while Δ == 0
      Δ = rand(Poisson(3))
      sgn_Δ = rand([-1, 1])
      Δ = sgn_Δ * Δ
    end

    # Early return: Invalid update parameters
    if 360 < t + Δ || t + Δ <= 0
      log_q_ratio = -Inf
      Move_EI_track = [pos_id, t, 0, 4, Δ, NaN]
      return(combi_dfs, log_q_ratio, [0,0,0,0], Move_EI_track)
    end

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of S to E events at time t
    num_SE_t =  track[t, 13]

    # Generate a number to be moved
    # prob = rand()
    # num_SE_move_t = rand(Binomial(convert(Int64, num_SE_t), prob))
    num_SE_move_t = num_SE_t

  #### Generate new states ####

    # Update the events
    track_prime[t , 13] -= num_SE_move_t
    track_prime[(t+Δ) , 13] += num_SE_move_t

    sgnΔ = sign(Δ)
    A = convert(Int, (1 - sgnΔ)/2)
    B = 1-A

    # Update the states

    # :cE_init, :cE_Moves, :cE_postM, :pcE_init
    res_prime[( (t+1+(A*Δ)):(t+(B*Δ)) ) , [4, 7, 10, 33]] .+= sgnΔ * num_SE_move_t
    # :cE_postEI, :cE_postDet, :cE_final, :pcE_final
    res_prime[( (t+(A*Δ)):(t-1+(B*Δ)) ) , [13, 16, 19, 36]] .+= sgnΔ * num_SE_move_t

    # :cI_init, :cI_Moves, :cI_postM, :pcI_init
    res_prime[( (t+1+(A*Δ)):(t+(B*Δ)) ) , [5, 8, 11, 34]] .-= sgnΔ * num_SE_move_t
    # :cI_postEI, :cI_postDet, :cI_final, :pcI_final
    res_prime[( (t+(A*Δ)):(t-1+(B*Δ)) ) , [14, 17, 20, 37]] .-= sgnΔ * num_SE_move_t

    #### Quick check for validity ####

    posi_check = res_prime[(t+(A*Δ)):(t+(B*Δ)), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      Move_EI_track = [pos_id, t, 0, 3, Δ, num_SE_move_t]
      return(combi_dfs, log_q_ratio, [0,0,0,0], Move_EI_track)
    end

    # Calculate changes to probabilities

    β_c = epi_params[1]
    F = epi_params[4]

    for τ in (t+(A*Δ)):(t+(B*Δ))

      cI_init = res_prime[τ , 5]
      cN_init = sum(res_prime[τ , 3:5])

      p_env_prev = pers_prime[τ , 3]

      c_exp_prob = 0.

      if cN_init > 0.
        c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
      end

      # Alternative to loop
      # 1 .- exp.( -(cI_init ./ cN_init) * β_c - F * p_env_prev)
      # gives ans of 1 when cN_init = 0

      pers_prime[τ , 5] = c_exp_prob
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcE_init = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 33]
    pcI_init = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 34]

    pcE_final = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 36]
    pcI_final = res_prime[(t+(A*Δ)):(t+(B*Δ)) , 37]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 33] = pcE_init
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 34] = pcI_init

      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 36] = pcE_final
      combi_dfs_prime[farm,1][(t+(A*Δ)):(t+(B*Δ)), 37] = pcI_final
    end

    # Calculate the log q ratio (proposal density ratio)
    log_q_ratio = 0

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = (t+(A*Δ))
    upper_t = (t+(B*Δ))
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    Move_EI_track = [pos_id, t, 0, 0, Δ, num_SE_move_t]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], Move_EI_track)
end

# propose_Move_EI_update(combi_dfs=combi_dfs, epi_params = epi_params_true)



###########################################
##### Adding/Removing an S -> E event #####
###########################################

function propose_AddRem_SE_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    # Extract all farms that have cS_postM_t > 0
    # AND prob_exp_t > 0 for all t

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        if combi_dfs[position_id, 1][timestep , 9] > 0 && combi_dfs[position_id, 3][timestep , 5] > 0
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = 360

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of S to E events at time t
    num_SE_t =  track[t, 12]

    # Generate a new number of events
    cS_postM_t = res[t, 9]
    prob = pers[t, 5]
    new_SE_t = rand(Binomial(convert(Int64, cS_postM_t), prob))

    Δ = new_SE_t - num_SE_t

    if Δ == 0
      log_q_ratio = -Inf
      # println("No diff")
      AddRem_SE_track = [pos_id, t, 0, 2, num_SE_t, new_SE_t, Δ, cS_postM_t, prob]
      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_SE_track)
    end

  #### Generate new states####

    # Update the events
    track_prime[t , 12] += Δ

    # Update the states

    # :cS_postEI, :cS_postDet, :cS_final, :pcS_final
    res_prime[ t:t , [12, 15, 18, 35]] .-= Δ
    # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final, :pcS_init, :pcS_final
    res_prime[((t+1):T) , [3, 6, 9, 12, 15, 18, 32, 35]] .-= Δ

    # :cE_postEI, :cE_postDet, :cE_final, :pcE_final
    res_prime[ t:t , [13, 16, 19, 36]] .+= Δ
    # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final, :pcE_init, :pcE_final
    res_prime[((t+1):T) , [4, 7, 10, 13, 16, 19, 33, 36]] .+= Δ

    #### Quick check for validity ####

    posi_check = res_prime[(t:T), [3,4,6,7,9,10,12,13,15,16,18,19]] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      # println("Invalid")
      AddRem_SE_track = [pos_id, t, 0, 3, num_SE_t, new_SE_t, Δ, cS_postM_t, prob]
      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_SE_track)
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcS_init = res_prime[(t:T) , 32]
    pcE_init = res_prime[(t:T) , 33]

    pcS_final = res_prime[(t:T) , 35]
    pcE_final = res_prime[(t:T) , 36]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in each state in the parish
      combi_dfs_prime[farm,1][(t:T), 32] = pcS_init
      combi_dfs_prime[farm,1][(t:T), 33] = pcE_init

      combi_dfs_prime[farm,1][(t:T), 35] = pcS_final
      combi_dfs_prime[farm,1][(t:T), 36] = pcE_final
    end

    # Calculate the log q ratio (proposal density ratio)
    num_SE_t_prime =  track_prime[t, 12]
    cS_postM_t_prime = res_prime[t, 9]
    prob_prime = pers_prime[t, 5]

    q_prime_given_cur = logpdf(Binomial(convert(Int64, cS_postM_t), prob), num_SE_t_prime)
    q_cur_given_prime = logpdf(Binomial(convert(Int64, cS_postM_t_prime), prob_prime), num_SE_t)

    log_q_ratio = q_cur_given_prime - q_prime_given_cur

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = t
    upper_t = T
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    AddRem_SE_track = [pos_id, t, 0, 0, num_SE_t, num_SE_t_prime, Δ, cS_postM_t, prob]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], AddRem_SE_track)
end

# propose_AddRem_SE_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#
# for i in 1:100
#   a,b,c = propose_AddRem_SE_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: Very unliklely to propose an update, 10 in 1000 maybe, and then 1 of those will be valid



###########################################
##### Adding/Removing an E -> I event #####
###########################################

function propose_AddRem_EI_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    # Extract all farms that have cE_postM_t > 0 for all t

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        if combi_dfs[position_id, 1][timestep , 10] > 0
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = 360

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of E to I events at time t
    num_EI_t =  track[t, 13]

    # Generate a new number of events
    cE_postM_t = res[t, 10]
    prob = pers[t, 7]
    new_EI_t = rand(Binomial(convert(Int64, cE_postM_t), prob))

    Δ = new_EI_t - num_EI_t

    if Δ == 0
      log_q_ratio = -Inf
      # println("No diff")
      AddRem_EI_track = [pos_id, t, 0, 2, num_EI_t, new_EI_t, Δ, cE_postM_t, prob]

      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_EI_track)
    end

  #### Generate new states####

    # Update the events
    track_prime[t , 12] += Δ

    # Update the states

    # :cE_postEI, :cE_postDet, :cE_final, :pcE_final
    res_prime[ t:t , [13, 16, 19, 36]] .-= Δ
    # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final, :pcE_init, :pcE_final
    res_prime[((t+1):T) , [4, 7, 10, 13, 16, 19, 33, 36]] .-= Δ

    # :cI_postEI, :cI_postDet, :cI_final, :pcI_final
    res_prime[ t:t , [14, 17, 20, 37]] .+= Δ
    # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final, :pcI_init, :pcI_final
    res_prime[((t+1):T) , [5, 8, 11, 14, 17, 20, 34, 37]] .+= Δ

    #### Quick check for validity ####

    posi_check = res_prime[(t:T), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      # println("Invalid")
      AddRem_EI_track = [pos_id, t, 0, 3, num_EI_t, new_EI_t, Δ, cE_postM_t, prob]

      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_EI_track)
    end

    # Calculate changes to probabilities

    β_c = epi_params[1]
    F = epi_params[4]

    for τ in (t:T)

      cI_init = res_prime[τ , 5]
      cN_init = sum(res_prime[τ , 3:5])

      p_env_prev = pers_prime[τ , 3]

      c_exp_prob = 0.

      if cN_init > 0.
        c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
      end

      # Alternative to loop
      # 1 .- exp.( -(cI_init ./ cN_init) * β_c - F * p_env_prev)
      # gives ans of 1 when cN_init = 0

      pers_prime[τ , 5] = c_exp_prob
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcE_init = res_prime[(t:T) , 33]
    pcI_init = res_prime[(t:T) , 34]

    pcE_final = res_prime[(t:T) , 36]
    pcI_final = res_prime[(t:T) , 37]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      combi_dfs_prime[farm,1][(t:T), 33] = pcE_init
      combi_dfs_prime[farm,1][(t:T), 34] = pcI_init

      combi_dfs_prime[farm,1][(t:T), 36] = pcE_final
      combi_dfs_prime[farm,1][(t:T), 37] = pcI_final
    end

    # Calculate the log q ratio (proposal density ratio)
    num_EI_t_prime =  track_prime[t, 13]
    cE_postM_t_prime = res_prime[t, 10]
    prob_prime = pers_prime[t, 7]

    q_prime_given_cur = logpdf(Binomial(convert(Int64, cE_postM_t), prob), num_EI_t_prime)
    q_cur_given_prime = logpdf(Binomial(convert(Int64, cE_postM_t_prime), prob_prime), num_EI_t)

    log_q_ratio = q_cur_given_prime - q_prime_given_cur

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = t
    upper_t = T
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    AddRem_EI_track = [pos_id, t, 0, 0, num_EI_t, new_EI_t, Δ, cE_postM_t, prob]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], AddRem_EI_track)
end

# propose_AddRem_EI_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#
# for i in 1:1000
#   a,b,c = propose_AddRem_EI_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: Very unliklely to propose an update, 50 in 1000 maybe, and then 2 of those will be valid



##############################################
##### Adding/Removing an Detection event #####
##############################################

function detection_permutations(cE_postEI, cI_postEI, Edet_cur, Idet_cur, epi_params)

  inf_det_prob = epi_params[6]
  exp_det_prob = inf_det_prob * epi_params[7]

  total_det_t = Edet_cur + Idet_cur

  permutations = fill(-10., 5000, 3)
  tra_nm = 1

  for nE in 0:total_det_t

    nI = (total_det_t - nE)

    prob_Edet = logpdf(Binomial(convert(Int64, cE_postEI), exp_det_prob), nE)
    prob_Idet = logpdf(Binomial(convert(Int64, cI_postEI), inf_det_prob), nI)

    perm_prob = prob_Edet + prob_Idet

    permutations[tra_nm, :] = [nE, nI, perm_prob]
    tra_nm += 1
  end

  permutations = permutations[findall(permutations[:,1] .> -1), :]

  choose_perm = wsample(1:size(permutations, 1), permutations[:, 3])

  chosen_perm = [permutations[choose_perm, :] ; permutations[convert(Int64, (Edet_cur+1)), 3]]

  return(chosen_perm)
end

function propose_AddRem_Det_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    # Extract all farms that have detections for all t
    # Could also condition on test occur
    # but we are only permuting the detections
    # so if no detections happened theres no changes possible

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        if sum(combi_dfs[position_id, 2][timestep , 18:19]) > 0.
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = 360

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of E and I detection events at time t
    num_Edet_t =  track[t, 18]
    num_Idet_t =  track[t, 19]

    # Generate a new number of events
    cE_postEI_t = res[t, 10]
    cI_postEI_t = res[t, 10]
    chosen_perm = detection_permutations(cE_postEI_t, cI_postEI_t, num_Edet_t, num_Idet_t, epi_params)

    new_Edet_t = chosen_perm[1]
    new_Idet_t = chosen_perm[2]

    ΔE = new_Edet_t - num_Edet_t
    ΔI = new_Idet_t - num_Idet_t

    if ΔE == 0
      log_q_ratio = -Inf
      AddRem_dets_track = [pos_id, t, 0, 2, num_Edet_t, num_Idet_t, new_Edet_t, new_Idet_t, ΔE, ΔI, cE_postEI_t, cI_postEI_t]
      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_dets_track)
    end

  #### Generate new states####

    # Update the events
    track_prime[t , 18] += ΔE
    track_prime[t , 19] += ΔI

    # Update the states

    # :cE_postDet, :cE_final, :pcE_final
    res_prime[ t:t , [16, 19, 36]] .-= ΔE
    # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final, :pcE_init, :pcE_final
    res_prime[((t+1):T) , [4, 7, 10, 13, 16, 19, 33, 36]] .-= ΔE

    # :cI_postDet, :cI_final, :pcI_final
    res_prime[ t:t , [17, 20, 37]] .-= ΔI
    # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final, :pcI_init, :pcI_final
    res_prime[((t+1):T) , [5, 8, 11, 14, 17, 20, 34, 37]] .-= ΔI

    #### Quick check for validity ####

    posi_check = res_prime[(t:T), [4,5,7,8,10,11,13,14,16,17,19,20]] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      AddRem_dets_track = [pos_id, t, 0, 3, num_Edet_t, num_Idet_t, new_Edet_t, new_Idet_t, ΔE, ΔI, cE_postEI_t, cI_postEI_t]
      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_dets_track)
    end

    # Calculate changes to probabilities

    β_c = epi_params[1]
    F = epi_params[4]

    for τ in (t:T)

      cI_init = res_prime[τ , 5]
      cN_init = sum(res_prime[τ , 3:5])

      p_env_prev = pers_prime[τ , 3]

      c_exp_prob = 0.

      if cN_init > 0.
        c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
      end

      # Alternative to loop
      # 1 .- exp.( -(cI_init ./ cN_init) * β_c - F * p_env_prev)
      # gives ans of 1 when cN_init = 0

      pers_prime[τ , 5] = c_exp_prob
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcE_init = res_prime[(t:T) , 33]
    pcI_init = res_prime[(t:T) , 34]

    pcE_final = res_prime[(t:T) , 36]
    pcI_final = res_prime[(t:T) , 37]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      combi_dfs_prime[farm,1][(t:T), 33] = pcE_init
      combi_dfs_prime[farm,1][(t:T), 34] = pcI_init

      combi_dfs_prime[farm,1][(t:T), 36] = pcE_final
      combi_dfs_prime[farm,1][(t:T), 37] = pcI_final
    end

    # Calculate the log q ratio (proposal density ratio)
    q_prime_given_cur = chosen_perm[3]
    q_cur_given_prime = chosen_perm[4]

    log_q_ratio = q_cur_given_prime - q_prime_given_cur

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = t
    upper_t = T
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    AddRem_dets_track = [pos_id, t, 0, 0, num_Edet_t, num_Idet_t, new_Edet_t, new_Idet_t, ΔE, ΔI, cE_postEI_t, cI_postEI_t]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], AddRem_dets_track)
end

# propose_AddRem_Det_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#
# for i in 1:1000
#   a,b,c = propose_AddRem_Det_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: About 10 in 1000 valid

##############################################
##### Adding/Removing an Detection event #####
##############################################

function log_pdf_mvhyper(K, k)
  # K is the number of objects of each type i
  # k is a vector of the number of succeses of each type i
  # N is the total number of objects
  # n is the number of draws total

  N = sum(K)
  n = sum(k)

  m = length(K) # number of groups

  x = fill(0., m) # results

  for i in 1:m
    Kck_i = binomial(BigInt(K[i]), BigInt(k[i]))

    x[i] = log(Kck_i)
  end

  top = sum(x)
  bot = binomial(BigInt(N), BigInt(n))

  prob = top-log(bot)

  return(prob)
end

function deaths_permutations(cStates_postDet, cDeaths, epi_params)

  total_deaths = sum(cDeaths)

  permutations = fill(-10., 5000, 4)
  trac_nm = 1

  for nS in 0:total_deaths
    for nE in 0:(total_deaths - nS)

      nI = (total_deaths - nS - nE)

      perm_prob = log_pdf_mvhyper(cStates_postDet, [nS, nE, nI])

      permutations[trac_nm, :] = [nS, nE, nI, perm_prob]
      trac_nm += 1
    end
  end

  permutations = permutations[findall(permutations[:,1] .> -1), :]

  choose_perm = wsample(1:size(permutations, 1), permutations[:, 4])

  chosen_perm = [permutations[choose_perm, :] ; log_pdf_mvhyper(cStates_postDet, cDeaths)]

  return(chosen_perm)
end

function propose_AddRem_Deaths_update(;combi_dfs, epi_params, arg...)

  #### Choose a farm and a timestep ####

    # Extract all farms that have deaths for any t and not all animals are S

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        total_deaths = sum(combi_dfs[position_id, 2][timestep , 21:23])
        total_animals = sum(combi_dfs[position_id, 1][timestep , 15:17])
        if  total_deaths > 0. && combi_dfs[position_id, 1][timestep , 15] < total_animals
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = 360

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # Number of death events at time t
    cDeaths_t =  track[t, 21:23]

    # Generate a new number of events
    cStates_postDet_t = res[t, 15:17]

    chosen_perm = deaths_permutations(cStates_postDet_t, cDeaths_t, epi_params)

    new_Sdeaths_t = chosen_perm[1]
    new_Edeaths_t = chosen_perm[2]
    new_Ideaths_t = chosen_perm[3]

    ΔS = new_Sdeaths_t - cDeaths_t[1]
    ΔE = new_Edeaths_t - cDeaths_t[2]
    ΔI = new_Ideaths_t - cDeaths_t[3]

    if (ΔS + ΔE + ΔI) == 0
      log_q_ratio = -Inf
      # println("No Diff")
      AddRem_dths_track = [pos_id, t, 0, 2, cDeaths_t[1], cDeaths_t[2], cDeaths_t[3], new_Sdeaths_t, new_Edeaths_t, new_Ideaths_t, ΔS, ΔE, ΔI, cStates_postDet_t[1], cStates_postDet_t[2], cStates_postDet_t[3]]

      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_dths_track)
    end

  #### Generate new states####

    # Update the events
    track_prime[t , 21] += ΔS
    track_prime[t , 22] += ΔE
    track_prime[t , 23] += ΔI

    # Update the states

    # :cS_final, :pcS_final
    res_prime[ t:t , [18, 35]] .-= ΔS
    # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final, :pcS_init, :pcS_final
    res_prime[((t+1):T) , [3, 6, 9, 12, 15, 18, 32, 35]] .-= ΔS

    # :cE_final, :pcE_final
    res_prime[ t:t , [19, 36]] .-= ΔE
    # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final, :pcE_init, :pcE_final
    res_prime[((t+1):T) , [4, 7, 10, 13, 16, 19, 33, 36]] .-= ΔE

    # :cI_final, :pcI_final
    res_prime[ t:t , [20, 37]] .-= ΔI
    # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final, :pcI_init, :pcI_final
    res_prime[((t+1):T) , [5, 8, 11, 14, 17, 20, 34, 37]] .-= ΔI

    #### Quick check for validity ####

    posi_check = res_prime[(t:T), 3:20] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      log_q_ratio = -Inf
      # println("Invalid")
      AddRem_dths_track = [pos_id, t, 0, 3, cDeaths_t[1], cDeaths_t[2], cDeaths_t[3], new_Sdeaths_t, new_Edeaths_t, new_Ideaths_t, ΔS, ΔE, ΔI, cStates_postDet_t[1], cStates_postDet_t[2], cStates_postDet_t[3]]

      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_dths_track)
    end

    # Calculate changes to probabilities

    β_c = epi_params[1]
    F = epi_params[4]

    for τ in (t:T)

      cI_init = res_prime[τ , 5]
      cN_init = sum(res_prime[τ , 3:5])

      p_env_prev = pers_prime[τ , 3]

      c_exp_prob = 0.

      if cN_init > 0.
        c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
      end

      # Alternative to loop
      # 1 .- exp.( -(cI_init ./ cN_init) * β_c - F * p_env_prev)
      # gives ans of 1 when cN_init = 0

      pers_prime[τ , 5] = c_exp_prob
    end

    # Update all farms in the parish with parish level states

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    pcS_init = res_prime[(t:T) , 32]
    pcE_init = res_prime[(t:T) , 33]
    pcI_init = res_prime[(t:T) , 34]

    pcS_final = res_prime[(t:T) , 35]
    pcE_final = res_prime[(t:T) , 36]
    pcI_final = res_prime[(t:T) , 37]

    for farm in farms_in_parish
      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      combi_dfs_prime[farm,1][(t:T), 32] = pcS_init
      combi_dfs_prime[farm,1][(t:T), 33] = pcE_init
      combi_dfs_prime[farm,1][(t:T), 34] = pcI_init

      combi_dfs_prime[farm,1][(t:T), 35] = pcS_final
      combi_dfs_prime[farm,1][(t:T), 36] = pcE_final
      combi_dfs_prime[farm,1][(t:T), 37] = pcI_final
    end

    # Calculate the log q ratio (proposal density ratio)
    q_prime_given_cur = chosen_perm[4]
    q_cur_given_prime = chosen_perm[5]

    log_q_ratio = q_cur_given_prime - q_prime_given_cur

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = t
    upper_t = T
    h_pos_ids = [pos_id]
    h_element_range = 1:13

    AddRem_dths_track = [pos_id, t, 0, 0, cDeaths_t[1], cDeaths_t[2], cDeaths_t[3], new_Sdeaths_t, new_Edeaths_t, new_Ideaths_t, ΔS, ΔE, ΔI, cStates_postDet_t[1], cStates_postDet_t[2], cStates_postDet_t[3]]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], AddRem_dths_track)
end

# propose_AddRem_Deaths_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#
# for i in 1:1000
#   a,b,c = propose_AddRem_Deaths_update(combi_dfs=combi_dfs, epi_params = epi_params_true)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: Haven't seen an accepted yet, all no difference, too few E and I compared to S


############################################
##### Adding/Removing a Movement event #####
############################################

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

function generate_new_movements(t, T, pos_id, combi_dfs, moves_record)

  ####################
  ### Extract Data ###
  ####################

  # :cS_Moves, :cE_Moves, :cI_Moves
  states = combi_dfs[pos_id, 1][t , 6:8]
  # :sus_off, :exp_off, :inf_off
  moves_off = combi_dfs[pos_id, 2][t , 6:8]

  total_moves = sum(moves_off)
  off_row_id = combi_dfs[pos_id, 1][t , 2]

  ############################################
  ### Extract the movement data for time t ###
  ############################################

  move_record_rows = findall((moves_record[:, 1] .== t) .& (moves_record[:, 2] .== off_row_id))

  # This extracts the movement data at time t for the off farm
  moves_data_cur = moves_record[move_record_rows, :]

  moves_data_prime = deepcopy(moves_data_cur)

  # Extract on_row_ids that recieved animals from this farm
  on_cphs_all = moves_data_cur[:, 3]

  # Extract movements to each farm
  m_farm = moves_data_cur[:, 10]


  ##################################
  ### Generate new movement data ###
  ##################################

  new_move_states = rng_mvhyper(states, total_moves)
  running_move_states = deepcopy(new_move_states)

  # FOR each of the farms that recieved animals
  for j in 1:size(on_cphs_all, 1)

    # IF there is more than one farm that recieved animals from this farm
    if (size(on_cphs_all, 1) > 1)
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

  combi_dfs_prime = deepcopy(combi_dfs)
  combi_dfs_prime[pos_id, 2][t , 6:8] = new_move_states

  moves_record_prime = deepcopy(moves_record)
  moves_record_prime[move_record_rows, :] .= moves_data_prime


  ############################################
  ### Calculate the farm level differences ###
  ############################################

  # :Δ_S_on, :Δ_E_on, :Δ_I_on, :on_row_id, :on_parish_idx
  differences = fill(-10., size(on_cphs_all, 1), 5)

  for j in 1:size(on_cphs_all, 1)

    on_row_id = on_cphs_all[j]

    moves_onto_j_cur = moves_data_cur[j, 7:9]
    moves_onto_j_prime = moves_data_prime[j, 7:9]

    Δ_moves = Array(moves_onto_j_prime) - Array(moves_onto_j_cur)

    differences[j, 1:4] = [Δ_moves ; on_row_id]

    # Find the parish

    for pos in 1:100
      if combi_dfs_prime[pos, 1][t, 2] == on_row_id
        differences[j, 5] = combi_dfs_prime[pos, 3][t, 13] # parish_idx
        break # found parish idx ⟹ move on, next farm
      end
    end

  end

  ###################################
  ### Early return: No difference ###
  ###################################

  if all(differences[:, 1:3] .== 0)
    # returns changed position as a Int instead of a Vector
    return(combi_dfs, moves_record, moves_data_cur, moves_data_prime, 2)
  end

  ##############################################
  ### Calculate the parish level differences ###
  ##############################################

  off_p_idx = combi_dfs_prime[pos_id, 3][t, 13]
  all_affected_parishes = [differences[:, 5]; off_p_idx]
  # Extract a list of the unique parishes that were changed
  unique_parishes = unique(all_affected_parishes)

  # :Δ_pS, :Δ_pE, :Δ_pI, :parish_idx
  parish_differences = fill(0., size(unique_parishes, 1), 4)

  pdi = 1
  for par in unique_parishes

    parish_differences[pdi, 4] = par

    for i in 1:size(differences, 1)
      if differences[i, 5] == par
        parish_differences[pdi, 1:3] += differences[i, 1:3]
      end
    end

    pdi += 1
  end

  for i in 1:size(parish_differences, 1)
    if parish_differences[i, 4] == off_p_idx
      Δ_off = Array(combi_dfs_prime[pos_id, 2][t , 6:8]) - Array(combi_dfs[pos_id, 2][t , 6:8])
      parish_differences[i, 1:3] -= Δ_off
      break
    end
  end


  ############################################
  ### Update the data with the differences ###
  ############################################

  combi_dfs_prime = update_data_after_move(t, T, pos_id, combi_dfs, combi_dfs_prime, differences, parish_differences)

  ###################################
  ### Find all affected positions ###
  ###################################

  changed_ids = [differences[(differences[:, 4] .> 0 ), 4] ; off_row_id] # just farms in cheshire
  changed_pos = fill(-1. , size(changed_ids ,1))
  for i in 1:size(changed_ids ,1)
    for pos in 1:100
      if combi_dfs_prime[pos, 1][t, 2] == changed_ids[i]
        changed_pos[i] = pos
        break
      end
    end
  end

  changed_pos = changed_pos[changed_pos .> 0]
  changed_pos = convert(Vector{Int64}, changed_pos)

  ####################################
  ### Early return: Invalid Update ###
  ####################################

  for chpos in changed_pos
    posi_check = combi_dfs_prime[chpos, 1][(t:T), 1:20] .>= 0

    if sum(sum.(eachrow(posi_check))) != prod(size(posi_check))
      # returns changed position as a Int instead of a Vector
      return(combi_dfs, moves_record, moves_data_cur, moves_data_prime, 3)
    end
  end

  return(combi_dfs_prime, moves_record_prime, moves_data_cur, moves_data_prime, changed_pos)
end


function update_data_after_move(t, T, pos_id, combi_dfs_cur, combi_dfs_prime, differences, parish_differences)

  ##########################################
  ### Update the farm that moved animals ###
  ##########################################

  Δ_off = Array(combi_dfs_prime[pos_id, 2][t , 6:8]) - Array(combi_dfs_cur[pos_id, 2][t , 6:8])

  # :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_dfs_prime[pos_id, 1][t:t, [9, 12, 15, 18]] .-= Δ_off[1]
  # :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_dfs_prime[pos_id, 1][t:t, [10, 13, 16, 19]] .-= Δ_off[2]
  # :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_dfs_prime[pos_id, 1][t:t, [11, 14, 17, 20]] .-= Δ_off[3]

  # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
  combi_dfs_prime[pos_id, 1][(t+1):T, [3, 6, 9, 12, 15, 18]] .-= Δ_off[1]
  # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
  combi_dfs_prime[pos_id, 1][(t+1):T, [4, 7, 10, 13, 16, 19]] .-= Δ_off[2]
  # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
  combi_dfs_prime[pos_id, 1][(t+1):T, [5, 8, 11, 14, 17, 20]] .-= Δ_off[3]


  ##############################################
  ### Update the farms that recieved animals ###
  ##############################################

  for i in 1:size(differences, 1)
    if differences[i, 5] > 0 # Check in Cheshire

      on_row_id = differences[i, 4]

      for pos in 1:100

        if combi_dfs_prime[pos, 1][t, 2] == on_row_id

          Δ_i = differences[i, 1:3]

          # :cS_postM, :cS_postEI, :cS_postDet, :cS_final
          combi_dfs_prime[pos, 1][t:t, [9, 12, 15, 18]] .+= Δ_i[1]
          # :cE_postM, :cE_postEI, :cE_postDet, :cE_final
          combi_dfs_prime[pos, 1][t:t, [10, 13, 16, 19]] .+= Δ_i[2]
          # :cI_postM, :cI_postEI, :cI_postDet, :cI_final
          combi_dfs_prime[pos, 1][t:t, [11, 14, 17, 20]] .+= Δ_i[3]


          # :cS_init, :cS_Moves, :cS_postM, :cS_postEI, :cS_postDet, :cS_final
          combi_dfs_prime[pos, 1][(t+1):T, [3, 6, 9, 12, 15, 18]] .+= Δ_i[1]
          # :cE_init, :cE_Moves, :cE_postM, :cE_postEI, :cE_postDet, :cE_final
          combi_dfs_prime[pos, 1][(t+1):T, [4, 7, 10, 13, 16, 19]] .+= Δ_i[2]
          # :cI_init, :cI_Moves, :cI_postM, :cI_postEI, :cI_postDet, :cI_final
          combi_dfs_prime[pos, 1][(t+1):T, [5, 8, 11, 14, 17, 20]] .+= Δ_i[3]

          # :sus_on, :exp_on, :inf_on
          combi_dfs_prime[pos, 2][t:t, 3] .+= Δ_i[1]
          combi_dfs_prime[pos, 2][t:t, 4] .+= Δ_i[2]
          combi_dfs_prime[pos, 2][t:t, 5] .+= Δ_i[3]

          break # Found the right farm ⟹ move on, next movement
        end
      end # end of pos cycle

    end
  end # end of for each move


  ################################
  ### Update the parish totals ###
  ################################

  for j in 1:size(parish_differences, 1)

    p_idx = parish_differences[j, 4]

    Δ_pj = parish_differences[j, 1:3]

    if p_idx > 0

      for pos in 1:100

        if combi_dfs_prime[pos, 3][t, 13] == p_idx

          # :pcS_final
          combi_dfs_prime[pos, 1][t:t, [35]] .+= Δ_pj[1]
          # :pcE_final
          combi_dfs_prime[pos, 1][t:t, [36]] .+= Δ_pj[2]
          # :pcI_final
          combi_dfs_prime[pos, 1][t:t, [37]] .+= Δ_pj[3]


          # :pcS_init, :pcS_final
          combi_dfs_prime[pos, 1][(t+1):T, [32,35]] .+= Δ_pj[1]
          # :pcE_init, :pcE_final
          combi_dfs_prime[pos, 1][(t+1):T, [33,36]] .+= Δ_pj[2]
          # :pcI_init, :pcI_final
          combi_dfs_prime[pos, 1][(t+1):T, [34,37]] .+= Δ_pj[3]

        end
      end
    end # end if p_idx_on > 0

  end # end of for each parish


  return(combi_dfs_prime)
end


function propose_AddRem_Move_update(;combi_dfs, epi_params, moves_record)

  ####################################
  ### Choose a farm and a timestep ###
  ####################################

    # Extract all farms that moved animals off at t
    # Ensuring moves > 0 and not all animals are S

    choosing_a_farm = []

    for timestep in 1:360
      for position_id in 1:100
        total_moves = sum(combi_dfs[position_id, 2][timestep , 6:8])
        total_animals = sum(combi_dfs[position_id, 1][timestep , 6:8])
        if total_moves > 0. && combi_dfs[position_id, 1][timestep , 3] < total_animals
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  ######################################
  ### Generate the update parameters ###
  ######################################

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = 360

  #######################
  ### Generate update ###
  #######################

    # Number of E to I events at time t
    combi_dfs_prime, moves_record_prime, moves_data_cur, moves_data_prime, changed_pos = generate_new_movements(t, T, pos_id, combi_dfs, moves_record)

  #################################################
  ### Early return: No change or Invalid Update ###
  #################################################

    # If there is no change or the update is invalid for any farm
    # changed_pos will be returned as a Int instead of a Vector
    if typeof(changed_pos) == Int64
      cStates_Moves_cur = combi_dfs[pos_id, 1][t , 6:8]
      cMovesOff_cur = combi_dfs[pos_id, 2][t , 6:8]
      cMovesOff_prime = combi_dfs_prime[pos_id, 2][t , 6:8]

      log_q_ratio = -Inf
      AddRem_Moves_track = [pos_id, t, 0, changed_pos,
                            cMovesOff_cur[1], cMovesOff_cur[2], cMovesOff_cur[3],
                            cMovesOff_prime[1], cMovesOff_prime[2], cMovesOff_prime[3],
                            cStates_Moves_cur[1], cStates_Moves_cur[2], cStates_Moves_cur[3]]

      return(combi_dfs, log_q_ratio, [0,0,0,0], moves_record, AddRem_Moves_track)
    end

    ##########################################
    ### Calculate changes to probabilities ###
    ##########################################

    β_c = epi_params[1]
    F = epi_params[4]

    for pos_inf in changed_pos
      for τ in (t:T)

        cI_init = combi_dfs_prime[pos_inf, 1][τ , 5]
        cN_init = sum(combi_dfs_prime[pos_inf, 1][τ , 3:5])

        p_env_prev = combi_dfs_prime[pos_inf, 3][τ , 3]

        c_exp_prob = 0.

        if cN_init > 0.
          c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
        end

        # Alternative to loop
        # 1 .- exp.( -(cI_init ./ cN_init) * β_c - F * p_env_prev)
        # gives ans of 1 when cN_init = 0

        combi_dfs_prime[pos_inf, 3][τ , 5] = c_exp_prob
      end
    end


    ##########################################################
    ### Calculate the log q ratio (proposal density ratio) ###
    ##########################################################

    cStates_Moves_cur = combi_dfs[pos_id, 1][t , 6:8]
    cMovesOff_cur = combi_dfs[pos_id, 2][t , 6:8]

    cStates_Moves_prime = combi_dfs_prime[pos_id, 1][t , 6:8]
    cMovesOff_prime = combi_dfs_prime[pos_id, 2][t , 6:8]

    q_prime_given_cur = log_pdf_mvhyper(cStates_Moves_cur, cMovesOff_prime)

    q_cur_given_prime = log_pdf_mvhyper(cStates_Moves_prime, cMovesOff_cur)

    for i in 1:size(moves_data_cur, 1)

      cStates_gen_i_cur = moves_data_cur[i, 4:6]
      Moves_i_cur = moves_data_cur[i, 7:9]

      cStates_gen_i_prime = moves_data_prime[i, 4:6]
      Moves_i_prime = moves_data_prime[i, 7:9]

      q_prime_given_cur += log_pdf_mvhyper(cStates_gen_i_prime, Moves_i_prime)
      q_cur_given_prime += log_pdf_mvhyper(cStates_gen_i_cur, Moves_i_cur)

    end

    log_q_ratio = q_cur_given_prime - q_prime_given_cur

    ##############
    ### Return ###
    ##############

    lower_t = t
    upper_t = T
    h_pos_ids = changed_pos
    h_element_range = 1:14

    AddRem_Moves_track = [pos_id, t, 0, 0,
                          cMovesOff_cur[1], cMovesOff_cur[2], cMovesOff_cur[3],
                          cMovesOff_prime[1], cMovesOff_prime[2], cMovesOff_prime[3],
                          cStates_Moves_cur[1], cStates_Moves_cur[2], cStates_Moves_cur[3]]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], moves_record_prime, AddRem_Moves_track)
end

# propose_AddRem_Move_update(combi_dfs=combi_dfs, epi_params = epi_params_true, moves_record = record_of_movements)
#
# for i in 1:1000
#   a,b,c,d = propose_AddRem_Move_update(combi_dfs=combi_dfs, epi_params = epi_params_true, moves_record = record_of_movements)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: 10% - 20% of updates are valid and novel



#########################################
##### Adding/Removing a p-env event #####
#########################################

# p_env at time t in parish is made up of some amount of pressure from time t-1
# plus some new pressure. We will update both components at once.
# It can be thought of adding or removing pressure.

# Need to change ϵ to 1-ϵ when updated sim

function propose_AddRem_penv_update(;combi_dfs, epi_params, scaling)

  #### Choose a farm and a timestep ####

    # Extract all farms that have p_env_t > 0 or pI_t > 0 for any t
    # Bigger parishes more likely to be updated

    choosing_a_farm = []

    for timestep in 1:359 # Not 360 because changing 360 makes almost no difference
      for position_id in 1:100
        if combi_dfs[position_id, 3][timestep , 3] > 0 || sum(combi_dfs[position_id, 1][timestep , [34,40]]) > 0
          push!(choosing_a_farm, [position_id, timestep])
        end
      end
    end

    choose_farm_time = rand(1:size(choosing_a_farm, 1))

  #### Generate the update parameters ####

    pos_id = choosing_a_farm[choose_farm_time][1]
    t = choosing_a_farm[choose_farm_time][2]
    T = t+1

  #### Extract the farm data ####

    res = deepcopy(combi_dfs[pos_id, 1])
    track = deepcopy(combi_dfs[pos_id, 2])
    pers = deepcopy(combi_dfs[pos_id, 3])

    combi_dfs_prime = deepcopy(combi_dfs)

    res_prime = deepcopy(res)
    track_prime = deepcopy(track)
    pers_prime = deepcopy(pers)

  #### Generate update ####

    # p_env events at time t
    p_env_prev_t = pers[t, 3]

    remaining_pressure_t = track[t, 28]
    new_pressure_t = track[t, 29]

    # States at time t (pcI and pbI)
    pI_init_t = sum(res[t, [34,40]])

    # Generate a new number of events
    ϵ = epi_params[5]
    c_idxx = convert(Int64, pers[1, 12])
    p_idxx = convert(Int64, pers[1, 13])
    scaling_cp = scaling[c_idxx][p_idxx]

    propose_r_pressure_t = rand(Binomial(round(Int, scaling_cp * p_env_prev_t), ϵ)) # should be 1 - ϵ
    propose_n_pressure_t = rand(Poisson(pI_init_t))

    Δr = propose_r_pressure_t - remaining_pressure_t
    Δn = propose_n_pressure_t - new_pressure_t

    if (Δr+Δn) == 0
      log_q_ratio = -Inf
      AddRem_penv_track = [pos_id, t, 0, 2, remaining_pressure_t, new_pressure_t, propose_r_pressure_t, propose_n_pressure_t, Δr, Δn, p_env_prev_t, pI_init_t]

      return(combi_dfs, log_q_ratio, [0,0,0,0], AddRem_penv_track)
    end

  #### Generate new states####

    # Update the events
    track_prime[t , 28] += Δr
    track_prime[t , 29] += Δn

    Δpenv = (Δr + Δn)/scaling_cp

    # Update the states

    # :p_env_cur
    pers_prime[t, 4] += Δpenv
    # :p_env_prev
    pers_prime[T, 3] += Δpenv


    #### Extract all the farms in the parish

    parish_idx = pers_prime[1, 13]
    farms_in_parish = fill(false, size(combi_dfs, 1))
    for pos in 1:size(combi_dfs, 1)
      p_idx = combi_dfs_prime[pos, 3][1,13]
      farms_in_parish[pos] = (p_idx == parish_idx)
    end
    farms_in_parish = findall(farms_in_parish .== true)

    # Calculate changes to probabilities

    β_c = epi_params[1]
    β_b = epi_params[2]
    F = epi_params[4]

    r_pres_prime = track_prime[t , 28]
    n_pres_prime = track_prime[t , 29]
    penv_cur_prime = pers_prime[t, 4]
    penv_prev_prime = pers_prime[T, 3]


    for farm in farms_in_parish
      # Record for each holding in the parish, the updated values of p_env
      # and the p_env events
      combi_dfs_prime[farm,2][t, 28] = r_pres_prime
      combi_dfs_prime[farm,2][t, 29] = n_pres_prime
      combi_dfs_prime[farm,3][t, 4] = penv_cur_prime
      combi_dfs_prime[farm,3][T, 3] = penv_prev_prime

      # Recalculate the probabilities of infection on each farm

      cI_init = combi_dfs_prime[farm,1][T , 5]
      cN_init = sum(combi_dfs_prime[farm,1][T , 3:5])

      c_exp_prob = 0.

      if cN_init > 0.
        c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * penv_prev_prime)
      end

      combi_dfs_prime[farm,3][T, 5] = c_exp_prob

      bI_init = combi_dfs_prime[farm,1][T , 23]
      bN_init = sum(combi_dfs_prime[farm,1][T , 21:23])

      b_exp_prob = 0.

      if bN_init > 0.
        b_exp_prob = 1 - exp( - bI_init/bN_init * β_b - F * penv_prev_prime)
      end

      combi_dfs_prime[farm,3][T, 6] = b_exp_prob
    end

    # Calculate the log q ratio (proposal density ratio)

    q_prime_r_press_given_cur = logpdf(Binomial(round(Int, scaling_cp * p_env_prev_t), ϵ), propose_r_pressure_t) # should be 1 - ϵ
    q_prime_n_press_given_cur = logpdf(Poisson(pI_init_t), propose_n_pressure_t)

    q_cur_r_press_given_prime = logpdf(Binomial(round(Int, scaling_cp * p_env_prev_t), ϵ), remaining_pressure_t) # should be 1 - ϵ
    q_cur_n_press_given_prime = logpdf(Poisson(pI_init_t), new_pressure_t)

    log_q_ratio = (q_cur_r_press_given_prime + q_cur_n_press_given_prime) - (q_prime_r_press_given_cur + q_prime_n_press_given_cur)

    # Replace arrays in combi_dfs

    combi_dfs_prime[pos_id, 1] = res_prime
    combi_dfs_prime[pos_id, 2] = track_prime
    combi_dfs_prime[pos_id, 3] = pers_prime

    lower_t = t
    upper_t = T
    h_pos_ids = farms_in_parish
    h_element_range = [3,8]

    AddRem_penv_track = [pos_id, t, 0, 0, remaining_pressure_t, new_pressure_t, propose_r_pressure_t, propose_n_pressure_t, Δr, Δn, p_env_prev_t, pI_init_t]

  return(combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range], AddRem_penv_track)
end

# propose_AddRem_penv_update(combi_dfs=combi_dfs, epi_params = epi_params_true, scaling = area_of_parish)
#
# for i in 1:1000
#   a,b,c = propose_AddRem_penv_update(combi_dfs=combi_dfs, epi_params = epi_params_true, scaling = area_of_parish)
#   if b != -Inf println(i, " ", b) end
# end

# NOTE: All valid, but can't update t = 360




#########################################
##### Updating Infection Parameters #####
#########################################

# TODO FIXME CHANGED XXX IDEA HACK NOTE REVIEW NB BUG QUESTION COMBAK TEMP

# NOTE: Either remove or move to end the "sample" column in the results to make indexing easier
# removes the block parameter and just uses params_oi

function update_combi_dfs(combi_dfs, epi_params)

  combi_dfs_prime = deepcopy(combi_dfs)

  β_c = epi_params[1]
  β_b = epi_params[2]
  γ = epi_params[3]
  F = epi_params[4]

  for farm in 1:100
    for τ in 1:360

      p_env_prev_prime = combi_dfs_prime[farm,3][τ, 3]

      # Calculate cattle exposure probability

      cStates_init = combi_dfs_prime[farm,1][τ , 3:5]

      c_exp_prob = calc_exp_prop(States_init = cStates_init,
                                 p_env_prev = p_env_prev_prime,
                                 β = β_c, F = F)

      combi_dfs_prime[farm,3][τ, 5] = c_exp_prob

      # Calculate badger exposure probability

      bStates_init = combi_dfs_prime[farm,1][τ , 21:23]

      b_exp_prob = calc_exp_prop(States_init = bStates_init,
                                 p_env_prev = p_env_prev_prime,
                                 β = β_b, F = F)

      combi_dfs_prime[farm,3][τ, 6] = b_exp_prob

      # Calculate infection probability

      inf_prob = 1 - exp(-γ)

      combi_dfs_prime[farm,3][τ, 7] = inf_prob
      combi_dfs_prime[farm,3][τ, 8] = inf_prob

    end
  end

  return(combi_dfs_prime)
end


function tune_λ(;λ_cur, it, n_tune, other_res, acc_sum_idx)

    log_λ = log(λ_cur)

    ## Calculate acceptance rate of last batch ##
    batch_start = ((n_tune - 1) * 25) + 1
    batch_end = (n_tune*25) - 1

    bot = batch_end - batch_start
    top = sum(other_res[batch_start:batch_end, acc_sum_idx])

    acc_prop = top/bot

    ## Update tuning parameter ##
    if acc_prop < 0.33
      delta_n = -1 * min(0.05, 1/sqrt(n_tune))
    else
      delta_n = min(0.05, 1/sqrt(n_tune))
    end

    log_λ = log_λ + delta_n

    λ_cur = exp(log_λ)

  return(λ_cur)
end # Don't need it

function Blk_update(;N_its, combi_dfs, results, other_res,
                      it, params_oi, log_epi_params_cur,
                      n_tune, acc_sum_idx, m, λ)

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      d = size(params_oi, 1)

      log_epi_params_draw = deepcopy(log_epi_params_cur)

      log_epi_params_oi_cur = log_epi_params_cur[params_oi]


  ######################
  ### Initial tuning ###
  ######################

      if it < min(5000, (N_its/10))

        if it == (n_tune*25)
          λ = tune_λ(λ_cur = λ, it = it, n_tune = n_tune, other_res = other_res, acc_sum_idx = acc_sum_idx)
        end

        fixed_scalar = (1/d) * (λ)^2

        log_epi_params_oi_draw = rand(MvNormal(log_epi_params_oi_cur, fixed_scalar))

        mixture = 0
      end

  ##########################
  ### Draw the paramters ###
  ##########################

      if it >= min(5000, (N_its/10))

        mixture = rand()

        ### Non-adaptive draw ###

        if mixture < 0.05

          fixed_scalar = (1/d) * (λ)^2

          log_epi_params_oi_draw = rand(MvNormal(log_epi_params_oi_cur, fixed_scalar))
        end


        ### Adaptive Draw ###

        if mixture >= 0.05

          log_emp_post = Array(results[1:(it-1), params_oi])

          Σ = cov(log_emp_post, dims = 1)

          scalarΣ = (m)^2 * Σ

          log_epi_params_oi_draw = rand(MvNormal(log_epi_params_oi_cur, scalarΣ))
        end

      end

    log_epi_params_draw[params_oi] = log_epi_params_oi_draw

  return(mixture, λ, log_epi_params_draw)
end


function calc_exp_prop(;States_init, p_env_prev, β, F)

  I = States_init[3]
  N = sum(States_init)

  exp_prob = 0.

  if N > 0.
    exp_prob = 1 - exp( - I/N * β - F * p_env_prev)
  end

  return(exp_prob)
end


function propose_infection_params_update(;N_its, combi_dfs, res, other_res,
                      it, log_epi_params_cur,
                      n_tune, m, λ)

  # Generate an update for the infection parameters

    mixture, λ, log_epi_params_draw =
      Blk_update(N_its = N_its, combi_dfs = combi_dfs, results = res,
                other_res = other_res, it = it, params_oi = [1,2,3,4],
                log_epi_params_cur = log_epi_params_cur,
                n_tune = n_tune, acc_sum_idx = 5,
                m = m, λ = λ)

  # Calculate changes to probabilities

    combi_dfs_prime = deepcopy(combi_dfs)

    epi_params_draw = inverselogtransform(log_epi_params_draw)

    β_c = epi_params_draw[1]
    β_b = epi_params_draw[2]
    γ = epi_params_draw[3]
    F = epi_params_draw[4]

    for farm in 1:100
      for τ in 1:360

        p_env_prev_prime = combi_dfs_prime[farm,3][τ, 3]

        # Calculate cattle exposure probability

        cStates_init = combi_dfs_prime[farm,1][τ , 3:5]

        c_exp_prob = calc_exp_prop(States_init = cStates_init,
                                   p_env_prev = p_env_prev_prime,
                                   β = β_c, F = F)

        combi_dfs_prime[farm,3][τ, 5] = c_exp_prob

        # Calculate badger exposure probability

        bStates_init = combi_dfs_prime[farm,1][τ , 21:23]

        b_exp_prob = calc_exp_prop(States_init = bStates_init,
                                   p_env_prev = p_env_prev_prime,
                                   β = β_b, F = F)

        combi_dfs_prime[farm,3][τ, 6] = b_exp_prob

        # Calculate infection probability

        inf_prob = 1 - exp(-γ)

        combi_dfs_prime[farm,3][τ, 7] = inf_prob
        combi_dfs_prime[farm,3][τ, 8] = inf_prob

      end
    end

    # Calculate the log q ratio (proposal density ratio)

    log_q_ratio = 0

    # Scope

    lower_t = 1
    upper_t = 360
    h_pos_ids = 1:100
    h_element_range = [3,4,8,9]

  return(mixture, λ, log_epi_params_draw, combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range])
end

# propose_infection_params_update(N_its = 100000, combi_dfs = combi_dfs,
#                                 results = res, other_res = other_res,
#                                 it = 1, log_epi_params_cur = log_epi_params_cur,
#                                 n_tune = 1, m = m_inf, λ = λ_inf)



#########################################
##### Updating Infection Parameters #####
#########################################



function propose_detection_params_update(;N_its, combi_dfs, res, other_res,
                                          it, log_epi_params_cur,
                                          n_tune, m, λ)

  # Generate an update for the infection parameters

    mixture, λ, log_epi_params_draw =
      Blk_update(N_its = N_its, combi_dfs = combi_dfs, results = res,
                other_res = other_res, it = it, params_oi = [6,7],
                log_epi_params_cur = log_epi_params_cur,
                n_tune = n_tune, acc_sum_idx = 6,
                m = m, λ = λ)

  # Calculate changes to probabilities

    combi_dfs_prime = deepcopy(combi_dfs)

    epi_params_draw = inverselogtransform(log_epi_params_draw)

    ρ = epi_params_draw[6]
    ρ_E = epi_params_draw[7]

    # Calculate the log q ratio (proposal density ratio)

    log_q_ratio = 0

    # Scope

    lower_t = 1
    upper_t = 360
    h_pos_ids = 1:100
    h_element_range = [5,6]

  return(mixture, λ, log_epi_params_draw, combi_dfs_prime, log_q_ratio, [lower_t, upper_t, h_pos_ids, h_element_range])
end

# propose_detection_params_update(N_its = 100000, combi_dfs = combi_dfs,
#                                 results = res, other_res = other_res,
#                                 it = 1, log_epi_params_cur = log_epi_params_cur,
#                                 n_tune = 1, m = m_det, λ = λ_det)











#################
### OLD STUFF ###
#################





function θ_bb_gibbs(obs, track)

    num_time_steps = size(obs, 1)

    births = fill(0., num_time_steps)
    N = fill(0., num_time_steps)

    for t in 1:num_time_steps

      bStates_postEI_t = obs[t,[:bS_postEI, :bE_postEI, :bI_postEI]]

      b_new_birth_t = track[t, :b_birth] # Number of new susceptible cattle births

      births[t] = b_new_birth_t
      N[t] = sum(bStates_postEI_t)
    end

    sum_births = sum(births)
    sum_N = sum(N)

    return([sum_births, sum_N])
end

function θ_bb_gibbs_all(combi_dfs, hpam1, hpam2)

    num_dfs = size(combi_dfs, 1)

    sum_births_all = fill(0., num_dfs)
    sum_N_all = fill(0., num_dfs)

    for i in 1:num_dfs

      obs_i = combi_dfs[i, 1]
      track_i = combi_dfs[i, 2]

      sums_i = θ_bb_gibbs(obs_i, track_i)

      births_all[i] = sums_i[1]
      N_all[i] = sums_i[2]
    end

    sum_births = sum(births_all)
    sum_N = sum(N_all)

    d_t = Gamma((sum_births + hpam1) ,1/(sum_N + 1/hpam2))

    θ_bb_draw = rand(d_t)

    return(θ_bb_draw)
end

function θ_bd_gibbs(obs, track)

    num_time_steps = size(obs, 1)

    deaths = fill(0., num_time_steps)
    N = fill(0., num_time_steps)

    for t in 1:num_time_steps

      bStates_postEI_t = obs[t,[:bS_postEI, :bE_postEI, :bI_postEI]]

      bDeaths_t = track[t,[:b_S_death, :b_E_death, :b_I_death]]

      deaths[t] = sum(bDeaths_t)
      N[t] = sum(bStates_postEI_t)
    end

    sum_deaths = sum(deaths)
    sum_N = sum(N)

    return([sum_deaths, sum_N])
end

function θ_bd_gibbs_all(combi_dfs, hpam1, hpam2)

    num_dfs = size(combi_dfs, 1)

    sum_deaths_all = fill(0., num_dfs)
    sum_N_all = fill(0., num_dfs)

    for i in 1:num_dfs

      obs_i = combi_dfs[i, 1]
      track_i = combi_dfs[i, 2]

      sums_i = θ_bd_gibbs(obs_i, track_i, hpam1, hpam2)

      sum_deaths_all[i] = sums_i[1]
      sum_N_all[i] = sums_i[2]
    end

    sum_deaths = sum(sum_deaths_all)
    sum_N = sum(sum_N_all)

    d_t = Beta((sum_deaths + hpam1) ,((sum_N - sum_deaths) + hpam2))

    θ_bd_draw = rand(d_t)

    return(θ_bd_draw)
end
