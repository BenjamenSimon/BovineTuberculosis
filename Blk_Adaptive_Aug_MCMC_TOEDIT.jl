using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames
using DataFramesMeta
using LinearAlgebra
using CSV

function InitialiseAll(;combi_dfs, epi_params_true, epi_params_dists, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # Current values
  β_c_cur, β_b_cur, γ_cur, F_cur, ϵ_cur, ρ_cur, ρ_E_cur, θ_bb_cur, θ_bd_cur = fill(missing, 9)

  epi_params_draw = fill(missing, 9)

  # True values
  β_c_true, β_b_true, γ_true, F_true, ϵ_true, ρ_true, ρ_E_true, θ_bb_true, θ_bd_true = epi_params_true

  # Prior distributions
  d_β_c, d_β_b, d_γ, d_F, d_ϵ, d_ρ, d_ρ_E = epi_params_dists


  arb = true
  while arb == true
    # Initialise parameters
    if ismissing(β_c_true)
      β_c_cur = rand(d_β_c)
    else
      β_c_cur = β_c_true
    end
    if ismissing(β_b_true)
      β_b_cur = rand(d_β_b)
    else
      β_b_cur = β_b_true
    end
    if ismissing(γ_true)
      γ_cur = rand(d_γ)
    else
      γ_cur = γ_true
    end
    if ismissing(F_true)
      F_cur = rand(d_F)
    else
      F_cur = F_true
    end
    if ismissing(ϵ_true)
      ϵ_cur = rand(d_ϵ)
    else
      ϵ_cur = ϵ_true
    end
    if ismissing(ρ_true)
      ρ_cur = rand(d_ρ)
    else
      ρ_cur = ρ_true
    end
    if ismissing(ρ_E_true)
      ρ_E_cur = rand(d_ρ_E)
    else
      ρ_E_cur = ρ_E_true
    end
    if ismissing(θ_bb_true)
      θ_bb_cur = rand(d_θ_bb)
    else
      θ_bb_cur = θ_bb_true
    end
    if ismissing(θ_bd_true)
      θ_bd_cur = rand(d_θ_bd)
    else
      θ_bd_cur = θ_bd_true
    end

    # Current draws
    epi_params_draw = [β_c_cur, β_b_cur, γ_cur, F_cur, ϵ_cur, ρ_cur, ρ_E_cur, θ_bb_cur, θ_bd_cur]

    # Update the probabilities

    combi_dfs_prime = update_combi_dfs(combi_dfs, epi_params_draw)

    #Calculate the log-likelihood
    llh_array_cur = Array{Array{BigFloat, 2}, 1}(undef,100)
    for k in 1:100
      llh_array_cur[k] = Array{BigFloat, 2}(undef, 360, 14)
    end

    p_env_llh_array_cur = Array{Array{BigFloat, 2}, 1}(undef,24)
    for k in 1:24
      p_env_llh_array_cur[k] = Array{BigFloat, 2}(undef, 360, 2)
    end

    llh_array_cur = calc_llh_array_allh(calc_llh_array_row_all, [1, 360, 1:100, 1:13], combi_dfs_prime, llh_array_cur, epi_params_draw)
    p_env_llh_array_cur = calc_p_env_llh_array_allp([1, 360, 1:100, 1:13], p_env_llh_array_cur, combi_dfs_prime, epi_params_draw, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

    init_llh = calc_llh_h_and_p([1, 360, 1:100, 1:13], llh_array_cur, p_env_llh_array_cur)

    if init_llh != -Inf
      break
    end
  end

  println("Initialised")


  return(epi_params_draw)
end

#####################
### Adaptive MCMC ###
#####################

#### Transforms ####

function logtransform(input)

  # A simple log transform on the posterior can cause problems
  # From Examples of Adaptive MCMC (2008) Roberts & Rosenthal
  # we instead us

  # transform = sign.(input) .* log.(1 .+ abs.(input))

  transform = log.(input)

  return(transform)
end

function inverselogtransform(input)

  # transform = sign.(input) .* (exp.(abs.(input)) .- 1)

  transform = exp.(input)

  return(transform)
end

#### Metropolis Hastings functions ####

function mult_MH_accept(;post, post_prime, log_q_ratio, log_epi_params_cur, log_epi_params_draw)
  alpha =  (sum(log_epi_params_draw) + post_prime) - (sum(log_epi_params_cur) + post) + log_q_ratio

  return( min(0, alpha) )
end

function mh_accept_ratio(post, post_prime, log_q_ratio)

  log_α_ratio =  min( 0,  (post_prime) - (post) + log_q_ratio )

  return(log_α_ratio)
end

# MH function for general data augmentation
function metropolis_hastings_step_aug(proposal_func, posterior_func, combi_dfs_cur, llh_array_cur, epi_params_cur, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # Propose an update
  combi_dfs_prime, log_q_ratio, update_scope, update_tracker = proposal_func(combi_dfs=combi_dfs_cur, epi_params = epi_params_cur, scaling = scaling)

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == false
    is_accepted = false
    log_α_ratio = -Inf
    post = 0
    post_prime = -Inf
    return(combi_dfs_cur, llh_array_cur, p_env_llh_array_cur, [is_accepted, log_α_ratio, post, post_prime], update_tracker)
  end

  # Calculate new llh arrays and posteriors
  llh_array_prime, p_env_llh_array_prime, post, post_prime = posterior_func(combi_dfs_prime, update_scope, llh_array_cur, epi_params_cur, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post, post_prime ,log_q_ratio)
  is_accepted = false

  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    combi_dfs_cur = deepcopy(combi_dfs_prime)
    llh_array_cur = deepcopy(llh_array_prime)
    p_env_llh_array_cur = deepcopy(p_env_llh_array_prime)

    update_tracker[3:4] = [1.,1.]

    is_accepted = true
  end

  return(combi_dfs_cur, llh_array_cur, p_env_llh_array_cur, [is_accepted, log_α_ratio, post, post_prime], update_tracker)
end

# MH function for movements data augmentation
function metropolis_hastings_step_moves(combi_dfs_cur, llh_array, epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling, moves_record)

  # Propose an update
  combi_dfs_prime, log_q_ratio, update_scope, moves_record_prime, update_tracker = propose_AddRem_Move_update(combi_dfs=combi_dfs_cur, epi_params = epi_params_cur, moves_record = moves_record)

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == false
    is_accepted = false
    log_α_ratio = -Inf
    post = 0
    post_prime = -Inf
    return(combi_dfs_cur, llh_array, p_env_llh_array, [is_accepted, log_α_ratio, post, post_prime], moves_record, update_tracker)
  end

  # Calculate new llh arrays and posteriors
  llh_array_prime, p_env_llh_array_prime, post, post_prime = AddRem_Moves_posterior(combi_dfs_prime, update_scope, llh_array, epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling, moves_record_prime)

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post, post_prime ,log_q_ratio)

  is_accepted = false

  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    combi_dfs = deepcopy(combi_dfs_prime)
    llh_array = deepcopy(llh_array_prime)
    p_env_llh_array = deepcopy(p_env_llh_array_prime)
    moves_record = deepcopy(moves_record_prime)

    update_tracker[3:4] = [1.,1.]

    is_accepted = true
  end

  return(combi_dfs_cur, llh_array, p_env_llh_array, [is_accepted, log_α_ratio, post, post_prime], moves_record, update_tracker)
end

# MH function for parameters
function metropolis_hastings_step_params(proposal_func, posterior_func, combi_dfs_cur,
                                      llh_array_cur, p_env_llh_array_cur, log_prior_cur, log_prior_dists,
                                      N_its, res, other_res,
                                      it, log_epi_params_cur,
                                      n_tune, m, λ,
                                      first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # Propose an update
  mixture, λ_new, log_epi_params_draw, combi_dfs_prime, log_q_ratio, update_scope = proposal_func(N_its = N_its, combi_dfs = combi_dfs_cur, res = res, other_res = other_res,
                                                             it = it, log_epi_params_cur = log_epi_params_cur,
                                                             n_tune = n_tune, m = m, λ = λ)

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == false
    is_accepted = false
    log_α_ratio = -Inf
    post = 0
    post_prime = -Inf
    return(combi_dfs_cur, llh_array_cur, p_env_llh_array_cur, [is_accepted, log_α_ratio, post, post_prime], mixture, λ_new)
  end

  # Calculate new llh arrays and posteriors

  epi_params_draw = inverselogtransform(log_epi_params_draw)

  llh_array_prime, p_env_llh_array_prime, log_prior_prime, post, post_prime =
       posterior_func(combi_dfs_prime, update_scope, llh_array_cur, epi_params_draw, log_prior_cur, log_prior_dists, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # Calculate MH acceptance probability
  log_α_ratio = mult_MH_accept(post = post, post_prime = post_prime, log_q_ratio = log_q_ratio, log_epi_params_cur = log_epi_params_cur, log_epi_params_draw = log_epi_params_draw)
  is_accepted = false

  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    combi_dfs_cur = deepcopy(combi_dfs_prime)
    llh_array_cur = deepcopy(llh_array_prime)
    p_env_llh_array_cur = deepcopy(p_env_llh_array_prime)

    log_epi_params_cur = deepcopy(log_epi_params_draw)
    log_prior_cur = deepcopy(log_prior_prime)

    is_accepted = true
  end

  return(combi_dfs_cur, llh_array_cur, p_env_llh_array_cur, log_epi_params_cur, log_prior_cur, [is_accepted, log_α_ratio, post, post_prime], mixture, λ_new)
end


#### Adaptive MCMC function ####

function Blk_Adaptive_MCMC_post_all(;N_its, combi_dfs,
                                  epi_params_dists, epi_params_true,
                                  infer_block, λs_init, ms_init,
                                  scaling, first_h_of_p_pos, p_idx_first_h_of_p,
                                  moves_record, data_aug_infer)

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      epi_params_cur = InitialiseAll(;combi_dfs = combi_dfs,
                                      epi_params_true = epi_params_true,
                                      epi_params_dists = epi_params_dists,
                                      first_h_of_p_pos = first_h_of_p_pos,
                                      p_idx_first_h_of_p = p_idx_first_h_of_p,
                                      scaling = scaling)

      combi_dfs = update_combi_dfs(combi_dfs, epi_params_cur)

      epi_params_draw = deepcopy(epi_params_cur)

      log_epi_params_cur = logtransform(epi_params_cur)

  ######################
  ### Results matrix ###
  ######################

      res = Array{Float64}(undef, N_its, 10)
      # :β_c, :β_b, :γ, :F, :ϵ, :ρ, :ρ_E, :θ_bb, :θ_bd, :sample

      other_res = Array{Float64}(undef, N_its, 13)
      # :λ_inf, :λ_det, :m_inf, :m_det,
      # :acc_inf, :acc_det,
      # :log_α_ratio_inf, :log_α_ratio_det,
      # :post_inf, :post_det,
      # :post_prime_inf, :post_prime_det,
      # :sample

      aug_res = Array{Float64}(undef, N_its, 32)
      # [is_accepted, log_α_ratio, post, post_prime]
      # for
      # Move SE
      # Move EI
      # Add/Rem SE
      # Add/Rem EI
      # Add/Rem Det
      # Add/Rem Death
      # Add/Rem penv
      # Add/Rem Moves

      move_SE_tracker = Array{Float64}(undef, N_its, 6)
      # :pos_id, :t, :is_accepted, :reason,
      # :Δ_time, :num_moved

      move_EI_tracker = Array{Float64}(undef, N_its, 6)
      # :pos_id, :t, :is_accepted, :reason,
      # :Δ_time, :num_moved

      AddRem_SE_tracker = Array{Float64}(undef, N_its, 9)
      # :pos_id, :t, :is_accepted, :reason,
      # :SE_before, :SE_after, :Δ_diff, :cS, :prob

      AddRem_EI_tracker = Array{Float64}(undef, N_its, 9)
      # :pos_id, :t, :is_accepted, :reason,
      # :EI_before, :EI_after, :Δ_diff, :cE, :prob

      AddRem_Dets_tracker = Array{Float64}(undef, N_its, 12)
      # :pos_id, :t, :is_accepted, :reason,
      # :Edet_before, :Idet_before, :Edet_after, :Idet_after,
      # :ΔE_diff, :ΔI_diff, :cE, :cI

      AddRem_Deaths_tracker = Array{Float64}(undef, N_its, 16)
      # :pos_id, :t, :is_accepted, :reason,
      # :Sdths_before, :Edths_before, :Idths_before,
      # :Sdths_after, :Edths_after, :Idths_after,
      # :ΔS_diff, :ΔE_diff, :ΔI_diff, :cS, :cE, :cI

      AddRem_penv_tracker = Array{Float64}(undef, N_its, 12)
      # :pos_id, :t, :is_accepted, :reason,
      # :r_pres_before, :n_pres_before, :r_pres_after, :n_pres_after,
      # :Δr_diff, :Δn_diff, :p_env_prev, :pI

      AddRem_Moves_tracker = Array{Float64}(undef, N_its, 13)
      # :pos_id, :t, :is_accepted, :reason,
      # :Soff_before, :Eoff_before, :Ioff_before,
      # :Soff_after, :Eoff_after, :Ioff_after,
      # :cS, :cE, :cI

  ##########################
  ### Functional objects ###
  ##########################

      it = 1

      inf_log_prior_cur = Infection_prior(epi_params_dists, epi_params_cur)
      det_log_prior_cur = Detection_prior(epi_params_dists, epi_params_cur)


      llh_array = Array{Array{BigFloat, 2}, 1}(undef,100)
      for k in 1:100
        llh_array[k] = Array{BigFloat, 2}(undef, 360, 14)
      end

      p_env_llh_array = Array{Array{BigFloat, 2}, 1}(undef,24)
      for k in 1:24
        p_env_llh_array[k] = Array{BigFloat, 2}(undef, 360, 2)
      end

      llh_array = calc_llh_array_allh(calc_llh_array_row_everything, [1, 360, 1:100, 1:14], combi_dfs, llh_array, epi_params_cur, moves_record)
      p_env_llh_array = calc_p_env_llh_array_allp([1, 360, 1:100, 1:14], p_env_llh_array, combi_dfs, epi_params_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

      λ_inf, λ_det = λs_init
      m_inf, m_det = ms_init
      Δ_inf, Δ_det = ms_init/100

      n_tune = 1

  ###########################
  ### ~~ THE ALGORITHM ~~ ###
  ###########################

      while it <= N_its

        ##############################
        ### Gibbs sampler for θ_bb ###
        ##############################

        if infer_block[3] == true
            # Draw θ_bb
            θ_bb_cur = θ_bb_gibbs_all(combi_dfs, (0.25/52), 1)

            epi_params_cur[9] = θ_bb_cur
        end


        ##############################
        ### Gibbs sampler for θ_bd ###
        ##############################

        if infer_block[4] == true
            # Draw θ_bd
            θ_bd_cur = θ_bd_gibbs_all(combi_dfs, 0.000304037, 0.0629356)

            epi_params_cur[10] = θ_bd_cur
        end


        ########################################
        ### MH step for Infection parameters ###
        ########################################

        mh_res_inf = [-Inf, -Inf, Inf, Inf]

        if infer_block[1] == true

          combi_dfs, llh_array, p_env_llh_array, log_epi_params_cur, inf_log_prior_cur, mh_res_inf, mixture, λ_inf =
            metropolis_hastings_step_params(propose_infection_params_update, Infection_params_posterior,
                              combi_dfs, llh_array, p_env_llh_array, inf_log_prior_cur, epi_params_dists,
                              N_its, res, other_res, it, log_epi_params_cur,
                              n_tune, m_inf, λ_inf, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


          if mixture >= 0.05 # mixture is set to 0 while tuning λ
            if mh_res_inf[1] == 0
              m_inf = m_inf - (Δ_inf/((it)^0.5))
            else
              m_inf = m_inf + 2.3*(Δ_inf/((it)^0.5))
            end
          end

        end


        ########################################
        ### MH step for Detection parameters ###
        ########################################

        mh_res_det = [-Inf, -Inf, Inf, Inf]

        if infer_block[2] == true

          combi_dfs, llh_array, p_env_llh_array, log_epi_params_cur, det_log_prior_cur, mh_res_det, mixture, λ_det =
            metropolis_hastings_step_params(propose_detection_params_update, Detection_params_posterior,
                              combi_dfs, llh_array, p_env_llh_array, det_log_prior_cur, epi_params_dists,
                              N_its, res, other_res, it, log_epi_params_cur,
                              n_tune, m_det, λ_det, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

            if mixture >= 0.05 # mixture is set to 0 while tuning λ
              if mh_res_det[1] == 0
                m_det = m_det - (Δ_det/((it)^0.5))
              else
                m_det = m_det + 2.3*(Δ_det/((it)^0.5))
              end
            end

        end


        ###############################
        ### Data Augmentation Steps ###
        ###############################

        epi_params_cur = inverselogtransform(log_epi_params_cur)

        ###################################
        ### Move S→E event through time ###
        ###################################

        mh_res_move_SE = [-Inf, -Inf, Inf, Inf]
        move_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[1] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_move_SE, move_SE_track =
            metropolis_hastings_step_aug(propose_Move_SE_update, Move_SE_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

        end

        ###################################
        ### Move E→I event through time ###
        ###################################

        mh_res_move_EI = [-Inf, -Inf, Inf, Inf]
        move_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[2] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_move_EI, move_EI_track =
            metropolis_hastings_step_aug(propose_Move_EI_update, Move_EI_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

        end

        ############################
        ### Add/Remove S→E event ###
        ############################

        mh_res_AddRem_SE = [-Inf, -Inf, Inf, Inf]
        AddRem_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[3] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_SE, AddRem_SE_track =
            metropolis_hastings_step_aug(propose_AddRem_SE_update, AddRem_SE_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)
        end

        ############################
        ### Add/Remove E→I event ###
        ############################

        mh_res_AddRem_EI = [-Inf, -Inf, Inf, Inf]
        AddRem_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[4] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_EI, AddRem_EI_track =
            metropolis_hastings_step_aug(propose_AddRem_EI_update, AddRem_EI_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)
        end

        ############################
        ### Add/Remove Det event ###
        ############################

        mh_res_AddRem_Det = [-Inf, -Inf, Inf, Inf]
        AddRem_dets_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[5] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_Det, AddRem_dets_track =
            metropolis_hastings_step_aug(propose_AddRem_Det_update, AddRem_Det_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)
        end

        ##############################
        ### Add/Remove Death event ###
        ##############################

        mh_res_AddRem_Death = [-Inf, -Inf, Inf, Inf]
        AddRem_deaths_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[6] == true

        combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_Death, AddRem_deaths_track =
          metropolis_hastings_step_aug(propose_AddRem_Deaths_update, AddRem_Death_posterior, combi_dfs, llh_array,
                                      epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)
        end

        ########################
        ### Add/Remove p_env ###
        ########################


        mh_res_AddRem_penv = [-Inf, -Inf, Inf, Inf]
        AddRem_penv_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[7] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_penv, AddRem_penv_track =
            metropolis_hastings_step_aug(propose_AddRem_penv_update, AddRem_penv_posterior, combi_dfs, llh_array,
                                        epi_params_cur, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, scaling)
        end

        #################################
        ### Add/Remove Movement event ###
        #################################

        mh_res_AddRem_moves = [-Inf, -Inf, Inf, Inf]
        AddRem_moves_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[8] == true

          combi_dfs, llh_array, p_env_llh_array, mh_res_AddRem_moves, moves_record, AddRem_moves_track =
            metropolis_hastings_step_moves(combi_dfs, llh_array, epi_params_cur, p_env_llh_array,
                                            first_h_of_p_pos, p_idx_first_h_of_p, scaling, moves_record)
        end


        ##########################
        ### Record the results ###
        ##########################

            # Record parameters
            res[it,:] = [log_epi_params_cur; it]

            # Record other
            other_res[it,:] = [λ_inf, λ_det, m_inf, m_det, mh_res_inf[1], mh_res_det[1], mh_res_inf[2], mh_res_det[2], mh_res_inf[3], mh_res_det[3], mh_res_inf[4], mh_res_det[4], it]

            # Record data aug
            aug_res[it, :] = [mh_res_move_SE ; mh_res_move_EI ; mh_res_AddRem_SE ; mh_res_AddRem_EI ; mh_res_AddRem_Det ; mh_res_AddRem_Death ; mh_res_AddRem_penv ; mh_res_AddRem_moves]

            # Record update tracking data
            move_SE_tracker[it, :] = move_SE_track
            move_EI_tracker[it, :] = move_EI_track
            AddRem_SE_tracker[it, :] = AddRem_SE_track
            AddRem_EI_tracker[it, :] = AddRem_EI_track
            AddRem_Dets_tracker[it, :] = AddRem_dets_track
            AddRem_Deaths_tracker[it, :] = AddRem_deaths_track
            AddRem_penv_tracker[it, :] = AddRem_penv_track
            AddRem_Moves_tracker[it, :] = AddRem_moves_track


            # Update count
            if it == (n_tune*25)
              println("it = ", it)
              n_tune = n_tune + 1
            end
            it = it + 1

            # if it == 10000
            #   break
            # end

      end #end of while

  res[:, 1:9] = inverselogtransform(res[:, 1:9])

  res = DataFrame(res, :auto)
  rename!(res, [:β_c, :β_b, :γ, :F, :ϵ, :ρ, :ρ_E, :θ_bb, :θ_bd, :sample])

  other_res = DataFrame(other_res, :auto)
  rename!(other_res, [:λ_inf, :λ_det, :m_inf, :m_det,
                      :acc_inf, :acc_det,
                      :log_α_ratio_inf, :log_α_ratio_det,
                      :post_inf, :post_det,
                      :post_prime_inf, :post_prime_det,
                      :sample])

  aug_res = DataFrame(aug_res, :auto)
  rename!(aug_res, [:is_accepted_move_SE, :log_α_ratio_move_SE, :post_move_SE, :post_prime_move_SE,
                      :is_accepted_move_EI, :log_α_ratio_move_EI, :post_move_EI, :post_prime_move_EI,
                      :is_accepted_AddRem_SE, :log_α_ratio_AddRem_SE, :post_AddRem_SE, :post_prime_AddRem_SE,
                      :is_accepted_AddRem_EI, :log_α_ratio_AddRem_EI, :post_AddRem_EI, :post_prime_AddRem_EI,
                      :is_accepted_AddRem_Det, :log_α_ratio_AddRem_Det, :post_AddRem_Det, :post_prime_AddRem_Det,
                      :is_accepted_AddRem_Death, :log_α_ratio_AddRem_Death, :post_AddRem_Death, :post_prime_AddRem_Death,
                      :is_accepted_AddRem_penv, :log_α_ratio_AddRem_penv, :post_AddRem_penv, :post_prime_AddRem_penv,
                      :is_accepted_AddRem_moves, :log_α_ratio_AddRem_moves, :post_AddRem_moves, :post_prime_AddRem_moves])

  move_SE_tracker = DataFrame(move_SE_tracker, :auto)
  rename!(move_SE_tracker, [:pos_id, :t, :is_accepted, :reason,
                            :Δ_time, :num_moved])

  move_EI_tracker = DataFrame(move_EI_tracker, :auto)
  rename!(move_EI_tracker, [:pos_id, :t, :is_accepted, :reason,
                            :Δ_time, :num_moved])

  AddRem_SE_tracker = DataFrame(AddRem_SE_tracker, :auto)
  rename!(AddRem_SE_tracker, [:pos_id, :t, :is_accepted, :reason,
                              :SE_before, :SE_after, :Δ_diff, :cS, :prob])

  AddRem_EI_tracker = DataFrame(AddRem_EI_tracker, :auto)
  rename!(AddRem_EI_tracker, [:pos_id, :t, :is_accepted, :reason,
                              :EI_before, :EI_after, :Δ_diff, :cE, :prob])

  AddRem_Dets_tracker = DataFrame(AddRem_Dets_tracker, :auto)
  rename!(AddRem_Dets_tracker, [:pos_id, :t, :is_accepted, :reason,
                                :Edet_before, :Idet_before, :Edet_after, :Idet_after,
                                :ΔE_diff, :ΔI_diff, :cE, :cI])

  AddRem_Deaths_tracker = DataFrame(AddRem_Deaths_tracker, :auto)
  rename!(AddRem_Deaths_tracker, [:pos_id, :t, :is_accepted, :reason,
                                  :Sdths_before, :Edths_before, :Idths_before,
                                  :Sdths_after, :Edths_after, :Idths_after,
                                  :ΔS_diff, :ΔE_diff, :ΔI_diff, :cS, :cE, :cI])

  AddRem_penv_tracker = DataFrame(AddRem_penv_tracker, :auto)
  rename!(AddRem_penv_tracker, [:pos_id, :t, :is_accepted, :reason,
                                :r_pres_before, :n_pres_before, :r_pres_after, :n_pres_after,
                                :Δr_diff, :Δn_diff, :p_env_prev, :pI])

  AddRem_Moves_tracker = DataFrame(AddRem_Moves_tracker, :auto)
  rename!(AddRem_Moves_tracker, [:pos_id, :t, :is_accepted, :reason,
                                 :Soff_before, :Eoff_before, :Ioff_before,
                                 :Soff_after, :Eoff_after, :Ioff_after,
                                 :cS, :cE, :cI])

  return(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker,
          AddRem_Dets_tracker, AddRem_Deaths_tracker, AddRem_penv_tracker, AddRem_Moves_tracker)
end
