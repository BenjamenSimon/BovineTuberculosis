
#######################
### DATA STRUCTURES ###
#######################

function create_results_arrays(N_its, epi_params_true)

  number_of_parameters = size(epi_params_true, 1)

  res = fill(-99., N_its, number_of_parameters)
  # :β_c, :β_b, :γ, :F, :ϵ, :ρ, :ρ_E, :θ_bb, :θ_bd

  other_res = fill(-99., N_its, 9)
  # :is_accepted_inf, :log_α_ratio_inf, :post_inf, :post_prime_inf,
  # :is_accepted_det, :log_α_ratio_det, :post_det, :post_prime_det,
  # :sample

  aug_res = fill(-99., N_its, 32)
  # [is_accepted, log_α_ratio, post, post_prime] for
  # Move SE, Move EI, Add/Rem SE, Add/Rem EI
  # Add/Rem Det, Add/Rem Death, Add/Rem penv, Add/Rem Moves

  tuning_res = fill(-99., N_its, 5)
  # :λ_inf, :λ_det, :m_inf, :m_det, :sample

  update_tracker = fill(-99., N_its, 100)
  # move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker,
  # AddRem_Dets_tracker, AddRem_Deaths_tracker, AddRem_penv_tracker, AddRem_Moves_tracker

  return(res, other_res, aug_res, tuning_res, update_tracker)
end


####################
### MH FUNCTIONS ###
####################

function mult_mh_accept_ratio(post, post_prime, log_q_ratio, log_params_cur, log_params_draw)
  alpha =  (sum(log_params_draw) + post_prime) - (sum(log_params_cur) + post) + log_q_ratio

  return( min(0, alpha) )
end

function mh_accept_ratio(post, post_prime, log_q_ratio)

  log_α_ratio =  min( 0,  (post_prime) - (post) + log_q_ratio )

  return(log_α_ratio)
end


# MH function for parameters

function metropolis_hastings_step_params(N_its, res, other_res, it,
                                          proposal_func, llh_update_func, posterior_func,
                                          params_cur, combi_array_cur,
                                          llh_array_cur, p_env_llh_array_cur,
                                          f_to_p_dict,
                                          n_tune, m, λ, d, covarM)



  # Propose an update

  params_draw, log_q_ratio, mixture, λ, combi_array_prime = proposal_func(N_its, res, other_res,
                                                                          it, params_cur, combi_array_cur, f_to_p_dict,
                                                                          n_tune, m, λ, d, covarM)

  # println("Proposed Params")

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == -Inf
    return(params_cur, llh_array_cur, p_env_llh_array_cur, combi_array_cur, [false, -Inf, post_cur, -Inf, 2], mixture, λ_new)
                               # [is_accepted, log_α_ratio, post_cur, post_prime, reason]
  end


  # Update llh arrays

  llh_array_prime, p_env_llh_array_prime = llh_update_func(scope, llh_array_cur, p_env_llh_array_cur, combi_array_cur, params_draw, f_to_p_dict)

  # println("Updated llh")

  # Calculate new posteriors

  post_cur = posterior_func(llh_array_cur, p_env_llh_array_cur, scope, log_prior_dists, params_cur)

  post_prime = posterior_func(llh_array_prime, p_env_llh_array_prime, scope, log_prior_dists, params_draw)


  # println("Calculated posterior")

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post_cur, post_prime, log_q_ratio)
  is_accepted = false
  reason_ = 0

  # Accept/Reject
  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    params_cur = deepcopy(params_draw)

    llh_array_cur = deepcopy(llh_array_prime)
    p_env_llh_array_cur = deepcopy(p_env_llh_array_prime)

    combi_array_cur = deepcopy(combi_array_prime)

    is_accepted = true
    reason_ = 1
  end

  return(params_cur, llh_array_cur, p_env_llh_array_cur, combi_array_cur, [is_accepted, log_α_ratio, post_cur, post_prime, reason_], mixture, λ_new)
end


# MH function for data augmentation

function metropolis_hastings_step_aug(proposal_func, posterior_func,
                                          params_cur, combi_array_cur,
                                          llh_array_cur, p_env_llh_array_cur,
                                          f_to_p_dict, ids_to_pos_dict,
                                          movements_record, movement_dict)



  # Propose an update

  combi_array_prime, log_q_ratio, scope, dataaug_track = proposal_func(combi_array_cur, params_cur, f_to_p_dict)

  # println("Proposed update")

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == -Inf
    return(combi_array_cur, llh_array_cur, p_env_llh_array_cur, [false, -Inf, post_cur, -Inf, 2], dataaug_track)
                                                              # [is_accepted, log_α_ratio, post_cur, post_prime, reason]
  end


  # Update llh arrays

  llh_array_prime, p_env_llh_array_prime = update_llh_array_ALL_excindvmoves(scope, llh_array_cur, p_env_llh_array_cur, combi_array_cur, movements_record, params_cur, movement_dict, f_to_p_dict)

  # println("Updated llh")


  # Calculate new posteriors

  post_cur = posterior_func(llh_array_cur, p_env_llh_array_cur, scope, params_cur)

  post_prime = posterior_func(llh_array_prime, p_env_llh_array_prime, scope, params_cur)


  # println("Calculated posteriors")

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post_cur, post_prime, log_q_ratio)
  is_accepted = false
  reason_ = 0

  # Accept/Reject
  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    llh_array_cur = deepcopy(llh_array_prime)
    p_env_llh_array_cur = deepcopy(p_env_llh_array_prime)

    combi_array_cur = deepcopy(combi_array_prime)

    is_accepted = true
    reason_ = 1
  end

  return(llh_array_cur, p_env_llh_array_cur, combi_array_cur, [is_accepted, log_α_ratio, post_cur, post_prime, reason_], dataaug_track)
end


# MH function for data augmentation

function metropolis_hastings_step_moves(params_cur, combi_array_cur,
                                          llh_array_cur, p_env_llh_array_cur,
                                          f_to_p_dict, ids_to_pos_dict,
                                          movements_record, movement_dict)



  # Propose an update

  combi_array_prime, log_q_ratio, scope, movements_record_prime, dataaug_track = propose_AddRem_Movements(combi_array_cur, params_cur, movement_record, movement_dict, f_to_p_dict, ids_to_pos_dict)

  # println("Proposed update")

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == -Inf
    return(combi_array_cur, llh_array_cur, p_env_llh_array_cur, [false, -Inf, post_cur, -Inf, 2], movements_record, dataaug_track)
                                                              # [is_accepted, log_α_ratio, post_cur, post_prime, reason]
  end


  # Update llh arrays

  llh_array_prime, p_env_llh_array_prime = update_llh_array_ALL(scope, llh_array_cur, p_env_llh_array_cur, combi_array_prime, movements_record_prime, params_cur, movement_dict, f_to_p_dict)

  # println("Updated llh")


  # Calculate new posteriors

  post_cur = generic_posterior_dataaug(llh_array_cur, p_env_llh_array_cur, scope, params_cur)

  post_prime = generic_posterior_dataaug(llh_array_prime, p_env_llh_array_prime, scope, params_cur)


  # println("Calculated posteriors")

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post_cur, post_prime, log_q_ratio)
  is_accepted = false
  reason_ = 0

  # Accept/Reject
  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    llh_array_cur = deepcopy(llh_array_prime)
    p_env_llh_array_cur = deepcopy(p_env_llh_array_prime)

    combi_array_cur = deepcopy(combi_array_prime)
    movements_record = deepcopy(movements_record_prime)

    is_accepted = true
    reason_ = 1
  end

  return(llh_array_cur, p_env_llh_array_cur, combi_array_cur, [is_accepted, log_α_ratio, post_cur, post_prime, reason_], movements_record, dataaug_track)
end


#####################
### Adaptive MCMC ###
#####################


function Blk_Adaptive_RWM_MCMC(;N_its, combi_array, moves_record,
                                params_init, tuning,
                                dict_of_movements, f_to_p_dict)

  ##################
  ### Extraction ###
  ##################

    λ_inf = tuning[1]
    λ_det = tuning[2]

    m_inf = tuning[3]
    m_det = tuning[4]

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

    params_cur = Initialise(params_init, ...)

    covarM_inf = CovMatrix()
    covarM_det = CovMatrix()

  ########################
  ### Results matrices ###
  ########################

    res, other_res, aug_res, tuning_res, update_tracker = create_results_arrays(N_its, param_names)

  ##########################
  ### Functional objects ###
  ##########################

    it = 1

    Δm_inf = m_init/100
    Δm_det = m_init/100

    n_tune = 1


  ##################
  ### Likelihood ###
  ##################

  llh_array_init = zeros(size(combi_array[1], 1), 360, 13)

  # 1. m_off_llh, 2. m_on_out_llh, 3. c_exp_llh, 4. c_inf_llh, 5. exp_det_llh,
  # 6. inf_det_llh, 7. c_dth_llh, 8. b_exp_llh, 9. b_inf_llh, 10. b_bths_llh,
  # 11. bS_dths_llh, 12. bE_dths_llh, 13. bI_dths_llh

  p_env_llh_array_init = zeros(size(combi_array[4], 1), 360, 2)

  llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL([1, 360, 1:size(combi_array[1], 1), 1:13],
                                                            llh_array_init, p_env_llh_array_init,
                                                            combi_array, record_of_movements,
                                                            params_init, dict_of_movements, f_to_p_dict)


  ###########################
  ### ~~ THE ALGORITHM ~~ ###
  ###########################

    progbar = Progress(N_its, dt=10,
                        desc = "Running $N_its iterations for $region:",
                        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                        barlen=10)

    while it <= N_its

      # println(" ", "it = ", it)

      #######################################
      ### MH step for Epidemic parameters ###
      #######################################

        mh_res_inf = [-Inf, -Inf, Inf, Inf, Inf]
                     # is_accepted, log_α_ratio, post_cur, post_prime, reason_

        if infer_block[1] == true

          if it > 1
            fit!(covarM_inf, res[(it-1), [1,2,3,4,5]])
          end

          params_cur, llh_array_cur, p_env_llh_array_cur, combi_array_cur, mh_res_inf, mixture, λ_inf =
                    metropolis_hastings_step_params(N_its, res, other_res, it,
                                                    propose_epidemic_params, update_llh_array_EPIDEMIC, Infection_params_posterior,
                                                    params_cur, combi_array_cur,
                                                    llh_array_cur, p_env_llh_array_cur,
                                                    n_tune, m_inf, λ_inf, 5, covarM_inf)


          if mixture >= 0.05 # mixture is set to 0 while tuning λ
            if mh_res_inf[1] == 0
              m_inf = m_inf - (Δm_inf/((it)^0.5))
            else
              m_inf = m_inf + 2.3*(Δm_inf/((it)^0.5))
            end
          end

        end # if infer_block[1] == true


      ########################################
      ### MH step for Detection parameters ###
      ########################################

        mh_res_det = [-Inf, -Inf, Inf, Inf, Inf]
                     # is_accepted, log_α_ratio, post_cur, post_prime, reason_

        if infer_block[2] == true

          if it > 1
            fit!(covarM_det, res[(it-1), [6,7]])
          end

          params_cur, llh_array_cur, p_env_llh_array_cur, combi_array_cur, mh_res_det, mixture, λ_det =
                    metropolis_hastings_step_params(N_its, res, other_res, it,
                                                    propose_detection_params, update_llh_array_DETECTION, Detection_params_posterior,
                                                    params_cur, combi_array_cur,
                                                    llh_array_cur, p_env_llh_array_cur,
                                                    n_tune, m_det, λ_det, 2, covarM_det)


          if mixture >= 0.05 # mixture is set to 0 while tuning λ
            if mh_res_det[1] == 0
              m_det = m_det - (Δm_det/((it)^0.5))
            else
              m_det = m_det + 2.3*(Δm_det/((it)^0.5))
            end
          end

        end # if infer_block[2] == true


      ###################################
      ### Move S→E event through time ###
      ###################################

        mh_res_move_SE = [-Inf, -Inf, Inf, Inf]
        move_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[1] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_move_SE, move_SE_track =
                                                          metropolis_hastings_step_aug(propose_Move_SE, generic_posterior_SE_dataaug,
                                                                                        params_cur, combi_array_cur,
                                                                                        llh_array_cur, p_env_llh_array_cur,
                                                                                        f_to_p_dict, ids_to_pos_dict,
                                                                                        movements_record, movement_dict)
        end


      ###################################
      ### Move E→I event through time ###
      ###################################

        mh_res_move_EI = [-Inf, -Inf, Inf, Inf]
        move_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[2] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_move_EI, move_EI_track =
                                                        metropolis_hastings_step_aug(propose_Move_EI, generic_posterior_dataaug,
                                                                                      params_cur, combi_array_cur,
                                                                                      llh_array_cur, p_env_llh_array_cur,
                                                                                      f_to_p_dict, ids_to_pos_dict,
                                                                                      movements_record, movement_dict)
        end


      ############################
      ### Add/Remove S→E event ###
      ############################

        mh_res_AddRem_SE = [-Inf, -Inf, Inf, Inf]
        AddRem_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[3] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_SE, AddRem_SE_track =
                                                  metropolis_hastings_step_aug(propose_AddRem_SE, generic_posterior_SE_dataaug,
                                                                                params_cur, combi_array_cur,
                                                                                llh_array_cur, p_env_llh_array_cur,
                                                                                f_to_p_dict, ids_to_pos_dict,
                                                                                movements_record, movement_dict)
        end


      ############################
      ### Add/Remove E→I event ###
      ############################

        mh_res_AddRem_EI = [-Inf, -Inf, Inf, Inf]
        AddRem_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[4] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_EI, AddRem_EI_track =
                                                    metropolis_hastings_step_aug(propose_AddRem_EI, generic_posterior_dataaug,
                                                                                  params_cur, combi_array_cur,
                                                                                  llh_array_cur, p_env_llh_array_cur,
                                                                                  f_to_p_dict, ids_to_pos_dict,
                                                                                  movements_record, movement_dict)
        end


      ############################
      ### Add/Remove Det event ###
      ############################

        mh_res_AddRem_Det = [-Inf, -Inf, Inf, Inf]
        AddRem_dets_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[5] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_Det, AddRem_dets_track =
                                                  metropolis_hastings_step_aug(propose_AddRem_Det, generic_posterior_dataaug,
                                                                                params_cur, combi_array_cur,
                                                                                llh_array_cur, p_env_llh_array_cur,
                                                                                f_to_p_dict, ids_to_pos_dict,
                                                                                movements_record, movement_dict)
        end


      ##############################
      ### Add/Remove Death event ###
      ##############################

        mh_res_AddRem_Death = [-Inf, -Inf, Inf, Inf]
        AddRem_deaths_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[6] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_Death, AddRem_deaths_track =
                                                metropolis_hastings_step_aug(propose_AddRem_Deaths, generic_posterior_dataaug,
                                                                              params_cur, combi_array_cur,
                                                                              llh_array_cur, p_env_llh_array_cur,
                                                                              f_to_p_dict, ids_to_pos_dict,
                                                                              movements_record, movement_dict)
        end


      ########################
      ### Add/Remove p_env ###
      ########################

        mh_res_AddRem_penv = [-Inf, -Inf, Inf, Inf]
        AddRem_penv_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[7] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_penv, AddRem_penv_trackk =
                                              metropolis_hastings_step_aug(propose_AddRem_penv, generic_posterior_dataaug,
                                                                            params_cur, combi_array_cur,
                                                                            llh_array_cur, p_env_llh_array_cur,
                                                                            f_to_p_dict, ids_to_pos_dict,
                                                                            movements_record, movement_dict)
        end


      #################################
      ### Add/Remove Movement event ###
      #################################

        mh_res_AddRem_moves = [-Inf, -Inf, Inf, Inf]
        AddRem_moves_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[8] == true

          combi_array_cur, llh_array_cur, p_env_llh_array_cur, mh_res_AddRem_moves, moves_record, AddRem_moves_track =
                                                  metropolis_hastings_step_moves(params_cur, combi_array_cur,
                                                                                  llh_array_cur, p_env_llh_array_cur,
                                                                                  f_to_p_dict, ids_to_pos_dict,
                                                                                  moves_record, movement_dict)
        end



      ##########################
      ### Record the results ###
      ##########################

        # Record parameters
        res[it,:] = params_cur

        # Record other
        other_res[it,:] =

        # Record aug
        aug_res[it,:] = # NOTE: WRITE FUNCTION TO UPDATE

        # Record tuning params
        tuning_res[it, :] = [λ_inf, m_inf, λ_det, m_det, it]

        # Record tracking
        updater_tracker[it,:] = # NOTE: WRITE FUNCTION TO UPDATE

        # if rem(it, 1000) == 0
        #   println("it = ", it)
        # end
        next!(progbar)

        # Update count
        if it == (n_tune*25)
          n_tune = n_tune + 1
        end
        it = it + 1

    end #end of while

    res, other_res, aug_res, tuning_res, update_tracker = rename_results_arrays(res, other_res, aug_res, tuning_res, update_tracker)

  return(res, other_res, aug_res, tuning_res, update_tracker)
end
