using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames
using DataFramesMeta
using LinearAlgebra
using RData


function Move_SE_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur,args...)

  # For any update to the data, I need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  post = calc_llh_h(scope, llh_array_cur)

  post_prime = calc_llh_h(scope, llh_array_prime)

  return(llh_array_prime, p_env_llh_array_cur, post, post_prime)
end

# Move_SE_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true)
# Move_SE_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, 0., "blah")


function Move_EI_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # For any update to the data;
  # we need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # if the proportion of infecteds has changed;
  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# Move_EI_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)

function AddRem_SE_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, args...)

  # For any update to the data, I need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  post = calc_llh_h(scope, llh_array_cur)

  post_prime = calc_llh_h(scope, llh_array_prime)

  return(llh_array_prime, p_env_llh_array_cur, post, post_prime)
end

# AddRem_SE_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true)

function AddRem_EI_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # For any update to the data;
  # we need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # if the proportion of infecteds has changed;
  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# AddRem_EI_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)


function AddRem_Det_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # For any update to the data;
  # we need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # if the proportion of infecteds has changed;
  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# AddRem_Det_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)

function AddRem_Death_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # For any update to the data;
  # we need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # if the proportion of infecteds has changed;
  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# AddRem_Death_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)

function AddRem_Moves_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling, record_of_movements)

  # For any update to the data;
  # we need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh_inc_moves(calc_llh_array_row_all, scope, combi_dfs_prime, llh_array_prime, epi_pars, record_of_movements)

  # if the proportion of infecteds has changed;
  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# AddRem_Moves_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish, record_of_movements)


function AddRem_penv_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_penv, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime)

  return(llh_array_prime, p_env_llh_array_prime, post, post_prime)
end

# AddRem_penv_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)




function Infection_prior(epi_params_dists, epi_params_cur)

  β_c_prior = logpdf(epi_params_dists[1], epi_params_cur[1])
  β_b_prior = logpdf(epi_params_dists[2], epi_params_cur[2])
  γ_prior = logpdf(epi_params_dists[3], epi_params_cur[3])
  F_prior = logpdf(epi_params_dists[4], epi_params_cur[4])

  log_prior = β_c_prior + β_b_prior + γ_prior + F_prior

  return(log_prior)
end

function Infection_params_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, log_prior, log_prior_dists, p_env_llh_array_cur, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # For any update to the data, I need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_infection_params, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  # we also need to calculate the probability of the parish environmental reservoir

  p_env_llh_array_prime = deepcopy(p_env_llh_array_cur)
  p_env_llh_array_prime = calc_p_env_llh_array_allp(scope, p_env_llh_array_prime, combi_dfs_prime, epi_pars, first_h_of_p_pos, p_idx_first_h_of_p, scaling)


  log_prior_prime = Infection_prior(log_prior_dists, epi_pars)

  post = calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur) + log_prior

  post_prime = calc_llh_h_and_p(scope, llh_array_prime, p_env_llh_array_prime) + log_prior_prime

  return(llh_array_prime, p_env_llh_array_prime, log_prior_prime, post, post_prime)
end

# Infection_params_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, p_env_llh_array, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish, 0., epi_params_dists)

function Detection_prior(epi_params_dists, epi_params_cur)

  ρ_prior = logpdf(epi_params_dists[6], epi_params_cur[6])
  ρ_E_prior = logpdf(epi_params_dists[7], epi_params_cur[7])

  log_prior = ρ_prior + ρ_E_prior

  return(log_prior)
end

function Detection_params_posterior(combi_dfs_prime, scope, llh_array_cur, epi_pars, log_prior, log_prior_dists, args...)

  # For any update to the data, I need to calculate the likelihood of all events within the scope of the update

  llh_array_prime = deepcopy(llh_array_cur)
  llh_array_prime = calc_llh_array_allh(calc_llh_array_row_detection_params, scope, combi_dfs_prime, llh_array_prime, epi_pars)

  log_prior_prime = Detection_prior(log_prior_dists, epi_pars)

  post = calc_llh_h(scope, llh_array_cur) + log_prior

  post_prime = calc_llh_h(scope, llh_array_prime) + log_prior_prime

  return(llh_array_prime, p_env_llh_array, log_prior_prime, post, post_prime)
end

# Detection_params_posterior(combi_dfs, [1, 360, [1,2,3], 1:13], llh_array, epi_params_true, 0., epi_params_dists)




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
