using Distributions
using InvertedIndices
using Random
using DataFrames
using DataFramesMeta
using BenchmarkTools

                          #############################################
    #~~~~~~~~~~~~~~~~~~~~~### Likelihood for one farm one time step ###~~~~~~~~~~~~~~~~~~~~~~~#
                          #############################################

res = combi_array[1]
track = combi_array[2]
pers = combi_array[3]

res[1, 1, ["cS_postM", "cE_postM", "cI_postM"]]
track[2, 3, ["c_new_exp", "c_new_inf"]]
pers[3, 4, ["c_exp_prob", "c_inf_prob"]]


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60


@benchmark begin
  llh_sum = 0.

  for t in 2:360
    for farm in 1:100

      cStates_postM_t = res[farm, t, ["cS_postM", "cE_postM", "cI_postM"]]
      cEandI_t = track[farm, t, ["c_new_exp", "c_new_inf"]]
      cEandI_probs_t = pers[farm, t, ["c_exp_prob", "c_inf_prob"]]

      c_exp_llh = exposures(States_postM = cStates_postM_t, new_EandI = cEandI_t, exp_prob = cEandI_probs_t[1])

      llh_sum = llh_sum + c_exp_llh
    end
  end
end

llh_sum


@benchmark begin
  llh_sum = 0.
  prob = 0.

  for t in 2:360
    for farm in 1:100

      d = Binomial(convert(Int64, res[farm, t, "cS_postM"]), pers[farm, t, "c_exp_prob"])

      prob = logpdf(d, track[farm, t, "c_new_exp"])

      llh_sum = llh_sum + prob
    end
  end
end

llh_sum


Threads.nthreads()

@benchmark begin
  llh_sum = Threads.Atomic{Float64}(0)

  Threads.@threads for t in 2:360
    for farm in 1:100

      d = Binomial(convert(Int64, res[farm, t, "cS_postM"]), pers[farm, t, "c_exp_prob"])

      Threads.atomic_add!(llh_sum, logpdf(d, track[farm, t, "c_new_exp"]))
    end
  end
end


function llh_func()

  llh_sum = Threads.Atomic{Float64}(0)

  Threads.@threads for t in 2:360
    for farm in 1:100

      d = Binomial(convert(Int64, res[farm, t, "cS_postM"]), pers[farm, t, "c_exp_prob"])

      Threads.atomic_add!(llh_sum, logpdf(d, track[farm, t, "c_new_exp"]))
    end
  end

  return(llh_sum[])
end

llh_func()


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
### Functions to calculate movement likelihoods ###
###################################################

function moves_MHG(;States_Moves = [90000, 5000, 5000], Moves)

  prob = log_pdf_mvhyper(States_Moves, Moves)

  return(prob)
end


#################################################################
### Functions to calculate exposure and infection likelihoods ###
#################################################################

function exposures(;States_postM, new_EandI, exp_prob)

  d = Binomial(convert(Int64, States_postM[1]), exp_prob)

  prob = logpdf(d, new_EandI[1])

  return(prob)
end

function infections(;States_postM, new_EandI, inf_prob)

  d = Binomial(convert(Int64, States_postM[2]), inf_prob)

  prob = logpdf(d, new_EandI[2])

  return(prob)
end

####################################################
### Functions to calculate Detection likelihoods ###
####################################################

function exposed_detections(;States_postEI, dets, exp_det_prob)

  d = Binomial(convert(Int64, States_postEI[2]), exp_det_prob)

  prob = logpdf(d, dets[1])

  return(prob)
end

function infectious_detections(;States_postEI, dets, inf_det_prob)

  d = Binomial(convert(Int64, States_postEI[3]), inf_det_prob)

  prob = logpdf(d, dets[2])

  return(prob)
end


########################################################
### Functions to calculate Birth & Death likelihoods ###
########################################################

function cattle_death(;States_postDet, Deaths)

  prob = log_pdf_mvhyper(States_postDet, Deaths)

  return(prob)
end

function badger_births(;States_postEI, b_births, θ_bb)

  bN = sum(States_postEI)

  rate = θ_bb * bN

  d = Poisson(rate)

  prob = logpdf(d, b_births)

  return(prob)
end


function badger_S_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(convert(Int64, States_postEI[1]), θ_bd)

  prob = logpdf(d, Deaths[1])

  return(prob)
end

function badger_E_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(convert(Int64, States_postEI[2]), θ_bd)

  prob = logpdf(d, Deaths[2])

  return(prob)
end

function badger_I_deaths(;States_postEI, Deaths, θ_bd)

  d = Binomial(convert(Int64, States_postEI[3]), θ_bd)

  prob = logpdf(d, Deaths[3])

  return(prob)
end


#########################################
### Constructing the likelihood array ###
#########################################

# llh_array = Array{Array{BigFloat, 2}, 1}(undef,100)
# for k in 1:100
#   llh_array[k] = Array{BigFloat, 2}(undef, 360, 14)
# end

# 1. m_off_llh, 2. m_on_out_llh, 3. c_exp_llh, 4. c_inf_llh, 5. exp_det_llh,
# 6. inf_det_llh, 7. c_dth_llh, 8. b_exp_llh, 9. b_inf_llh, 10. b_bths_llh,
# 11. bS_dths_llh, 12. bE_dths_llh, 13. bI_dths_llh, 14. indv_moves_off_llh

###~~~~~~~~~~~
### Likelihood of everything
### ~~~~~~~~~

function calc_llh_array_row_indv_moves_off(;t, movement_records_i)

  # Extract Data for individual movements off

  indv_moves_off_llh = 0.

  for j in 1:size(movement_records_i, 1)
    cStates_gen_j_cur = movement_records_i[j, 4:6]
    Moves_j_cur = movement_records_i[j, 7:9]

    indv_moves_off_llh += moves_MHG(States_Moves = cStates_gen_j_cur, Moves = Moves_j_cur)
  end

  return(indv_moves_off_llh)
end

function calc_llh_array_row_everything(;t, res, track, pers, epi_params, movement_records)

  # Extract Data

  cStates_Moves_t = res[t,[:cS_Moves, :cE_Moves, :cI_Moves]]
  cStates_postM_t = res[t,[:cS_postM, :cE_postM, :cI_postM]]
  cStates_postEI_t = res[t,[:cS_postEI, :cE_postEI, :cI_postEI]]
  cStates_postDet_t = res[t,[:cS_postDet, :cE_postDet, :cI_postDet]]

  bStates_init_t = res[t,[:bS_init, :bE_init, :bI_init]]
  bStates_postEI_t = res[t,[:bS_postEI, :bE_postEI, :bI_postEI]]

  Moves_off_t = track[t,[:sus_off, :exp_off, :inf_off]]
  Moves_on_out_t = track[t,[:S_on_out, :E_on_out, :I_on_out]]

  cEandI_t = track[t,[:c_new_exp, :c_new_inf]]
  cEandI_probs_t = pers[t,[:c_exp_prob, :c_inf_prob]]

  bEandI_t = track[t,[:b_new_exp, :b_new_inf]]
  bEandI_probs_t = pers[t,[:b_exp_prob, :b_inf_prob]]

  dets_t = track[t,[:E_detected, :I_detected]]
  test_occur_t = track[t,:test_occur]

  cDeaths_t = track[t,[:c_S_death, :c_E_death, :c_I_death]]

  b_birth_t = track[t, :b_birth]
  bDeaths_t = track[t,[:b_S_death, :b_E_death, :b_I_death]]

  inf_det_prob = epi_params[6]
  exp_det_prob = inf_det_prob * epi_params[7]

  θ_bb = epi_params[8]
  θ_bd = epi_params[9]


  # Calculate llh components

  m_off_llh = moves_MHG(States_Moves = cStates_Moves_t, Moves = Moves_off_t)

  m_on_out_llh = moves_MHG(States_Moves = [90000, 5000, 5000], Moves = Moves_on_out_t)

  c_exp_llh = exposures(States_postM = cStates_postM_t, new_EandI = cEandI_t, exp_prob = cEandI_probs_t[1])
  c_inf_llh = infections(States_postM = cStates_postM_t, new_EandI = cEandI_t, inf_prob = cEandI_probs_t[2])

  exp_det_llh = test_occur_t * exposed_detections(States_postEI = cStates_postEI_t, dets = dets_t, exp_det_prob = exp_det_prob)
  inf_det_llh = test_occur_t * infectious_detections(States_postEI = cStates_postEI_t, dets = dets_t, inf_det_prob = inf_det_prob)
  # where test_occur = 1 if its a test week, 0 otherwise

  c_dth_llh = cattle_death(States_postDet = cStates_postDet_t, Deaths = cDeaths_t)

  b_exp_llh = exposures(;States_postM = bStates_init_t, new_EandI = bEandI_t, exp_prob = bEandI_probs_t[1])
  b_inf_llh = infections(;States_postM = bStates_init_t, new_EandI = bEandI_t, inf_prob = bEandI_probs_t[2])

  b_bths_llh = badger_births(;States_postEI = bStates_postEI_t, b_births = b_birth_t, θ_bb = θ_bb)
  bS_dths_llh = badger_S_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)
  bE_dths_llh = badger_E_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)
  bI_dths_llh = badger_I_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)

  off_row_id = res[1, 2]
  movement_records_i = movement_records[findall((movement_records[:,1] .== t) .& (movement_records[:,2] .== off_row_id)), :]

  indv_moves_off_llh = calc_llh_array_row_indv_moves_off(t = t, movement_records_i = movement_records_i)

  return(1:14, [m_off_llh, m_on_out_llh, c_exp_llh, c_inf_llh, exp_det_llh, inf_det_llh, c_dth_llh, b_exp_llh, b_inf_llh, b_bths_llh, bS_dths_llh, bE_dths_llh, bI_dths_llh, indv_moves_off_llh])
end

# calc_llh_array_row_everything(t = 1, res = combi_dfs[1,1], track = combi_dfs[1,2], pers = combi_dfs[1,3], epi_params = epi_params_true, movement_records = record_of_movements)

function calc_llh_array_h(t_start, t_end, row_func, pos_id, llh_array_cur, combi_dfs, epi_params, movement_records)

  res = combi_dfs[pos_id, 1]
  track = combi_dfs[pos_id, 2]
  pers = combi_dfs[pos_id, 3]

  for t in t_start:t_end
    changes_range, changes_t = row_func(t = t, res = res, track = track, pers = pers, epi_params = epi_params, movement_records = movement_records)

    llh_array_cur[pos_id][t, changes_range] = changes_t
  end

  return(llh_array_cur)
end

# calc_llh_array_h(1, 5, calc_llh_array_row_everything, 1, llh_array, combi_dfs, epi_params_true, record_of_movements)

function calc_llh_array_allh(row_func, scope, combi_dfs, llh_array_cur, epi_params, movement_records)

  t_start = scope[1]
  t_end = scope[2]
  pos_ids = scope[3]

  for pos in pos_ids
    calc_llh_array_h(t_start, t_end, row_func, pos, llh_array_cur, combi_dfs, epi_params, movement_records)
  end

  return(llh_array_cur)
end

# calc_llh_array_allh(calc_llh_array_row_everything, [1, 5, [1,3]], combi_dfs, llh_array, epi_params_true, record_of_movements)


###~~~~~~~~~~~
### Likelihood of all elements except individual movements
### ~~~~~~~~~

function calc_llh_array_row_all(;t, res, track, pers, epi_params)

  # Extract Data

  cStates_Moves_t = res[t,[:cS_Moves, :cE_Moves, :cI_Moves]]
  cStates_postM_t = res[t,[:cS_postM, :cE_postM, :cI_postM]]
  cStates_postEI_t = res[t,[:cS_postEI, :cE_postEI, :cI_postEI]]
  cStates_postDet_t = res[t,[:cS_postDet, :cE_postDet, :cI_postDet]]

  bStates_init_t = res[t,[:bS_init, :bE_init, :bI_init]]
  bStates_postEI_t = res[t,[:bS_postEI, :bE_postEI, :bI_postEI]]

  Moves_off_t = track[t,[:sus_off, :exp_off, :inf_off]]
  Moves_on_out_t = track[t,[:S_on_out, :E_on_out, :I_on_out]]

  cEandI_t = track[t,[:c_new_exp, :c_new_inf]]
  cEandI_probs_t = pers[t,[:c_exp_prob, :c_inf_prob]]

  bEandI_t = track[t,[:b_new_exp, :b_new_inf]]
  bEandI_probs_t = pers[t,[:b_exp_prob, :b_inf_prob]]

  dets_t = track[t,[:E_detected, :I_detected]]
  test_occur_t = track[t,:test_occur]

  cDeaths_t = track[t,[:c_S_death, :c_E_death, :c_I_death]]

  b_birth_t = track[t, :b_birth]
  bDeaths_t = track[t,[:b_S_death, :b_E_death, :b_I_death]]

  inf_det_prob = epi_params[6]
  exp_det_prob = inf_det_prob * epi_params[7]

  θ_bb = epi_params[8]
  θ_bd = epi_params[9]


  # Calculate llh components

  m_off_llh = moves_MHG(States_Moves = cStates_Moves_t, Moves = Moves_off_t)

  m_on_out_llh = moves_MHG(States_Moves = [90000, 5000, 5000], Moves = Moves_on_out_t)

  c_exp_llh = exposures(States_postM = cStates_postM_t, new_EandI = cEandI_t, exp_prob = cEandI_probs_t[1])
  c_inf_llh = infections(States_postM = cStates_postM_t, new_EandI = cEandI_t, inf_prob = cEandI_probs_t[2])

  exp_det_llh = test_occur_t * exposed_detections(States_postEI = cStates_postEI_t, dets = dets_t, exp_det_prob = exp_det_prob)
  inf_det_llh = test_occur_t * infectious_detections(States_postEI = cStates_postEI_t, dets = dets_t, inf_det_prob = inf_det_prob)
  # where test_occur = 1 if its a test week, 0 otherwise

  c_dth_llh = cattle_death(States_postDet = cStates_postDet_t, Deaths = cDeaths_t)

  b_exp_llh = exposures(;States_postM = bStates_init_t, new_EandI = bEandI_t, exp_prob = bEandI_probs_t[1])
  b_inf_llh = infections(;States_postM = bStates_init_t, new_EandI = bEandI_t, inf_prob = bEandI_probs_t[2])

  b_bths_llh = badger_births(;States_postEI = bStates_postEI_t, b_births = b_birth_t, θ_bb = θ_bb)
  bS_dths_llh = badger_S_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)
  bE_dths_llh = badger_E_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)
  bI_dths_llh = badger_I_deaths(;States_postEI = bStates_postEI_t, Deaths = bDeaths_t, θ_bd = θ_bd)

  return(1:13, [m_off_llh, m_on_out_llh, c_exp_llh, c_inf_llh, exp_det_llh, inf_det_llh, c_dth_llh, b_exp_llh, b_inf_llh, b_bths_llh, bS_dths_llh, bE_dths_llh, bI_dths_llh])
end

# test_range, test = calc_llh_array_row(t = 1, res = combi_dfs[1,1], track = combi_dfs[1,2], pers = combi_dfs[1,3], epi_params = epi_params_true)
# llh_array[1][1, test_range] = test

function calc_llh_array_h(t_start, t_end, row_func, pos_id, llh_array_cur, combi_dfs, epi_params)

  res = combi_dfs[pos_id, 1]
  track = combi_dfs[pos_id, 2]
  pers = combi_dfs[pos_id, 3]

  for t in t_start:t_end
    changes_range, changes_t = row_func(t = t, res = res, track = track, pers = pers, epi_params = epi_params)

    llh_array_cur[pos_id][t, changes_range] = changes_t
  end

  return(llh_array_cur)
end

# calc_llh_array_h(1, 5, calc_llh_array_row, 1, llh_array, combi_dfs, epi_params_true)

function calc_llh_array_allh(row_func, scope, combi_dfs, llh_array_cur, epi_params)

  t_start = scope[1]
  t_end = scope[2]
  pos_ids = scope[3]

  for pos in pos_ids
    calc_llh_array_h(t_start, t_end, row_func, pos, llh_array_cur, combi_dfs, epi_params)
  end

  return(llh_array_cur)
end

# calc_llh_array_allh(calc_llh_array_row_all, [1, 5, [1,2]], combi_dfs, llh_array, epi_params_true)


###~~~~~~~~~~~
### Likelihood of all elements and of individual movements at one timepoint only
### ~~~~~~~~~

function calc_llh_array_allh_inc_moves(row_func, scope, combi_dfs, llh_array_cur, epi_params, movement_records)

  t_start = scope[1]
  t_end = scope[2]
  pos_ids = scope[3]

  for pos in pos_ids
    calc_llh_array_h(t_start, t_end, row_func, pos, llh_array_cur, combi_dfs, epi_params)
  end

  off_row_id = combi_dfs[pos_ids[end],1][1, 2]
  movement_records_i = movement_records[findall((movement_records[:,1] .== t_start) .& (movement_records[:,2] .== off_row_id)), :]

  indv_moves_off_llh = calc_llh_array_row_indv_moves_off(t = t_start, movement_records_i = movement_records_i)

  llh_array_cur[pos_ids[end]][t_start, 14] = indv_moves_off_llh

  return(llh_array_cur)
end

# calc_llh_array_allh_inc_moves(calc_llh_array_row_all, [1, 5, [1,2]], combi_dfs, llh_array, epi_params_true, record_of_movements)

###~~~~~~~~~~~
### Likelihoods for specific selections of elements
### ~~~~~~~~~

function calc_llh_array_row_penv(;t, res, track, pers, epi_params)

  # Extract Data

  cStates_postM_t = res[t,[:cS_postM, :cE_postM, :cI_postM]]

  bStates_init_t = res[t,[:bS_init, :bE_init, :bI_init]]

  cEandI_t = track[t,[:c_new_exp, :c_new_inf]]
  cEandI_probs_t = pers[t,[:c_exp_prob, :c_inf_prob]]

  bEandI_t = track[t,[:b_new_exp, :b_new_inf]]
  bEandI_probs_t = pers[t,[:b_exp_prob, :b_inf_prob]]

  # Calculate llh components


  c_exp_llh = exposures(States_postM = cStates_postM_t, new_EandI = cEandI_t, exp_prob = cEandI_probs_t[1])

  b_exp_llh = exposures(;States_postM = bStates_init_t, new_EandI = bEandI_t, exp_prob = bEandI_probs_t[1])

  return([3,8], [c_exp_llh, b_exp_llh])
end
# Technically only need to calculate farm level for t+1 and parish stuff for t

function calc_llh_array_row_infection_params(;t, res, track, pers, epi_params)

  # Extract Data

  cStates_postM_t = res[t,[:cS_postM, :cE_postM, :cI_postM]]

  bStates_init_t = res[t,[:bS_init, :bE_init, :bI_init]]

  cEandI_t = track[t,[:c_new_exp, :c_new_inf]]
  cEandI_probs_t = pers[t,[:c_exp_prob, :c_inf_prob]]

  bEandI_t = track[t,[:b_new_exp, :b_new_inf]]
  bEandI_probs_t = pers[t,[:b_exp_prob, :b_inf_prob]]

  # Calculate llh components

  c_exp_llh = exposures(States_postM = cStates_postM_t, new_EandI = cEandI_t, exp_prob = cEandI_probs_t[1])
  c_inf_llh = infections(States_postM = cStates_postM_t, new_EandI = cEandI_t, inf_prob = cEandI_probs_t[2])

  b_exp_llh = exposures(;States_postM = bStates_init_t, new_EandI = bEandI_t, exp_prob = bEandI_probs_t[1])
  b_inf_llh = infections(;States_postM = bStates_init_t, new_EandI = bEandI_t, inf_prob = bEandI_probs_t[2])

  return([3,4,8,9], [c_exp_llh, c_inf_llh, b_exp_llh, b_inf_llh])
end

function calc_llh_array_row_detection_params(;t, res, track, pers, epi_params)

  # Extract Data

  cStates_postEI_t = res[t,[:cS_postEI, :cE_postEI, :cI_postEI]]

  dets_t = track[t,[:E_detected, :I_detected]]
  test_occur_t = track[t,:test_occur]

  inf_det_prob = epi_params[6]
  exp_det_prob = inf_det_prob * epi_params[7]


  # Calculate llh components

  exp_det_llh = test_occur_t * exposed_detections(States_postEI = cStates_postEI_t, dets = dets_t, exp_det_prob = exp_det_prob)
  inf_det_llh = test_occur_t * infectious_detections(States_postEI = cStates_postEI_t, dets = dets_t, inf_det_prob = inf_det_prob)
  # where test_occur = 1 if its a test week, 0 otherwise

  return([5,6], [exp_det_llh, inf_det_llh])
end





################################################
### Constructing the parish likelihood array ###
################################################

# cph_to_formatted_reduced = load("Data/cph_to_formatted_chesh.rds")[ids_oi, :]
# first_h_of_p = combine(first, groupby(cph_to_formatted_reduced, [:county_idx, :parish_idx]))[:, :row_id]
#
# p_idx_first_h_of_p = convert.(Int64, combine(first, groupby(cph_to_formatted_reduced, [:county_idx, :parish_idx]))[:, :parish_idx])
#
# first_h_of_p_pos = findall(in.(ids_oi, Ref(first_h_of_p)))

# p_env_llh_array = Array{Array{BigFloat, 2}, 1}(undef,24)
# for k in 1:24
#   p_env_llh_array[k] = Array{BigFloat, 2}(undef, 360, 2)
# end

function calc_p_env_llh_array_row(;t, res, track, pers, epi_params, scaling_cp)

  # Extract Data

  pcStates_init_t = res[t,[:pcS_init, :pcE_init, :pcI_init]]
  pbStates_init_t = res[t,[:pbS_init, :pbE_init, :pbI_init]]

  pI_init_t = pcStates_init_t[3] + pbStates_init_t[3]

  remaining_pressure_t = track[t, :remaining_pressure]
  new_pressure_t = track[t, :new_pressure]

  p_env_prev_t = pers[t, :p_env_prev]

  ϵ = epi_params[5]

  # Calculate llh components

  new_pressure_llh = logpdf(Poisson(pI_init_t), new_pressure_t)

  remaining_pressure_llh = logpdf(Binomial(round(Int, scaling_cp * p_env_prev_t), ϵ), remaining_pressure_t)

  return([new_pressure_llh, remaining_pressure_llh])
end

# calc_p_env_llh_array_row(t = 1, res = combi_dfs[1,1], track = combi_dfs[1,2], pers = combi_dfs[1,3], epi_params = epi_params_true, scaling_cp = 5)

function calc_p_env_llh_array_p(;t_start, t_end, p_pos_id, p_env_llh_array_cur, res, track, pers, epi_params, scaling_cp)

  for t in t_start:t_end
    changes_t = calc_p_env_llh_array_row(t = t, res = res, track = track, pers = pers, epi_params = epi_params, scaling_cp = scaling_cp)

    p_env_llh_array_cur[p_pos_id][t, 1:2] = changes_t
  end

  return(p_env_llh_array_cur)
end

# calc_p_env_llh_array_p(t_start = 1, t_end = 100, p_pos_id = 1, p_env_llh_array_cur = p_env_llh_array,
#                       res = combi_dfs[1,1], track = combi_dfs[1,2], pers = combi_dfs[1,3], epi_params = epi_params_true, scaling_cp = 5)

function calc_p_env_llh_array_allp(scope, p_env_llh_array_cur, combi_dfs, epi_params, first_h_of_p_pos, p_idx_first_h_of_p, scaling)

  # first_h_of_p = row id of the first farm in each parish
  # p_idx_first_h_of_p = is the parish id of each parish in the data set of farms of interest
  # first_h_of_p_pos = the position of the first h in the parish in the list of farms of interest

  t_start = scope[1]
  t_end = scope[2]

  for (index, value) in enumerate(first_h_of_p_pos)

    scalingcp = scaling[1][p_idx_first_h_of_p[index]]

    res_i = combi_dfs[value, 1]
    track_i = combi_dfs[value, 2]
    pers_i = combi_dfs[value, 3]

    calc_p_env_llh_array_p(t_start = t_start, t_end = t_end, p_pos_id = index,
                            p_env_llh_array_cur = p_env_llh_array_cur,
                            res = res_i, track = track_i, pers = pers_i,
                            epi_params = epi_params, scaling_cp = scalingcp)
  end

  return(p_env_llh_array_cur)
end

# calc_p_env_llh_array_allp([1, 360], p_env_llh_array, combi_dfs, epi_params_true, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)




##########################################
### Calculating the desired likelihood ###
##########################################

# calc_llh_array_allh(calc_llh_array_row_everything, [1, 360, 1:100, 1:14], combi_dfs, llh_array, epi_params_true, record_of_movements)
# calc_p_env_llh_array_allp([1, 360, 1:100, 1:13], p_env_llh_array, combi_dfs, epi_params_true, first_h_of_p_pos, p_idx_first_h_of_p, area_of_parish)

function calc_llh_h(scope, llh_array_cur)

  t_start = scope[1]
  t_end = scope[2]
  pos_ids = scope[3]
  element_range = scope[4]

  llh = 0

  for pos in pos_ids
    llh += sum(llh_array_cur[pos][t_start:t_end, element_range])
  end

  return(llh)
end

function calc_llh_h_and_p(scope, llh_array_cur, p_env_llh_array_cur)

  llh = calc_llh_h(scope, llh_array_cur)

  t_start = scope[1]
  t_end = scope[2]
  #p_pos_ids = scope[5]

  p_env_llh = 0

  for p_pos in 1:size(p_env_llh_array_cur, 1) #p_pos_ids
    p_env_llh += sum(p_env_llh_array_cur[p_pos][t_start:t_end, :])
  end

  return(llh + p_env_llh)
end

# calc_llh_h([1, 360, [1,2,3], 1:13], llh_array)
#
# calc_llh_h_and_p([1, 360, [1,2,3], 1:13], llh_array, p_env_llh_array)
