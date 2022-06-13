using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames
using DataFramesMeta
using LinearAlgebra
using CSV
using RData

include("BlkMCMC_Aug_Set1_setup.jl")
include("Req_array_likelihood.jl")
include("Req_new_DataAug_posteriors.jl")
include("Req_new_DataAug_updater.jl")
include("Blk_Adaptive_Aug_MCMC.jl")


λs_init = [0.7, 0.2]
ms_init = [(2.38/(5^0.5)), (2.38/(2^0.5))]

### Infection parameters MCMC ###

infer_inf = [true, false, false, false]
random_start = [missing, missing, missing, missing, ϵ_tr, ρ_tr, ρ_E_tr, θ_bb_tr, θ_bd_tr]
data_aug_inf = [true, false, false, false, false, false, false, false]

Random.seed!(3)

@time begin
  res_infect, other_res_infect, aug_res_infect,
   mSE_track_inf, mEI_track_inf,
   arSE_track_inf, arEI_track_inf,
   arDets_track_inf, arDeaths_track_inf,
   arPenv_track_inf, arMoves_track_inf = Blk_Adaptive_MCMC_post_all(N_its= 100000, combi_dfs = combi_dfs,
                                            epi_params_dists = epi_params_dists,
                                            epi_params_true = random_start,
                                            infer_block = infer_inf,
                                            λs_init = λs_init, ms_init = ms_init,
                                            scaling = area_of_parish,
                                            first_h_of_p_pos = first_h_of_p_pos,
                                            p_idx_first_h_of_p = p_idx_first_h_of_p,
                                            moves_record = record_of_movements,
                                            data_aug_infer = data_aug_inf)
end



plot(res_infect[:,:β_c], title = "β_c")
  plot!([β_c_tr], seriestype = :hline)
plot(res_infect[:,:β_b], title = "β_b")
  plot!([β_b_tr], seriestype = :hline)
plot(res_infect[:,:γ], title = "γ")
  plot!([γ_tr], seriestype = :hline)
plot(res_infect[:,:F], title = "F")
  plot!([F_tr], seriestype = :hline)
plot(res_infect[:,:ϵ], title = "ϵ")
  plot!([ϵ_tr], seriestype = :hline)

plot(other_res_infect[:,:λ_inf], title = "λ_inf")
plot(other_res_infect[:,:m_inf], title = "m_inf")

acc_prop_inf = cumsum(other_res_infect[:,:acc_inf])./other_res_infect[:,:sample]
acc_prop_mh_Move_SE = cumsum(aug_res_infect[:,:is_accepted_move_SE])./other_res_infect[:,:sample]
acc_prop_mh_Move_EI = cumsum(aug_res_infect[:,:is_accepted_move_EI])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_SE = cumsum(aug_res_infect[:,:is_accepted_AddRem_SE])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_EI = cumsum(aug_res_infect[:,:is_accepted_AddRem_EI])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_Det = cumsum(aug_res_infect[:,:is_accepted_AddRem_Det])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_Death = cumsum(aug_res_infect[:,:is_accepted_AddRem_Death])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_penv = cumsum(aug_res_infect[:,:is_accepted_AddRem_penv])./other_res_infect[:,:sample]
acc_prop_mh_AddRem_Moves = cumsum(aug_res_infect[:,:is_accepted_AddRem_moves])./other_res_infect[:,:sample]

plot(acc_prop_inf, title = "acc_inf")

n_tune_acc_rate = fill(-Inf, 2000)

for n_tune in 1:2000
  batch_start = ((n_tune - 1)*25) + 1
  batch_end = (n_tune*25) - 1

  bot = batch_end - batch_start
  top = sum(other_res_infect[batch_start:batch_end, 5])

  acc_prop = top/bot

  n_tune_acc_rate[n_tune] = acc_prop
end

n_tune_accept_rate = repeat(n_tune_acc_rate, inner = 25)

plot(n_tune_accept_rate[1:5000], label = "acc rate of batch")
  plot!(other_res_infect[1:5000], :λ_inf], label = "λ_inf")

plot(acc_prop_inf[1:5000], title = "acceptance rates")
  plot!(acc_prop_mh_Move_SE[1:5000], label = "Move_SE")
  plot!(acc_prop_mh_Move_EI[1:5000], label = "Move_EI")
  plot!(acc_prop_mh_AddRem_SE[1:5000], label = "AddRem_SE")
  plot!(acc_prop_mh_AddRem_EI[1:5000], label = "AddRem_EI")
  plot!(acc_prop_mh_AddRem_Det[1:5000], label = "AddRem_Det")
  plot!(acc_prop_mh_AddRem_Death[1:5000], label = "AddRem_Death")
  plot!(acc_prop_mh_AddRem_penv[1:5000], label = "AddRem_penv")
  plot!(acc_prop_mh_AddRem_Moves[1:5000], label = "AddRem_Moves")



plot(other_res_infect[200:end, :post_inf], title = "Posteriors")
  plot!(other_res_infect[200:end, :post_prime_inf], label = "post_prime")

histogram(other_res_infect[:,:post_prime_inf], title = "β_c")

histogram(res_infect[:,:β_c], title = "β_c")
  plot!([β_c_tr], seriestype = :vline)
histogram(res_infect[:,:β_b], title = "β_b")
  plot!([β_b_tr], seriestype = :vline)
histogram(res_infect[:,:γ], title = "γ")
  plot!([γ_tr], seriestype = :vline)
histogram(res_infect[:,:F], title = "F")
  plot!([F_tr], seriestype = :vline)
histogram(res_infect[:,:ϵ], title = "ϵ")
  plot!([ϵ_tr], seriestype = :vline)

describe(res_infect[:, :])

CSV.write("Inference/Set 1/res_infect.csv", res_infect, writeheader=true)
CSV.write("Inference/Set 1/other_res_infect.csv", other_res_infect, writeheader=true)
CSV.write("Inference/Set 1/aug_res_infect.csv", aug_res_infect, writeheader=true)

CSV.write("Inference/Set 1/mSE_track_inf.csv", mSE_track_inf, writeheader=true)
CSV.write("Inference/Set 1/mEI_track_inf.csv", mEI_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arSE_track_inf.csv", arSE_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arEI_track_inf.csv", arEI_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arDets_track_inf.csv", arDets_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arDeaths_track_inf.csv", arDeaths_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arPenv_track_inf.csv", arPenv_track_inf, writeheader=true)
CSV.write("Inference/Set 1/arMoves_track_inf.csv", arMoves_track_inf, writeheader=true)


### Detection Parameters MCMC ###

infer_det = [false, true, false, false]
random_start = [β_c_tr, β_b_tr, γ_tr, F_tr, ϵ_tr, missing, missing, θ_bb_tr, θ_bd_tr]

@time begin
  res_detect, other_res_detect = Blk_Adaptive_MCMC_post_all(N_its = 100000,
                                            combi_dfs = combi_dfs,
                                            epi_params_dists = epi_params_dists,
                                            epi_params_true = random_start,
                                            infer_block = infer_det,
                                            λs_init = λs_init, ms_init = ms_init,
                                            scaling = area_of_parish,
                                            first_h_of_p_pos = first_h_of_p_pos,
                                            p_idx_first_h_of_p = p_idx_first_h_of_p)
end


plot(res_detect[:,:ρ], title = "ρ")
  plot!([ρ_tr], seriestype = :hline)
plot(res_detect[:,:ρ_E], title = "ρ_E")
  plot!([ρ_E_tr], seriestype = :hline)

plot(other_res_detect[:,:λ_det], title = "λ_det")
plot(other_res_detect[:,:m_det], title = "m_det")

acc_prop_det = cumsum(other_res_detect[:,:acc_det])./other_res_detect[:,:sample]

plot(acc_prop_det, title = "acc_det")

histogram(res_detect[:,:ρ], title = "ρ")
  plot!([ρ_tr], seriestype = :vline)
histogram(res_detect[:,:ρ_E], title = "ρ_E")
  plot!([ρ_E_tr], seriestype = :vline)

describe(res_detect[10000:100000, :])

CSV.write("Inference/Set 1/res_detect.csv", res_detect, writeheader=true)
CSV.write("Inference/Set 1/other_res_detect.csv", other_res_detect, writeheader=true)
