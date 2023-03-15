
include("VisualisationFunctions.jl")

##############################
#### 1. Badger parameters ####
##############################

r1 = CSV.read("Inference/Badger/Results_1_[Bad]/res_1_[Bad].csv", DataFrame)
or1 = CSV.read("Inference/Badger/Results_1_[Bad]/other_res_1_[Bad].csv", DataFrame)
ar1 = CSV.read("Inference/Badger/Results_1_[Bad]/aug_res_1_[Bad].csv", DataFrame)
tr1 = CSV.read("Inference/Badger/Results_1_[Bad]/tuning_res_1_[Bad].csv", DataFrame)
ut1 = CSV.read("Inference/Badger/Results_1_[Bad]/update_tracker_1_[Bad].csv", DataFrame)

number_of_samples = size(r1, 1)

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52

param_names = ["θ_bb", "θ_bd"]

describe(r1[:, :])

### SAMPLES ###

chain_full = Chains(Array(r1[:, 8:9]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r1)[1:5000, 8:9], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r1)[5001:end, 8:9], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_bad_posterior_compare(or1, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_diagnostic_likelihoods_tuning.pdf"))

post_plot2 = plot_bad_posterior_compare(or1, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_diagnostic_likelihoods_post_tuning.pdf"))

post_plot3 = plot_bad_posterior_compare(or1, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_diagnostic_likelihoods_full.pdf"))


acc_ = calc_acc_bad(or1)

acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

savefig(acc_plot, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_cumulative_acceptance_rate.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_bad_scaled(or1, tr1, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_batch_acceptance_rate_tuning.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_bad_scaled(or1, tr1, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_batch_acceptance_rate_post_tuning.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_bad(or1, tr1, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_batch_acceptance_rate_full.pdf"))


### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r1[5001:end, 8:9]), label = ["θ_bb", "θ_bd"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Badger/Results_1_[Bad]/", "BadPars", "_Multi_Corr_bbg.pdf"))

### PLOT SETTINGS ###

Plots.resetfontsizes()


println(describe(r1[:, 8:9]))
println(describe(or1))
println(describe(tr1))



###########################
#### 2. All parameters ####
###########################

r2 = CSV.read("Inference/Badger/Results_2_[Inf][Det][Bad]/res_2_[Inf][Det][Bad].csv", DataFrame)
or2 = CSV.read("Inference/Badger/Results_2_[Inf][Det][Bad]/other_res_2_[Inf][Det][Bad].csv", DataFrame)
ar2 = CSV.read("Inference/Badger/Results_2_[Inf][Det][Bad]/aug_res_2_[Inf][Det][Bad].csv", DataFrame)
tr2 = CSV.read("Inference/Badger/Results_2_[Inf][Det][Bad]/tuning_res_2_[Inf][Det][Bad].csv", DataFrame)
ut2 = CSV.read("Inference/Badger/Results_2_[Inf][Det][Bad]/update_tracker_2_[Inf][Det][Bad].csv", DataFrame)

number_of_samples = size(r2, 1)

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015
F_tr = 0.004
ϵ_tr = 0.05

ρ_tr = 0.75
ρ_E_tr = 0.2

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52

param_names = ["β_c", "β_b", "γ", "F", "ϵ", "ρ", "ρ_E", "θ_bb", "θ_bd"]

param_names_inf = ["β_c", "β_b", "γ", "F", "ϵ"]

param_names_det = ["ρ", "ρ_E"]

param_names_bad = ["θ_bb", "θ_bd"]

describe(r2[:, :])

### SAMPLES ###

chain_full = Chains(Array(r2[:, 1:9]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r2)[1:5000, 1:9], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r2)[5001:end, 1:9], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or2, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

post_plot2 = plot_inf_posterior_compare(or2, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

post_plot3 = plot_inf_posterior_compare(or2, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


post_plot4 = plot_det_posterior_compare(or2, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot4, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

post_plot5 = plot_det_posterior_compare(or2, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot5, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

post_plot6 = plot_det_posterior_compare(or2, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot6, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))


post_plot7 = plot_bad_posterior_compare(or2, 1, 5000, string("Badger Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot7, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_tuning_bad.pdf"))

post_plot8 = plot_bad_posterior_compare(or2, 5001, number_of_samples, string("Badger Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot8, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_post_tuning_bad.pdf"))

post_plot9 = plot_bad_posterior_compare(or2, 1, number_of_samples, string("Badger Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot9, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_diagnostic_likelihoods_full_bad.pdf"))



acc_inf = calc_acc_inf(or2)

acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

savefig(acc_plot1, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or2, tr2, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or2, tr2, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf(or2, tr2, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




acc_det = calc_acc_det(or2)

acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

savefig(acc_plot5, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or2, tr2, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot6, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or2, tr2, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot7, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

acc_plot8 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot8, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))





acc_bad = calc_acc_bad(or2)

acc_plot9 = plot_acc(acc_bad, string("Badger Parameters Cumulative Acceptance Rates"))

savefig(acc_plot9, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_cumulative_acceptance_rate_bad.pdf"))


acc_plot10 = plot_acc_rates_against_λ_and_m_bad_scaled(or2, tr2, 1, 5000, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot10, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_tuning_bad.pdf"))

acc_plot11 = plot_acc_rates_against_λ_and_m_bad_scaled(or2, tr2, 5001, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot11, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_post_tuning_bad.pdf"))

acc_plot12 = plot_acc_rates_against_λ_and_m_bad(or2, tr2, 1, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot12, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_batch_acceptance_rate_full_bad.pdf"))




### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r2[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_Corr_bbg.pdf"))

corr_plot2 = corrplot(Array(r2[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot2, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_Corr_bbf.pdf"))

corr_plot3 = corrplot(Array(r2[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot3, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_Corr_inf.pdf"))


corr_plot4 = corrplot(Array(r2[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot4, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_Corr_det.pdf"))


corr_plot5 = corrplot(Array(r2[5001:end, 8:9]), label = ["θ_bb", "θ_bd"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot5, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_Corr_bad.pdf"))




multi_kde_plot1 = @df r2[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot1, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_KDE_bbg.pdf"))

multi_kde_plot2 = @df r2[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot2, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_KDE_bbf.pdf"))

multi_kde_plot3 = @df r2[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot3, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_KDE_inf.pdf"))


multi_kde_plot4 = @df r2[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot4, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_KDE_det.pdf"))


multi_kde_plot5 = @df r2[5001:end, 8:9] cornerplot([:θ_bb, :θ_bd], label = ["θ_bb", "θ_bd"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot5, string("Visualisation/Badger/Results_2_[Inf][Det][Bad]/", "allpars", "_Multi_KDE_bad.pdf"))

### PLOT SETTINGS ###

Plots.resetfontsizes()


println(describe(r2[:, :]))
println(describe(or2))
println(describe(tr2))




#########################################
#### 3. All parameters and SE and EI ####
#########################################

r3 = CSV.read("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", DataFrame)
or3 = CSV.read("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/other_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", DataFrame)
ar3 = CSV.read("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/aug_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", DataFrame)
tr3 = CSV.read("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/tuning_res_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", DataFrame)
ut3 = CSV.read("Inference/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/update_tracker_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI].csv", DataFrame)

number_of_samples = size(r3, 1)

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015
F_tr = 0.004
ϵ_tr = 0.05

ρ_tr = 0.75
ρ_E_tr = 0.2

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52

param_names = ["β_c", "β_b", "γ", "F", "ϵ", "ρ", "ρ_E", "θ_bb", "θ_bd"]

param_names_inf = ["β_c", "β_b", "γ", "F", "ϵ"]

param_names_det = ["ρ", "ρ_E"]

param_names_bad = ["θ_bb", "θ_bd"]

describe(r3[:, :])

### SAMPLES ###

chain_full = Chains(Array(r3[:, 1:9]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r3)[1:5000, 1:9], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r3)[5001:end, 1:9], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or3, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

post_plot2 = plot_inf_posterior_compare(or3, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

post_plot3 = plot_inf_posterior_compare(or3, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


post_plot4 = plot_det_posterior_compare(or3, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot4, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

post_plot5 = plot_det_posterior_compare(or3, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot5, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

post_plot6 = plot_det_posterior_compare(or3, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot6, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))


post_plot7 = plot_bad_posterior_compare(or3, 1, 5000, string("Badger Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot7, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_tuning_bad.pdf"))

post_plot8 = plot_bad_posterior_compare(or3, 5001, number_of_samples, string("Badger Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot8, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_post_tuning_bad.pdf"))

post_plot9 = plot_bad_posterior_compare(or3, 1, number_of_samples, string("Badger Parameters Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot9, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_diagnostic_likelihoods_full_bad.pdf"))



acc_inf = calc_acc_inf(or3)

acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

savefig(acc_plot1, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf(or3, tr3, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




acc_det = calc_acc_det(or3)

acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

savefig(acc_plot5, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot6, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot7, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

acc_plot8 = plot_acc_rates_against_λ_and_m_det(or3, tr3, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot8, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))





acc_bad = calc_acc_bad(or3)

acc_plot9 = plot_acc(acc_bad, string("Badger Parameters Cumulative Acceptance Rates"))

savefig(acc_plot9, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_cumulative_acceptance_rate_bad.pdf"))


acc_plot10 = plot_acc_rates_against_λ_and_m_bad_scaled(or3, tr3, 1, 5000, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot10, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_tuning_bad.pdf"))

acc_plot11 = plot_acc_rates_against_λ_and_m_bad_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot11, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_post_tuning_bad.pdf"))

acc_plot12 = plot_acc_rates_against_λ_and_m_bad(or3, tr3, 1, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot12, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_batch_acceptance_rate_full_bad.pdf"))




### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r3[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_Corr_bbg.pdf"))

corr_plot2 = corrplot(Array(r3[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot2, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_Corr_bbf.pdf"))

corr_plot3 = corrplot(Array(r3[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot3, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_Corr_inf.pdf"))


corr_plot4 = corrplot(Array(r3[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot4, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_Corr_det.pdf"))


corr_plot5 = corrplot(Array(r3[5001:end, 8:9]), label = ["θ_bb", "θ_bd"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot5, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_Corr_bad.pdf"))




multi_kde_plot1 = @df r3[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot1, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_KDE_bbg.pdf"))

multi_kde_plot2 = @df r3[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot2, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_KDE_bbf.pdf"))

multi_kde_plot3 = @df r3[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot3, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_KDE_inf.pdf"))


multi_kde_plot4 = @df r3[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot4, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_KDE_det.pdf"))


multi_kde_plot5 = @df r3[5001:end, 8:9] cornerplot([:θ_bb, :θ_bd], label = ["θ_bb", "θ_bd"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot5, string("Visualisation/Badger/Results_3_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI]/", "allpars", "_Multi_KDE_bad.pdf"))

### PLOT SETTINGS ###

Plots.resetfontsizes()


println(describe(r3[:, :]))
println(describe(or3))
println(describe(tr3))




#####################################################
#### 4. All parameters and SE and EI and bb + bd ####
#####################################################

r4 = CSV.read("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", DataFrame)
        or4 = CSV.read("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/other_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", DataFrame)
        ar4 = CSV.read("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/aug_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", DataFrame)
        tr4 = CSV.read("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/tuning_res_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", DataFrame)
        ut4 = CSV.read("Inference/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/update_tracker_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth].csv", DataFrame)

        number_of_samples = size(r4, 1)

        β_c_tr = 0.002
        β_b_tr = 0.004
        γ_tr = 0.015
        F_tr = 0.004
        ϵ_tr = 0.05

        ρ_tr = 0.75
        ρ_E_tr = 0.2

        θ_bb_tr = 0.25/52
        θ_bd_tr = 0.25/52

        param_names = ["β_c", "β_b", "γ", "F", "ϵ", "ρ", "ρ_E", "θ_bb", "θ_bd"]

        param_names_inf = ["β_c", "β_b", "γ", "F", "ϵ"]

        param_names_det = ["ρ", "ρ_E"]

        param_names_bad = ["θ_bb", "θ_bd"]

        describe(r4[:, :])

        ### SAMPLES ###

        chain_full = Chains(Array(r4[:, 1:9]), param_names)
        chainplot1 = plot(chain_full)

        savefig(chainplot1, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_chain_full.pdf"))

        chain_mix1 = Chains(Array(r4)[1:5000, 1:9], param_names)
        chainplot2 = plot(chain_mix1)

        savefig(chainplot2, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_chain_tuning.pdf"))

        chain_mix2 = Chains(Array(r4)[5001:end, 1:9], param_names)
        chainplot3 = plot(chain_mix2)

        savefig(chainplot3, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_chain_post_tuning.pdf"))

        ### PLOT SETTINGS ###

        Plots.scalefontsizes(3)

        ### AUXILARY ###

        post_plot1 = plot_inf_posterior_compare(or4, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot1, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

        post_plot2 = plot_inf_posterior_compare(or4, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot2, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

        post_plot3 = plot_inf_posterior_compare(or4, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot3, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


        post_plot4 = plot_det_posterior_compare(or4, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot4, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

        post_plot5 = plot_det_posterior_compare(or4, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot5, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

        post_plot6 = plot_det_posterior_compare(or4, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot6, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))


        post_plot7 = plot_bad_posterior_compare(or4, 1, 5000, string("Badger Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot7, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_tuning_bad.pdf"))

        post_plot8 = plot_bad_posterior_compare(or4, 5001, number_of_samples, string("Badger Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot8, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_post_tuning_bad.pdf"))

        post_plot9 = plot_bad_posterior_compare(or4, 1, number_of_samples, string("Badger Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot9, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_diagnostic_likelihoods_full_bad.pdf"))



        acc_inf = calc_acc_inf(or4)

        acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot1, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


        acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot2, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

        acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot3, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

        acc_plot4 = plot_acc_rates_against_λ_and_m_inf(or4, tr4, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot4, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




        acc_det = calc_acc_det(or4)

        acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot5, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


        acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or4, tr4, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot6, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

        acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or4, tr4, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot7, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

        acc_plot8 = plot_acc_rates_against_λ_and_m_det(or4, tr4, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot8, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))





        acc_bad = calc_acc_bad(or4)

        acc_plot9 = plot_acc(acc_bad, string("Badger Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot9, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_cumulative_acceptance_rate_bad.pdf"))


        acc_plot10 = plot_acc_rates_against_λ_and_m_bad_scaled(or4, tr4, 1, 5000, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot10, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_tuning_bad.pdf"))

        acc_plot11 = plot_acc_rates_against_λ_and_m_bad_scaled(or4, tr4, 5001, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot11, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_post_tuning_bad.pdf"))

        acc_plot12 = plot_acc_rates_against_λ_and_m_bad(or4, tr4, 1, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot12, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_batch_acceptance_rate_full_bad.pdf"))




        ### CORR AND MUTLI KDE PLOTS ###

        corr_plot1 = corrplot(Array(r4[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot1, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_Corr_bbg.pdf"))

        corr_plot2 = corrplot(Array(r4[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot2, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_Corr_bbf.pdf"))

        corr_plot3 = corrplot(Array(r4[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot3, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_Corr_inf.pdf"))


        corr_plot4 = corrplot(Array(r4[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot4, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_Corr_det.pdf"))


        corr_plot5 = corrplot(Array(r4[5001:end, 8:9]), label = ["θ_bb", "θ_bd"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot5, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_Corr_bad.pdf"))




        multi_kde_plot1 = @df r4[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot1, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_KDE_bbg.pdf"))

        multi_kde_plot2 = @df r4[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot2, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_KDE_bbf.pdf"))

        multi_kde_plot3 = @df r4[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot3, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_KDE_inf.pdf"))


        multi_kde_plot4 = @df r4[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot4, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_KDE_det.pdf"))


        multi_kde_plot5 = @df r4[5001:end, 8:9] cornerplot([:θ_bb, :θ_bd], label = ["θ_bb", "θ_bd"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot5, string("Visualisation/Badger/Results_4_[Inf][Det][Bad][mbSE][mbEI][arbSE][arbEI][arbDeath][arbBirth]/", "allpars", "_Multi_KDE_bad.pdf"))

        ### PLOT SETTINGS ###

        Plots.resetfontsizes()


        println(describe(r4[:, :]))
        println(describe(or4))
        println(describe(tr4))



################
#### 5. All ####
################

r5 = CSV.read("Inference/Badger/Results_5_[All]/res_5_[All].csv", DataFrame)
        or5 = CSV.read("Inference/Badger/Results_5_[All]/other_res_5_[All].csv", DataFrame)
        ar5 = CSV.read("Inference/Badger/Results_5_[All]/aug_res_5_[All].csv", DataFrame)
        tr5 = CSV.read("Inference/Badger/Results_5_[All]/tuning_res_5_[All].csv", DataFrame)
        ut5 = CSV.read("Inference/Badger/Results_5_[All]/update_tracker_5_[All].csv", DataFrame)

        number_of_samples = size(r5, 1)

        β_c_tr = 0.002
        β_b_tr = 0.004
        γ_tr = 0.015
        F_tr = 0.004
        ϵ_tr = 0.05

        ρ_tr = 0.75
        ρ_E_tr = 0.2

        θ_bb_tr = 0.25/52
        θ_bd_tr = 0.25/52

        param_names = ["β_c", "β_b", "γ", "F", "ϵ", "ρ", "ρ_E", "θ_bb", "θ_bd"]

        param_names_inf = ["β_c", "β_b", "γ", "F", "ϵ"]

        param_names_det = ["ρ", "ρ_E"]

        param_names_bad = ["θ_bb", "θ_bd"]

        describe(r5[:, :])

        ### SAMPLES ###

        chain_full = Chains(Array(r5[:, 1:9]), param_names)
        chainplot1 = plot(chain_full)

        savefig(chainplot1, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_chain_full.pdf"))

        chain_mix1 = Chains(Array(r5)[1:5000, 1:9], param_names)
        chainplot2 = plot(chain_mix1)

        savefig(chainplot2, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_chain_tuning.pdf"))

        chain_mix2 = Chains(Array(r5)[5001:end, 1:9], param_names)
        chainplot3 = plot(chain_mix2)

        savefig(chainplot3, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_chain_post_tuning.pdf"))

        ### PLOT SETTINGS ###

        Plots.scalefontsizes(3)

        ### AUXILARY ###

        post_plot1 = plot_inf_posterior_compare(or5, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot1, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

        post_plot2 = plot_inf_posterior_compare(or5, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot2, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

        post_plot3 = plot_inf_posterior_compare(or5, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot3, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


        post_plot4 = plot_det_posterior_compare(or5, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot4, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

        post_plot5 = plot_det_posterior_compare(or5, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot5, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

        post_plot6 = plot_det_posterior_compare(or5, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot6, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))


        post_plot7 = plot_bad_posterior_compare(or5, 1, 5000, string("Badger Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

        savefig(post_plot7, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_tuning_bad.pdf"))

        post_plot8 = plot_bad_posterior_compare(or5, 5001, number_of_samples, string("Badger Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

        savefig(post_plot8, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_post_tuning_bad.pdf"))

        post_plot9 = plot_bad_posterior_compare(or5, 1, number_of_samples, string("Badger Parameters Posterior Comparison Full"), :bottomright, 0.5)

        savefig(post_plot9, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_diagnostic_likelihoods_full_bad.pdf"))



        acc_inf = calc_acc_inf(or5)

        acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot1, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


        acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or5, tr5, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot2, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

        acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or5, tr5, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot3, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

        acc_plot4 = plot_acc_rates_against_λ_and_m_inf(or5, tr5, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot4, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




        acc_det = calc_acc_det(or5)

        acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot5, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


        acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or5, tr5, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot6, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

        acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or5, tr5, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot7, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

        acc_plot8 = plot_acc_rates_against_λ_and_m_det(or5, tr5, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot8, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))





        acc_bad = calc_acc_bad(or5)

        acc_plot9 = plot_acc(acc_bad, string("Badger Parameters Cumulative Acceptance Rates"))

        savefig(acc_plot9, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_cumulative_acceptance_rate_bad.pdf"))


        acc_plot10 = plot_acc_rates_against_λ_and_m_bad_scaled(or5, tr5, 1, 5000, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

        savefig(acc_plot10, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_tuning_bad.pdf"))

        acc_plot11 = plot_acc_rates_against_λ_and_m_bad_scaled(or5, tr5, 5001, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

        savefig(acc_plot11, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_post_tuning_bad.pdf"))

        acc_plot12 = plot_acc_rates_against_λ_and_m_bad(or5, tr5, 1, number_of_samples, number_of_samples, string("Badger Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

        savefig(acc_plot12, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_batch_acceptance_rate_full_bad.pdf"))




        ### CORR AND MUTLI KDE PLOTS ###

        corr_plot1 = corrplot(Array(r5[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot1, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_Corr_bbg.pdf"))

        corr_plot2 = corrplot(Array(r5[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot2, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_Corr_bbf.pdf"))

        corr_plot3 = corrplot(Array(r5[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot3, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_Corr_inf.pdf"))


        corr_plot4 = corrplot(Array(r5[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot4, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_Corr_det.pdf"))


        corr_plot5 = corrplot(Array(r5[5001:end, 8:9]), label = ["θ_bb", "θ_bd"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(corr_plot5, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_Corr_bad.pdf"))




        multi_kde_plot1 = @df r5[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot1, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_KDE_bbg.pdf"))

        multi_kde_plot2 = @df r5[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot2, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_KDE_bbf.pdf"))

        multi_kde_plot3 = @df r5[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot3, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_KDE_inf.pdf"))


        multi_kde_plot4 = @df r5[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot4, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_KDE_det.pdf"))


        multi_kde_plot5 = @df r5[5001:end, 8:9] cornerplot([:θ_bb, :θ_bd], label = ["θ_bb", "θ_bd"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

        savefig(multi_kde_plot5, string("Visualisation/Badger/Results_5_[All]/", "allpars", "_Multi_KDE_bad.pdf"))

        ### PLOT SETTINGS ###

        Plots.resetfontsizes()


        println(describe(r5[:, :]))
        println(describe(or5))
        println(describe(tr5))
