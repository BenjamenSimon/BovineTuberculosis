
include("VisualisationFunctions.jl")

#################################################
#### 1. Infection parameters - Strong Priors ####
#################################################

r1 = CSV.read("Inference/Test1/res_infpars.csv", DataFrame)
or1 = CSV.read("Inference/Test1/other_res_infpars.csv", DataFrame)
ar1 = CSV.read("Inference/Test1/aug_res_infpars.csv", DataFrame)
tr1 = CSV.read("Inference/Test1/tuning_res_infpars.csv", DataFrame)
ut1 = CSV.read("Inference/Test1/update_tracker_infpars.csv", DataFrame)

number_of_samples = size(r1, 1)

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015
F_tr = 0.004
ϵ_tr = 0.05

param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

describe(r1[:, :])

### SAMPLES ###

chain_full = Chains(Array(r1[:, 1:5]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Test1/", "InfPars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r1)[1:5000, 1:5], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Test1/", "InfPars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r1)[5001:end, 1:5], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Test1/", "InfPars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or1, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Test1/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

post_plot2 = plot_inf_posterior_compare(or1, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Test1/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

post_plot3 = plot_inf_posterior_compare(or1, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Test1/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


acc_ = calc_acc_inf(or1)

acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

savefig(acc_plot, string("Visualisation/Test1/", "InfPars", "_cumulative_acceptance_rate.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Test1/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Test1/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Test1/", "InfPars", "_batch_acceptance_rate_full.pdf"))


### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r1[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Test1/", "InfPars", "_Multi_Corr_bbg.pdf"))

corr_plot2 = corrplot(Array(r1[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot2, string("Visualisation/Test1/", "InfPars", "_Multi_Corr_bbf.pdf"))

corr_plot3 = corrplot(Array(r1[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot3, string("Visualisation/Test1/", "InfPars", "_Multi_Corr_all.pdf"))


multi_kde_plot1 = @df r1[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot1, string("Visualisation/Test1/", "InfPars", "_Multi_KDE_bbg.pdf"))

multi_kde_plot2 = @df r1[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot2, string("Visualisation/Test1/", "InfPars", "_Multi_KDE_bbf.pdf"))

multi_kde_plot3 = @df r1[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot3, string("Visualisation/Test1/", "InfPars", "_Multi_KDE_all.pdf"))


### PLOT SETTINGS ###

Plots.resetfontsizes()


println(describe(r1[:, 1:5]))
println(describe(or1))
println(describe(tr1))




#################################
#### 2. Infection parameters ####
#################################

r2 = CSV.read("Inference/Test2/res_infpars.csv", DataFrame)
or2 = CSV.read("Inference/Test2/other_res_infpars.csv", DataFrame)
ar2 = CSV.read("Inference/Test2/aug_res_infpars.csv", DataFrame)
tr2 = CSV.read("Inference/Test2/tuning_res_infpars.csv", DataFrame)
ut2 = CSV.read("Inference/Test2/update_tracker_infpars.csv", DataFrame)

number_of_samples = size(r2, 1)

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015
F_tr = 0.004
ϵ_tr = 0.05

param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

describe(r2[:, :])

### SAMPLES ###

chain_full = Chains(Array(r2[:, 1:5]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Test2/", "InfPars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r2)[1:5000, 1:5], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Test2/", "InfPars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r2)[5001:end, 1:5], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Test2/", "InfPars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or2, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Test2/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

post_plot2 = plot_inf_posterior_compare(or2, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Test2/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

post_plot3 = plot_inf_posterior_compare(or2, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Test2/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


acc_ = calc_acc_inf(or2)

acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

savefig(acc_plot, string("Visualisation/Test2/", "InfPars", "_cumulative_acceptance_rate.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or2, tr2, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Test2/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or2, tr2, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Test2/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or2, tr2, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Test2/", "InfPars", "_batch_acceptance_rate_full.pdf"))


### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r2[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Test2/", "InfPars", "_Multi_Corr_bbg.pdf"))

corr_plot2 = corrplot(Array(r2[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot2, string("Visualisation/Test2/", "InfPars", "_Multi_Corr_bbf.pdf"))

corr_plot3 = corrplot(Array(r2[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot3, string("Visualisation/Test2/", "InfPars", "_Multi_Corr_all.pdf"))


multi_kde_plot1 = @df r2[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot1, string("Visualisation/Test2/", "InfPars", "_Multi_KDE_bbg.pdf"))

multi_kde_plot2 = @df r2[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot2, string("Visualisation/Test2/", "InfPars", "_Multi_KDE_bbf.pdf"))

multi_kde_plot3 = @df r2[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot3, string("Visualisation/Test2/", "InfPars", "_Multi_KDE_all.pdf"))


println(describe(r2[:, 1:5]))
println(describe(or2))
println(describe(tr2))

### PLOT SETTINGS ###

Plots.resetfontsizes()
