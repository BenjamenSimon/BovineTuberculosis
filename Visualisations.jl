
include("VisualisationFunctions.jl")

#################################
#### 1. Infection parameters ####
#################################

r1 = CSV.read("Inference/Test/res_infpars.csv", DataFrame)
or1 = CSV.read("Inference/Test/other_res_infpars.csv", DataFrame)
ar1 = CSV.read("Inference/Test/aug_res_infpars.csv", DataFrame)
tr1 = CSV.read("Inference/Test/tuning_res_infpars.csv", DataFrame)
ut1 = CSV.read("Inference/Test/update_tracker_infpars.csv", DataFrame)

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

savefig(chainplot1, string("Visualisation/Test/", "InfPars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r1)[1:5000, 1:5], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Test/", "InfPars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r1)[5001:end, 1:5], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Test/", "InfPars", "_chain_post_tuning.pdf"))

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or1, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Test/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

post_plot2 = plot_inf_posterior_compare(or1, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Test/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

post_plot3 = plot_inf_posterior_compare(or1, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Test/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


acc_ = calc_acc_inf(or1)

acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

savefig(acc_plot, string("Visualisation/Test/", "InfPars", "_cumulative_acceptance_rate.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf(or1, tr1, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Test/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf(or1, tr1, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Test/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf(or1, tr1, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Test/", "InfPars", "_batch_acceptance_rate_full.pdf"))

println(describe(r1[:, 1:5]))
println(describe(or1))
println(describe(tr1))


#################################################
#### 2. Infection parameters - Strong Priors ####
#################################################

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

println(describe(r2[:, 1:5]))
println(describe(or2))
println(describe(tr2))

Plots.resetfontsizes()
