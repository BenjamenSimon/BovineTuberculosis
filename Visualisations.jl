
include("VisualisationFunctions.jl")

#################################################
#### 1. Infection parameters - Strong Priors ####
#################################################

r1 = CSV.read("Inference/Test 1/res_infpars.csv", DataFrame)
    or1 = CSV.read("Inference/Test 1/other_res_infpars.csv", DataFrame)
    ar1 = CSV.read("Inference/Test 1/aug_res_infpars.csv", DataFrame)
    tr1 = CSV.read("Inference/Test 1/tuning_res_infpars.csv", DataFrame)
    ut1 = CSV.read("Inference/Test 1/update_tracker_infpars.csv", DataFrame)

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

    savefig(chainplot1, string("Visualisation/Test 1/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r1)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Test 1/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r1)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Test 1/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or1, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Test 1/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or1, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Test 1/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or1, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Test 1/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or1)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Test 1/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Test 1/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Test 1/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Test 1/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r1[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Test 1/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r1[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Test 1/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r1[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Test 1/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r1[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Test 1/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r1[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Test 1/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r1[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Test 1/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r1[:, 1:5]))
    println(describe(or1))
    println(describe(tr1))


#################################################
#### 2. Detection parameters - Strong Priors ####
#################################################

r2 = CSV.read("Inference/Test 2/res_detpars.csv", DataFrame)
    or2 = CSV.read("Inference/Test 2/other_res_detpars.csv", DataFrame)
    ar2 = CSV.read("Inference/Test 2/aug_res_detpars.csv", DataFrame)
    tr2 = CSV.read("Inference/Test 2/tuning_res_detpars.csv", DataFrame)
    ut2 = CSV.read("Inference/Test 2/update_tracker_detpars.csv", DataFrame)

    number_of_samples = size(r2, 1)

    ρ_tr = 0.75
    ρ_E_tr = 0.2

    param_names = ["ρ", "ρ_E"]

    describe(r2[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r2[:, 6:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Test 2/", "detpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r2)[1:5000, 6:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Test 2/", "detpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r2)[5001:end, 6:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Test 2/", "detpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_det_posterior_compare(or2, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Test 2/", "detpars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_det_posterior_compare(or2, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Test 2/", "detpars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_det_posterior_compare(or2, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Test 2/", "detpars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_det(or2)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Test 2/", "detpars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Test 2/", "detpars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Test 2/", "detpars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Test 2/", "detpars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r2[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Test 2/", "detpars", "_Multi_Corr.pdf"))


    multi_kde_plot1 = @df r2[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Test 2/", "detpars", "_Multi_KDE.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r2[:, 6:7]))
    println(describe(or2[:, 6:10]))
    println(describe(tr2))


###########################################
#### 3. All parameters - Strong Priors ####
###########################################

r3 = CSV.read("Inference/Test 3/res_allpars.csv", DataFrame)
    or3 = CSV.read("Inference/Test 3/other_res_allpars.csv", DataFrame)
    ar3 = CSV.read("Inference/Test 3/aug_res_allpars.csv", DataFrame)
    tr3 = CSV.read("Inference/Test 3/tuning_res_allpars.csv", DataFrame)
    ut3 = CSV.read("Inference/Test 3/update_tracker_allpars.csv", DataFrame)

    number_of_samples = size(r3, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    ρ_tr = 0.75
    ρ_E_tr = 0.2

    param_names = ["β_c", "β_b", "γ", "F", "ϵ", "ρ", "ρ_E"]

    param_names_inf = ["β_c", "β_b", "γ", "F", "ϵ"]

    param_names_det = ["ρ", "ρ_E"]

    describe(r3[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r3[:, 1:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Test 3/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r3)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Test 3/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r3)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Test 3/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or3, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or3, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or3, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or3, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or3, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or3, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Test 3/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or3)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Test 3/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or3)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Test 3/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Test 3/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r3[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Test 3/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r3[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Test 3/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r3[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Test 3/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r3[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Test 3/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r3[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Test 3/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r3[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Test 3/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r3[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Test 3/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r3[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Test 3/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r3[:, :]))
    println(describe(or3))
    println(describe(tr3))




############################################################################
#### 4. Infection parameters - Strong Priors - MoveSE Data Augmentation ####
############################################################################


r4 = CSV.read("Inference/Test 4/res_infpars_mSE.csv", DataFrame)
or4 = CSV.read("Inference/Test 4/other_res_infpars_mSE.csv", DataFrame)
ar4 = CSV.read("Inference/Test 4/aug_res_infpars_mSE.csv", DataFrame)
tr4 = CSV.read("Inference/Test 4/tuning_res_infpars_mSE.csv", DataFrame)
ut4 = CSV.read("Inference/Test 4/update_tracker_infpars_mSE.csv", DataFrame)

number_of_samples = size(r4, 1)

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015
F_tr = 0.004
ϵ_tr = 0.05

param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

describe(r4[:, :])
describe(or4[:, :])
describe(ar4[:, :])
describe(tr4[:, :])

### SAMPLES ###

chain_full = Chains(Array(r4[:, 1:5]), param_names)
chainplot1 = plot(chain_full)

savefig(chainplot1, string("Visualisation/Test 4/", "InfPars", "_chain_full.pdf"))

chain_mix1 = Chains(Array(r4)[1:5000, 1:5], param_names)
chainplot2 = plot(chain_mix1)

savefig(chainplot2, string("Visualisation/Test 4/", "InfPars", "_chain_tuning.pdf"))

chain_mix2 = Chains(Array(r4)[5001:end, 1:5], param_names)
chainplot3 = plot(chain_mix2)

savefig(chainplot3, string("Visualisation/Test 4/", "InfPars", "_chain_post_tuning.pdf"))

### PLOT SETTINGS ###

Plots.scalefontsizes(3)

### AUXILARY ###

post_plot1 = plot_inf_posterior_compare(or4, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

savefig(post_plot1, string("Visualisation/Test 4/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

post_plot2 = plot_inf_posterior_compare(or4, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

savefig(post_plot2, string("Visualisation/Test 4/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

post_plot3 = plot_inf_posterior_compare(or4, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

savefig(post_plot3, string("Visualisation/Test 4/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


acc_ = calc_acc_inf(or4)

acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

savefig(acc_plot, string("Visualisation/Test 4/", "InfPars", "_cumulative_acceptance_rate.pdf"))


acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

savefig(acc_plot2, string("Visualisation/Test 4/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

savefig(acc_plot3, string("Visualisation/Test 4/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

savefig(acc_plot4, string("Visualisation/Test 4/", "InfPars", "_batch_acceptance_rate_full.pdf"))


### CORR AND MUTLI KDE PLOTS ###

corr_plot1 = corrplot(Array(r4[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot1, string("Visualisation/Test 4/", "InfPars", "_Multi_Corr_bbg.pdf"))

corr_plot2 = corrplot(Array(r4[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot2, string("Visualisation/Test 4/", "InfPars", "_Multi_Corr_bbf.pdf"))

corr_plot3 = corrplot(Array(r4[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(corr_plot3, string("Visualisation/Test 4/", "InfPars", "_Multi_Corr_all.pdf"))


multi_kde_plot1 = @df r4[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot1, string("Visualisation/Test 4/", "InfPars", "_Multi_KDE_bbg.pdf"))

multi_kde_plot2 = @df r4[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot2, string("Visualisation/Test 4/", "InfPars", "_Multi_KDE_bbf.pdf"))

multi_kde_plot3 = @df r4[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

savefig(multi_kde_plot3, string("Visualisation/Test 4/", "InfPars", "_Multi_KDE_all.pdf"))


### PLOT SETTINGS ###

Plots.resetfontsizes()


println(describe(r4[:, 1:5]))
println(describe(or4))
println(describe(tr4))
