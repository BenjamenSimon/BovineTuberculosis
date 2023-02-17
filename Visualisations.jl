
include("VisualisationFunctions.jl")

#################################################
#### 1. Infection parameters - Strong Priors ####
#################################################

r1 = CSV.read("Inference/Results_1_[Inf]/res_1_[Inf].csv", DataFrame)
    or1 = CSV.read("Inference/Results_1_[Inf]/other_res_1_[Inf].csv", DataFrame)
    ar1 = CSV.read("Inference/Results_1_[Inf]/aug_res_1_[Inf].csv", DataFrame)
    tr1 = CSV.read("Inference/Results_1_[Inf]/tuning_res_1_[Inf].csv", DataFrame)
    ut1 = CSV.read("Inference/Results_1_[Inf]/update_tracker_1_[Inf].csv", DataFrame)

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

    savefig(chainplot1, string("Visualisation/Results_1_[Inf]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r1)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_1_[Inf]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r1)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_1_[Inf]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or1, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_1_[Inf]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or1, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_1_[Inf]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or1, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_1_[Inf]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or1)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_1_[Inf]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_1_[Inf]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_1_[Inf]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or1, tr1, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_1_[Inf]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r1[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r1[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r1[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r1[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r1[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r1[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_1_[Inf]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r1[:, 1:5]))
    println(describe(or1))
    println(describe(tr1))


#################################################
#### 2. Detection parameters - Strong Priors ####
#################################################

r2 = CSV.read("Inference/Results_2_[Det]/res_2_[Det].csv", DataFrame)
    or2 = CSV.read("Inference/Results_2_[Det]/other_res_2_[Det].csv", DataFrame)
    ar2 = CSV.read("Inference/Results_2_[Det]/aug_res_2_[Det].csv", DataFrame)
    tr2 = CSV.read("Inference/Results_2_[Det]/tuning_res_2_[Det].csv", DataFrame)
    ut2 = CSV.read("Inference/Results_2_[Det]/update_tracker_2_[Det].csv", DataFrame)

    number_of_samples = size(r2, 1)

    ρ_tr = 0.75
    ρ_E_tr = 0.2

    param_names = ["ρ", "ρ_E"]

    describe(r2[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r2[:, 6:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_2_[Det]/", "detpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r2)[1:5000, 6:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_2_[Det]/", "detpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r2)[5001:end, 6:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_2_[Det]/", "detpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_det_posterior_compare(or2, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_2_[Det]/", "detpars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_det_posterior_compare(or2, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_2_[Det]/", "detpars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_det_posterior_compare(or2, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_2_[Det]/", "detpars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_det(or2)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_2_[Det]/", "detpars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_2_[Det]/", "detpars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_2_[Det]/", "detpars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_det(or2, tr2, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_2_[Det]/", "detpars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r2[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_2_[Det]/", "detpars", "_Multi_Corr.pdf"))


    multi_kde_plot1 = @df r2[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_2_[Det]/", "detpars", "_Multi_KDE.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r2[:, 6:7]))
    println(describe(or2[:, 6:10]))
    println(describe(tr2))


###########################################
#### 3. All parameters - Strong Priors ####
###########################################

r3 = CSV.read("Inference/Results_3_[Inf][Det]/res_3_[Inf][Det].csv", DataFrame)
    or3 = CSV.read("Inference/Results_3_[Inf][Det]/other_res_3_[Inf][Det].csv", DataFrame)
    ar3 = CSV.read("Inference/Results_3_[Inf][Det]/aug_res_3_[Inf][Det].csv", DataFrame)
    tr3 = CSV.read("Inference/Results_3_[Inf][Det]/tuning_res_3_[Inf][Det].csv", DataFrame)
    ut3 = CSV.read("Inference/Results_3_[Inf][Det]/update_tracker_3_[Inf][Det].csv", DataFrame)

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

    savefig(chainplot1, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r3)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r3)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or3, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or3, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or3, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or3, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or3, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or3, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or3)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or3, tr3, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or3)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or3, tr3, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r3[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r3[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r3[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r3[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r3[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r3[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r3[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r3[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Results_3_[Inf][Det]/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r3[:, :]))
    println(describe(or3))
    println(describe(tr3))




############################################################################
#### 4. Infection parameters - Strong Priors - MoveSE Data Augmentation ####
############################################################################


r4 = CSV.read("Inference/Results_4_[Inf][mSE]/res_4_[Inf][mSE].csv", DataFrame)
    or4 = CSV.read("Inference/Results_4_[Inf][mSE]/other_res_4_[Inf][mSE].csv", DataFrame)
    ar4 = CSV.read("Inference/Results_4_[Inf][mSE]/aug_res_4_[Inf][mSE].csv", DataFrame)
    tr4 = CSV.read("Inference/Results_4_[Inf][mSE]/tuning_res_4_[Inf][mSE].csv", DataFrame)
    ut4 = CSV.read("Inference/Results_4_[Inf][mSE]/update_tracker_4_[Inf][mSE].csv", DataFrame)

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

    savefig(chainplot1, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r4)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r4)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or4, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or4, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or4, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or4)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or4, tr4, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r4[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r4[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r4[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r4[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r4[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r4[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_4_[Inf][mSE]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r4[:, 1:5]))
    println(describe(or4))
    println(describe(tr4))



############################################################################
#### 5. Infection parameters - Strong Priors - MoveEI Data Augmentation ####
############################################################################


r5 = CSV.read("Inference/Results_5_[Inf][mEI]/res_5_[Inf][mEI].csv", DataFrame)
    or5 = CSV.read("Inference/Results_5_[Inf][mEI]/other_res_5_[Inf][mEI].csv", DataFrame)
    ar5 = CSV.read("Inference/Results_5_[Inf][mEI]/aug_res_5_[Inf][mEI].csv", DataFrame)
    tr5 = CSV.read("Inference/Results_5_[Inf][mEI]/tuning_res_5_[Inf][mEI].csv", DataFrame)
    ut5 = CSV.read("Inference/Results_5_[Inf][mEI]/update_tracker_5_[Inf][mEI].csv", DataFrame)

    number_of_samples = size(r5, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r5[:, :])
    describe(or5[:, :])
    describe(ar5[:, :])
    describe(tr5[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r5[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r5)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r5)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or5, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or5, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or5, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or5)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or5, tr5, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or5, tr5, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot5 = plot_acc_rates_against_λ_and_m_inf_scaled(or5, tr5, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot5, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r5[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r5[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r5[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r5[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r5[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r5[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_5_[Inf][mEI]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r5[:, 1:5]))
    println(describe(or5))
    println(describe(tr5))



###################################################################################
#### 6. Infection parameters - Strong Priors - Add/Remove SE Data Augmentation ####
###################################################################################


r6 = CSV.read("Inference/Results_6_[Inf][arSE]/res_6_[Inf][arSE].csv", DataFrame)
    or6 = CSV.read("Inference/Results_6_[Inf][arSE]/other_res_6_[Inf][arSE].csv", DataFrame)
    ar6 = CSV.read("Inference/Results_6_[Inf][arSE]/aug_res_6_[Inf][arSE].csv", DataFrame)
    tr6 = CSV.read("Inference/Results_6_[Inf][arSE]/tuning_res_6_[Inf][arSE].csv", DataFrame)
    ut6 = CSV.read("Inference/Results_6_[Inf][arSE]/update_tracker_6_[Inf][arSE].csv", DataFrame)

    number_of_samples = size(r6, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r6[:, :])
    describe(or6[:, :])
    describe(ar6[:, :])
    describe(tr6[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r6[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r6)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r6)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or6, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or6, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or6, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or6)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or6, tr6, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or6, tr6, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot6 = plot_acc_rates_against_λ_and_m_inf_scaled(or6, tr6, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot6, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r6[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r6[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r6[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r6[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r6[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r6[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_6_[Inf][arSE]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r6[:, 1:5]))
    println(describe(or6))
    println(describe(tr6))




###################################################################################
#### 7. Infection parameters - Strong Priors - Add/Remove EI Data Augmentation ####
###################################################################################


r7 = CSV.read("Inference/Results_7_[Inf][arEI]/res_7_[Inf][arEI].csv", DataFrame)
    or7 = CSV.read("Inference/Results_7_[Inf][arEI]/other_res_7_[Inf][arEI].csv", DataFrame)
    ar7 = CSV.read("Inference/Results_7_[Inf][arEI]/aug_res_7_[Inf][arEI].csv", DataFrame)
    tr7 = CSV.read("Inference/Results_7_[Inf][arEI]/tuning_res_7_[Inf][arEI].csv", DataFrame)
    ut7 = CSV.read("Inference/Results_7_[Inf][arEI]/update_tracker_7_[Inf][arEI].csv", DataFrame)

    number_of_samples = size(r7, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r7[:, :])
    describe(or7[:, :])
    describe(ar7[:, :])
    describe(tr7[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r7[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r7)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r7)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or7, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or7, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or7, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or7)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or7, tr7, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or7, tr7, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_inf_scaled(or7, tr7, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot7, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r7[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r7[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r7[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r7[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r7[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r7[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_7_[Inf][arEI]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r7[:, 1:5]))
    println(describe(or7))
    println(describe(tr7))





#############################################################################
#### 8. Infection parameters - Strong Priors - mSE mEI Data Augmentation ####
#############################################################################


r8 = CSV.read("Inference/Results_8_[Inf][mSE][mEI]/res_8_[Inf][mSE][mEI].csv", DataFrame)
    or8 = CSV.read("Inference/Results_8_[Inf][mSE][mEI]/other_res_8_[Inf][mSE][mEI].csv", DataFrame)
    ar8 = CSV.read("Inference/Results_8_[Inf][mSE][mEI]/aug_res_8_[Inf][mSE][mEI].csv", DataFrame)
    tr8 = CSV.read("Inference/Results_8_[Inf][mSE][mEI]/tuning_res_8_[Inf][mSE][mEI].csv", DataFrame)
    ut8 = CSV.read("Inference/Results_8_[Inf][mSE][mEI]/update_tracker_8_[Inf][mSE][mEI].csv", DataFrame)

    number_of_samples = size(r8, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r8[:, :])
    describe(or8[:, :])
    describe(ar8[:, :])
    describe(tr8[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r8[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r8)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r8)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or8, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or8, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or8, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or8)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or8, tr8, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or8, tr8, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_inf_scaled(or8, tr8, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r8[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r8[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r8[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r8[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r8[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r8[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_8_[Inf][mSE][mEI]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r8[:, 1:5]))
    println(describe(or8))
    println(describe(tr8))




###############################################################################
#### 9. Infection parameters - Strong Priors - arSE arEI Data Augmentation ####
###############################################################################


r9 = CSV.read("Inference/Results_9_[Inf][arSE][arEI]/res_9_[Inf][arSE][arEI].csv", DataFrame)
    or9 = CSV.read("Inference/Results_9_[Inf][arSE][arEI]/other_res_9_[Inf][arSE][arEI].csv", DataFrame)
    ar9 = CSV.read("Inference/Results_9_[Inf][arSE][arEI]/aug_res_9_[Inf][arSE][arEI].csv", DataFrame)
    tr9 = CSV.read("Inference/Results_9_[Inf][arSE][arEI]/tuning_res_9_[Inf][arSE][arEI].csv", DataFrame)
    ut9 = CSV.read("Inference/Results_9_[Inf][arSE][arEI]/update_tracker_9_[Inf][arSE][arEI].csv", DataFrame)

    number_of_samples = size(r9, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r9[:, :])
    describe(or9[:, :])
    describe(ar9[:, :])
    describe(tr9[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r9[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r9)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r9)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or9, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or9, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or9, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or9)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or9, tr9, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or9, tr9, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot9 = plot_acc_rates_against_λ_and_m_inf_scaled(or9, tr9, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot9, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r9[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r9[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r9[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r9[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r9[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r9[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_9_[Inf][arSE][arEI]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r9[:, 1:5]))
    println(describe(or9))
    println(describe(tr9))




###############################################################################
#### 10. Infection parameters - Strong Priors - mSE arSE Data Augmentation ####
###############################################################################


r10 = CSV.read("Inference/Results_10_[Inf][mSE][arSE]/res_10_[Inf][mSE][arSE].csv", DataFrame)
    or10 = CSV.read("Inference/Results_10_[Inf][mSE][arSE]/other_res_10_[Inf][mSE][arSE].csv", DataFrame)
    ar10 = CSV.read("Inference/Results_10_[Inf][mSE][arSE]/aug_res_10_[Inf][mSE][arSE].csv", DataFrame)
    tr10 = CSV.read("Inference/Results_10_[Inf][mSE][arSE]/tuning_res_10_[Inf][mSE][arSE].csv", DataFrame)
    ut10 = CSV.read("Inference/Results_10_[Inf][mSE][arSE]/update_tracker_10_[Inf][mSE][arSE].csv", DataFrame)

    number_of_samples = size(r10, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r10[:, :])
    describe(or10[:, :])
    describe(ar10[:, :])
    describe(tr10[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r10[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r10)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r10)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or10, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or10, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or10, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or10)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or10, tr10, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or10, tr10, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot10 = plot_acc_rates_against_λ_and_m_inf_scaled(or10, tr10, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot10, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r10[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r10[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r10[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r10[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r10[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r10[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_10_[Inf][mSE][arSE]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r10[:, 1:5]))
    println(describe(or10))
    println(describe(tr10))



###############################################################################
#### 11. Infection parameters - Strong Priors - mEI arEI Data Augmentation ####
###############################################################################


r11 = CSV.read("Inference/Results_11_[Inf][mEI][arEI]/res_11_[Inf][mEI][arEI].csv", DataFrame)
    or11 = CSV.read("Inference/Results_11_[Inf][mEI][arEI]/other_res_11_[Inf][mEI][arEI].csv", DataFrame)
    ar11 = CSV.read("Inference/Results_11_[Inf][mEI][arEI]/aug_res_11_[Inf][mEI][arEI].csv", DataFrame)
    tr11 = CSV.read("Inference/Results_11_[Inf][mEI][arEI]/tuning_res_11_[Inf][mEI][arEI].csv", DataFrame)
    ut11 = CSV.read("Inference/Results_11_[Inf][mEI][arEI]/update_tracker_11_[Inf][mEI][arEI].csv", DataFrame)

    number_of_samples = size(r11, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r11[:, :])
    describe(or11[:, :])
    describe(ar11[:, :])
    describe(tr11[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r11[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r11)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r11)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or11, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or11, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or11, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or11)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or11, tr11, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or11, tr11, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot11 = plot_acc_rates_against_λ_and_m_inf_scaled(or11, tr11, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot11, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r11[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r11[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r11[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r11[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r11[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r11[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_11_[Inf][mEI][arEI]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r11[:, 1:5]))
    println(describe(or11))
    println(describe(tr11))




############################################################################
#### 12. Detection parameters - Strong Priors - arDet Data Augmentation ####
############################################################################


r12 = CSV.read("Inference/Results_12_[Det][arDet]/res_12_[Det][arDet].csv", DataFrame)
    or12 = CSV.read("Inference/Results_12_[Det][arDet]/other_res_12_[Det][arDet].csv", DataFrame)
    ar12 = CSV.read("Inference/Results_12_[Det][arDet]/aug_res_12_[Det][arDet].csv", DataFrame)
    tr12 = CSV.read("Inference/Results_12_[Det][arDet]/tuning_res_12_[Det][arDet].csv", DataFrame)
    ut12 = CSV.read("Inference/Results_12_[Det][arDet]/update_tracker_12_[Det][arDet].csv", DataFrame)

    number_of_samples = size(r12, 1)

    ρ_tr = 0.75
    ρ_E_tr = 0.2

    param_names = ["ρ", "ρ_E"]

    describe(r12[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r12[:, 6:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r12)[1:5000, 6:7], param_names)
    chainplot12 = plot(chain_mix1)

    savefig(chainplot12, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r12)[5001:end, 6:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_det_posterior_compare(or12, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot12 = plot_det_posterior_compare(or12, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot12, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_det_posterior_compare(or12, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_det(or12)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_cumulative_acceptance_rate.pdf"))


    acc_plot12 = plot_acc_rates_against_λ_and_m_det(or12, tr12, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot12, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_det(or12, tr12, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_det(or12, tr12, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r12[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_Multi_Corr.pdf"))


    multi_kde_plot1 = @df r12[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_12_[Det][arDet]/", "detpars", "_Multi_KDE.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r12[:, 6:7]))
    println(describe(or12[:, 6:10]))
    println(describe(tr12))




##############################################################################
#### 13. Infection parameters - Strong Priors - arDeath Data Augmentation ####
##############################################################################


r13 = CSV.read("Inference/Results_13_[Inf][arDeath]/res_13_[Inf][arDeath].csv", DataFrame)
    or13 = CSV.read("Inference/Results_13_[Inf][arDeath]/other_res_13_[Inf][arDeath].csv", DataFrame)
    ar13 = CSV.read("Inference/Results_13_[Inf][arDeath]/aug_res_13_[Inf][arDeath].csv", DataFrame)
    tr13 = CSV.read("Inference/Results_13_[Inf][arDeath]/tuning_res_13_[Inf][arDeath].csv", DataFrame)
    ut13 = CSV.read("Inference/Results_13_[Inf][arDeath]/update_tracker_13_[Inf][arDeath].csv", DataFrame)

    number_of_samples = size(r13, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r13[:, :])
    describe(or13[:, :])
    describe(ar13[:, :])
    describe(tr13[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r13[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r13)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r13)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or13, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or13, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or13, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or13)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or13, tr13, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or13, tr13, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot13 = plot_acc_rates_against_λ_and_m_inf_scaled(or13, tr13, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot13, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r13[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r13[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r13[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r13[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r13[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r13[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_13_[Inf][arDeath]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r13[:, 1:5]))
    println(describe(or13))
    println(describe(tr13))



#############################################################################
#### 14. Infection parameters - Strong Priors - arPEnv Data Augmentation ####
#############################################################################


r14 = CSV.read("Inference/Results_14_[Inf][arPEnv]/res_14_[Inf][arPEnv].csv", DataFrame)
    or14 = CSV.read("Inference/Results_14_[Inf][arPEnv]/other_res_14_[Inf][arPEnv].csv", DataFrame)
    ar14 = CSV.read("Inference/Results_14_[Inf][arPEnv]/aug_res_14_[Inf][arPEnv].csv", DataFrame)
    tr14 = CSV.read("Inference/Results_14_[Inf][arPEnv]/tuning_res_14_[Inf][arPEnv].csv", DataFrame)
    ut14 = CSV.read("Inference/Results_14_[Inf][arPEnv]/update_tracker_14_[Inf][arPEnv].csv", DataFrame)

    number_of_samples = size(r14, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r14[:, :])
    describe(or14[:, :])
    describe(ar14[:, :])
    describe(tr14[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r14[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r14)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r14)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or14, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or14, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or14, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or14)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or14, tr14, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or14, tr14, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot14 = plot_acc_rates_against_λ_and_m_inf_scaled(or14, tr14, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot14, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r14[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r14[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r14[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r14[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r14[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r14[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_14_[Inf][arPEnv]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r14[:, 1:5]))
    println(describe(or14))
    println(describe(tr14))




##############################################################################
#### 15. Infection parameters - Strong Priors - arMoves Data Augmentation ####
##############################################################################


r15 = CSV.read("Inference/Results_15_[Inf][arMoves]/res_15_[Inf][arMoves].csv", DataFrame)
    or15 = CSV.read("Inference/Results_15_[Inf][arMoves]/other_res_15_[Inf][arMoves].csv", DataFrame)
    ar15 = CSV.read("Inference/Results_15_[Inf][arMoves]/aug_res_15_[Inf][arMoves].csv", DataFrame)
    tr15 = CSV.read("Inference/Results_15_[Inf][arMoves]/tuning_res_15_[Inf][arMoves].csv", DataFrame)
    ut15 = CSV.read("Inference/Results_15_[Inf][arMoves]/update_tracker_15_[Inf][arMoves].csv", DataFrame)

    number_of_samples = size(r15, 1)

    β_c_tr = 0.002
    β_b_tr = 0.004
    γ_tr = 0.015
    F_tr = 0.004
    ϵ_tr = 0.05

    param_names = ["β_c", "β_b", "γ", "F", "ϵ"]

    describe(r15[:, :])
    describe(or15[:, :])
    describe(ar15[:, :])
    describe(tr15[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r15[:, 1:5]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r15)[1:5000, 1:5], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r15)[5001:end, 1:5], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or15, 1, 5000, string("Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_diagnostic_likelihoods_tuning.pdf"))

    post_plot2 = plot_inf_posterior_compare(or15, 5001, number_of_samples, string("Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_diagnostic_likelihoods_post_tuning.pdf"))

    post_plot3 = plot_inf_posterior_compare(or15, 1, number_of_samples, string("Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_diagnostic_likelihoods_full.pdf"))


    acc_ = calc_acc_inf(or15)

    acc_plot = plot_acc(acc_, string("Cumulative Acceptance Rates"))

    savefig(acc_plot, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_cumulative_acceptance_rate.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or15, tr15, 1, 5000, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_batch_acceptance_rate_tuning.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or15, tr15, 5001, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_batch_acceptance_rate_post_tuning.pdf"))

    acc_plot15 = plot_acc_rates_against_λ_and_m_inf_scaled(or15, tr15, 1, number_of_samples, number_of_samples, string("Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot15, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_batch_acceptance_rate_full.pdf"))


    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r15[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r15[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r15[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_Corr_all.pdf"))


    multi_kde_plot1 = @df r15[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r15[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r15[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_15_[Inf][arMoves]/", "InfPars", "_Multi_KDE_all.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r15[:, 1:5]))
    println(describe(or15))
    println(describe(tr15))


####################################################################
#### 16. All parameters - Strong Priors - All Data Augmentation ####
####################################################################

r16 = CSV.read("Inference/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/res_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    or16 = CSV.read("Inference/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/other_res_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ar16 = CSV.read("Inference/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/aug_res_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    tr16 = CSV.read("Inference/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/tuning_res_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ut16 = CSV.read("Inference/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/update_tracker_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)

    number_of_samples = size(r16, 1)

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

    describe(r16[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r16[:, 1:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r16)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r16)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or16, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or16, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or16, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or16, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or16, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or16, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or16)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or16, tr16, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or16, tr16, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or16, tr16, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or16)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or16, tr16, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or16, tr16, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or16, tr16, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r16[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r16[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r16[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r16[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r16[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r16[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r16[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r16[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Results_16_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r16[:, :]))
    println(describe(or16))
    println(describe(tr16))


##################################################################
#### 17. All parameters - Weak Priors - All Data Augmentation ####
##################################################################

r17 = CSV.read("Inference/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/res_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    or17 = CSV.read("Inference/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/other_res_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ar17 = CSV.read("Inference/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/aug_res_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    tr17 = CSV.read("Inference/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/tuning_res_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ut17 = CSV.read("Inference/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/update_tracker_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)

    number_of_samples = size(r17, 1)

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

    describe(r17[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r17[:, 1:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r17)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r17)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or17, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or17, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or17, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or17, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or17, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or17, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or17)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or17, tr17, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or17, tr17, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or17, tr17, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or17)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or17, tr17, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or17, tr17, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or17, tr17, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r17[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r17[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r17[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r17[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r17[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r17[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r17[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r17[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Results_17_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r17[:, :]))
    println(describe(or17))
    println(describe(tr17))

##################################################################
#### 18. All parameters - Mid Priors - All Data Augmentation ####
##################################################################

r18 = CSV.read("Inference/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/res_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    or18 = CSV.read("Inference/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/other_res_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ar18 = CSV.read("Inference/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/aug_res_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    tr18 = CSV.read("Inference/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/tuning_res_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ut18 = CSV.read("Inference/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/update_tracker_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)

    number_of_samples = size(r18, 1)

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

    describe(r18[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r18[:, 1:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r18)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r18)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or18, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or18, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or18, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or18, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or18, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or18, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or18)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or18, tr18, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or18, tr18, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or18, tr18, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or18)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or18, tr18, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or18, tr18, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or18, tr18, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r18[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r18[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r18[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r18[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r18[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r18[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r18[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r18[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Results_18_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r18[:, :]))
    println(describe(or18))
    println(describe(tr18))

###################################################################
#### 19. All parameters - Wrong Priors - All Data Augmentation ####
###################################################################

r19 = CSV.read("Inference/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/res_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    or19 = CSV.read("Inference/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/other_res_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ar19 = CSV.read("Inference/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/aug_res_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    tr19 = CSV.read("Inference/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/tuning_res_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)
    ut19 = CSV.read("Inference/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/update_tracker_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves].csv", DataFrame)

    number_of_samples = size(r19, 1)

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

    describe(r19[:, :])

    ### SAMPLES ###

    chain_full = Chains(Array(r19[:, 1:7]), param_names)
    chainplot1 = plot(chain_full)

    savefig(chainplot1, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_full.pdf"))

    chain_mix1 = Chains(Array(r19)[1:5000, 1:7], param_names)
    chainplot2 = plot(chain_mix1)

    savefig(chainplot2, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_tuning.pdf"))

    chain_mix2 = Chains(Array(r19)[5001:end, 1:7], param_names)
    chainplot3 = plot(chain_mix2)

    savefig(chainplot3, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_chain_post_tuning.pdf"))

    ### PLOT SETTINGS ###

    Plots.scalefontsizes(3)

    ### AUXILARY ###

    post_plot1 = plot_inf_posterior_compare(or19, 1, 5000, string("Infection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot1, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_inf.pdf"))

    post_plot2 = plot_inf_posterior_compare(or19, 5001, number_of_samples, string("Infection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot2, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_inf.pdf"))

    post_plot3 = plot_inf_posterior_compare(or19, 1, number_of_samples, string("Infection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot3, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_inf.pdf"))


    post_plot4 = plot_det_posterior_compare(or19, 1, 5000, string("Detection Parameters Posterior Comparison Mixture 1"), :bottomright, 0.5)

    savefig(post_plot4, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_tuning_det.pdf"))

    post_plot5 = plot_det_posterior_compare(or19, 5001, number_of_samples, string("Detection Parameters Posterior Comparison Mixture 2"), :bottomright, 0.5)

    savefig(post_plot5, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_post_tuning_det.pdf"))

    post_plot6 = plot_det_posterior_compare(or19, 1, number_of_samples, string("Detection Parameters Posterior Comparison Full"), :bottomright, 0.5)

    savefig(post_plot6, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_diagnostic_likelihoods_full_det.pdf"))



    acc_inf = calc_acc_inf(or19)

    acc_plot1 = plot_acc(acc_inf, string("Infection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot1, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_inf.pdf"))


    acc_plot2 = plot_acc_rates_against_λ_and_m_inf_scaled(or19, tr19, 1, 5000, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot2, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_inf.pdf"))

    acc_plot3 = plot_acc_rates_against_λ_and_m_inf_scaled(or19, tr19, 5001, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot3, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_inf.pdf"))

    acc_plot4 = plot_acc_rates_against_λ_and_m_inf_scaled(or19, tr19, 1, number_of_samples, number_of_samples, string("Infection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot4, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_inf.pdf"))




    acc_det = calc_acc_det(or19)

    acc_plot5 = plot_acc(acc_det, string("Detection Parameters Cumulative Acceptance Rates"))

    savefig(acc_plot5, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_cumulative_acceptance_rate_det.pdf"))


    acc_plot6 = plot_acc_rates_against_λ_and_m_det_scaled(or19, tr19, 1, 5000, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 1"), :topright, false)

    savefig(acc_plot6, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_tuning_det.pdf"))

    acc_plot7 = plot_acc_rates_against_λ_and_m_det_scaled(or19, tr19, 5001, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Mixture 2"), :topright, false)

    savefig(acc_plot7, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_post_tuning_det.pdf"))

    acc_plot8 = plot_acc_rates_against_λ_and_m_det_scaled(or19, tr19, 1, number_of_samples, number_of_samples, string("Detection Parameters Acceptance Rates and Tuning Parameters Full"), :topright, true)

    savefig(acc_plot8, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_batch_acceptance_rate_full_det.pdf"))



    ### CORR AND MUTLI KDE PLOTS ###

    corr_plot1 = corrplot(Array(r19[5001:end, 1:3]), label = ["β_c", "β_b", "γ"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot1, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbg.pdf"))

    corr_plot2 = corrplot(Array(r19[5001:end, [1,2,4]]), label = ["β_c", "β_b", "F"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot2, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_bbf.pdf"))

    corr_plot3 = corrplot(Array(r19[5001:end, [1,2,3,4,5]]), label = ["β_c", "β_b", "γ", "F", "ϵ"], fc = :thermal, ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot3, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_inf.pdf"))


    corr_plot4 = corrplot(Array(r19[5001:end, 6:7]), label = ["ρ", "ρ_E"], fc = :thermal, ticks = true, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(corr_plot4, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_Corr_det.pdf"))



    multi_kde_plot1 = @df r19[5001:end, 1:3] cornerplot([:β_c, :β_b, :γ], label = ["β_c", "β_b", "γ"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot1, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbg.pdf"))

    multi_kde_plot2 = @df r19[5001:end, [1,2,4]] cornerplot([:β_c, :β_b, :F], label = ["β_c", "β_b", "F"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot2, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_bbf.pdf"))

    multi_kde_plot3 = @df r19[5001:end, [1,2,3,4,5]] cornerplot([:β_c, :β_b, :γ, :F, :ϵ], label = ["β_c", "β_b", "γ", "F", "ϵ"], ticks = false, size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot3, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_inf.pdf"))


    multi_kde_plot4 = @df r19[5001:end, 6:7] cornerplot([:ρ, :ρ_E], label = ["ρ", "ρ_E"], size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)

    savefig(multi_kde_plot4, string("Visualisation/Results_19_[Inf][Det][mSE][mEI][arSE][arEI][arDet][arDeath][arPEnv][arMoves]/", "allpars", "_Multi_KDE_det.pdf"))


    ### PLOT SETTINGS ###

    Plots.resetfontsizes()


    println(describe(r19[:, :]))
    println(describe(or19))
    println(describe(tr19))
