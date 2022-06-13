
include("LoadPackages.jl")

#####################
### Plot Settings ###
#####################

gr(size = (3000, 1800), bottom_margin = 20mm, left_margin = 25mm, right_margin = 55mm, dpi = 300)
Plots.scalefontsizes(3)

#################
### LOAD DATA ###
#################

results_df = CSV.read("Data/Set 1/results_df.csv", DataFrame)
track_df = CSV.read("Data/Set 1/track_df.csv", DataFrame)
pers_df = CSV.read("Data/Set 1/pers_df.csv", DataFrame)

##################
### Parameters ###
##################

β_c_tr = 0.0002
β_b_tr = 0.0004
γ_tr = 0.015

F_tr = 0.0004
ϵ_tr = 0.05

ρ_tr = 0.75
ρ_E_tr = 0.2

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52


#######################
### Epidemic Curves ###
#######################

### All Farms ###

temp_df1 = groupby(results_df[:, [:t, :cS_final, :cE_final, :cI_final]], [:t])
temp_df2 = combine(temp_df1, :cS_final => sum, :cE_final => sum, :cI_final => sum)
ys_3 = [temp_df2[:, :cS_final_sum], temp_df2[:, :cE_final_sum], temp_df2[:, :cI_final_sum]]

labels = ["Sus. Cattle" "Exp. Cattle" "Inf. Cattle"]

StatsPlots.plot(ys_3, title = "Epidemic Curves", color = [:green :blue :red],
                legend = :right, lw = 3, label = labels,
                xlabel = "Time", ylabel = "Number of Cattle")
StatsPlots.plot(ys_3[2:3], title = "Epidemic Curves", color = [:blue :red],
                legend = :topleft, lw = 3, label = labels[:, 2:3],
                xlabel = "Time", ylabel = "Number of Cattle")

#############################
### Infection Probability ###
#############################

### All Farms ###

temp_df3 = leftjoin(results_df[:, [:t, :row_id, :cS_final, :cE_final, :cI_final]],
                    pers_df[:, [:t, :row_id, :p_env_prev]], on = [:t, :row_id])
temp_df4 = leftjoin(temp_df3, track_df[:, [:t, :row_id, :c_new_exp, :c_new_inf]],
                    on = [:t, :row_id])
temp_df5 = @transform(temp_df4, cN_final = :cS_final .+ :cE_final .+ :cI_final)
temp_df6 = @transform(temp_df5, cProp_final = :cI_final ./ :cN_final)
replace!(temp_df6.cProp_final, NaN => 0)
temp_df7 = @transform(temp_df6, c_beta_cont = :cProp_final .* β_c_tr,
                                c_parishenv_cont = :p_env_prev .* F_tr)
temp_df8 = groupby(temp_df7, [:t])
temp_df9 = @combine(temp_df8, beta_comp_sum = sum(:c_beta_cont),
                              F_comp_sum = sum(:c_parishenv_cont),
                              new_E_sum = sum(:c_new_exp),
                              new_I_sum = sum(:c_new_inf))

ys_4 = [temp_df9[:, :beta_comp_sum] temp_df9[:, :F_comp_sum]]
ctg = repeat(["β_c", "F"], inner = 360)

groupedbar(ys_4, bar_position = :stack, bar_width=0.5,lw = 0, alpha = 0.5,
            group = ctg, xlabel = "Time Step", ylabel = "Contribution to Infectious Pressure")
    plot!(twinx(), [temp_df9[:, :new_E_sum]],
        linecolor = [:red], lw = 4, linealpha = 1,
        legend = :topleft, ylabel = "Number of Cattle",
        label = "New Exposed")

############################
### Number of detections ###
############################

### All Farms ###

temp_df10 = groupby(track_df[:, [:t, :E_detected, :I_detected]], [:t])
temp_df11 = combine(temp_df10, :E_detected => sum, :I_detected => sum)
ys_5 = [temp_df11[:, :E_detected_sum] temp_df11[:, :I_detected_sum]]

ctg = repeat(["E", "I"], inner = 360)
groupedbar(ys_5, bar_position = :dodge, bar_width=0.5, lw = 0, alpha = 0.8,
            group = ctg, xlabel = "Time Step", ylabel = "Number of Cattle")
