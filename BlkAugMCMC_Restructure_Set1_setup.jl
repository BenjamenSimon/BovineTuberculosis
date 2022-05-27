using NamedArrays
using JLD2

#################
### LOAD DATA ###
#################

combi_array = load("Restructure/Data/Set 1/combi_array.jld2")["array"]

f_to_p_dict = load("Restructure/Data/Set 1/f_to_p_dict.jld2")["dict"]


###########
### RUN ###
###########

β_c_tr = 0.002
β_b_tr = 0.004
γ_tr = 0.015

F_tr = 0.004
ϵ_tr = 0.05

ρ_tr = 0.75
ρ_E_tr = 0.2

θ_bb_tr = 0.25/52
θ_bd_tr = 0.25/52

epi_params_true = [β_c_tr, β_b_tr, γ_tr, F_tr, ϵ_tr, ρ_tr, ρ_E_tr, θ_bb_tr, θ_bd_tr]

d_β_c = Uniform(0, 0.04)
d_β_b = Uniform(0, 0.08)
d_γ = Uniform(0, 0.3)
d_F = Uniform(0, 0.08)
d_ϵ = Uniform(0, 0.1)
d_ρ = Uniform(0, 1)
d_ρ_E = Uniform(0, 1)

epi_params_dists = [d_β_c, d_β_b, d_γ, d_F, d_ϵ, d_ρ, d_ρ_E]
