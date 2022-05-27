using Distributions
using InvertedIndices
using Random
using DataFrames
using DataFramesMeta
using LinearAlgebra
using CSV
using RData
using NamedArrays
using JLD2

#################
### LOAD DATA ###
#################

results_df = CSV.read("Data/XtraM Sim/Set 1/results_df.csv", DataFrame)

track_df = CSV.read("Data/XtraM Sim/Set 1/track_df.csv", DataFrame)

pers_df = CSV.read("Data/XtraM Sim/Set 1/pers_df.csv", DataFrame)

record_of_movements = CSV.read("Data/XtraM Sim/Set 1/record_of_movements.csv", DataFrame)

cph_to_formatted = load("Data/cph_to_formatted_chesh.rds")



########################
### Create 3D Arrays ###
########################

#### RESULTS ####

init_res_3D = zeros(2172, 360, 34)

res_3D = NamedArray(init_res_3D)

setdimnames!(res_3D, ("Farm", "Time", "Data"))

setnames!(res_3D, [
  "t",
  "farm_uid",
  "pos",
  "cS_init",
  "cE_init",
  "cI_init",
  "cS_Moves",
  "cE_Moves",
  "cI_Moves",
  "cS_postM",
  "cE_postM",
  "cI_postM",
  "cS_postEI",
  "cE_postEI",
  "cI_postEI",
  "cS_postDet",
  "cE_postDet",
  "cI_postDet",
  "cS_final",
  "cE_final",
  "cI_final",
  "bS_init",
  "bE_init",
  "bI_init",
  "bS_postEI",
  "bE_postEI",
  "bI_postEI",
  "bS_final",
  "bE_final",
  "bI_final",
  "move_res_init",
  "move_res_final",
  # "pcS_init",
  # "pcE_init",
  # "pcI_init",
  # "pcS_final",
  # "pcE_final",
  # "pcI_final",
  # "pbS_init",
  # "pbE_init",
  # "pbI_init",
  # "pbS_final",
  # "pbE_final",
  # "pbI_final",
  "county",
  "parish",
  # "holding",
  # "county_idx",
  # "parish_idx",
  # "holding_idx",
], 3)

for farm in 1:2172
    res_3D[farm, :, [1:2 ; 4:34]] = Array(@subset(results_df, :row_id .== farm))[:, [1:31 ; 44:45]]
end


#### TRACK ####

init_track_3D = zeros(2172, 360, 28)

track_3D = NamedArray(init_track_3D)

setdimnames!(track_3D, ("Farm", "Time", "Data"))

setnames!(track_3D,
  [
    "t",
    "row_id",
    "pos",
    "sus_on",
    "exp_on",
    "inf_on",
    "sus_off",
    "exp_off",
    "inf_off",
    "S_on_out",
    "E_on_out",
    "I_on_out",
    "c_new_exp",
    "c_new_inf",
    "b_new_exp",
    "b_new_inf",
    "test_occur",
    "next_test_week",
    "E_detected",
    "I_detected",
    "c_birth",
    "c_S_death",
    "c_E_death",
    "c_I_death",
    "b_birth",
    "b_S_death",
    "b_E_death",
    "b_I_death",
  ],
  3,
)

for farm in 1:2172
    track_3D[farm, :, [1:2 ; 4:28]] = Array(@subset(track_df, :row_id .== farm))[:, 1:27]
end

track_3D[1,100,:]


#### PERSISTENTS ####

init_pers_3D = zeros(2172, 360, 9)

pers_3D = NamedArray(init_pers_3D)

setdimnames!(pers_3D, ("Farm", "Time", "Data"))

setnames!(pers_3D,
  [
    "t",
    "row_id",
    "pos",
    "p_env_prev",
    "p_env_cur",
    "c_exp_prob",
    "b_exp_prob",
    "c_inf_prob",
    "b_inf_prob",
  ],
  3,
)

for farm in 1:2172
    pers_3D[farm, :, [1:2 ; 4:9]] = Array(@subset(pers_df, :row_id .== farm))[:, 1:8]
end

pers_3D[1,100,:]


# #### PARISH ####
#
# first_h_of_p = combine(first, groupby(cph_to_formatted, [:county_idx, :parish_idx]))[:, :row_id]
#
# function calc_area_of_parish()
#
#   num_of_parishes = size(first_h_of_p, 1)
#
#   area_of_parish = fill(-1, num_of_parishes)
#
#   for i in 1:(num_of_parishes-1)
#     area_of_parish[i] = (first_h_of_p[(i+1)] - first_h_of_p[(i)]) * 1000
#   end
#
#   area_of_parish[num_of_parishes] = (size(cph_to_formatted, 1) - first_h_of_p[num_of_parishes] + 1) * 1000
#
#   return(area_of_parish)
# end
#
# area_of_parish = calc_area_of_parish()
#
#
# init_parish_3D = zeros(311, 360, 20)
#
# parish_3D = NamedArray(init_parish_3D)
#
# setdimnames!(parish_3D, ("Parish", "Time", "Data"))
#
# setnames!(parish_3D, [
#   "t",
#   "parish_uid",
#   "parish_pos",
#   "pcS_init",
#   "pcE_init",
#   "pcI_init",
#   "pcS_final",
#   "pcE_final",
#   "pcI_final",
#   "pbS_init",
#   "pbE_init",
#   "pbI_init",
#   "pbS_final",
#   "pbE_final",
#   "pbI_final",
#   "remaining_pressure",
#   "new_pressure",
#   "scaling",
#   "county",
#   "parish",
# ], 3)
#
# for parish in 1:311
#     parish_3D[parish, :, 2] .= parish
#     parish_3D[parish, :, [1 ; 4:15 ; 19:20]] = Array(@subset(results_df, :row_id .== first_h_of_p[parish]))[:, [1 ; 32:45]]
#     parish_3D[parish, :, 16:17] = Array(@subset(track_df, :row_id .== first_h_of_p[parish]))[:, 28:29]
#     parish_3D[parish, :, 18] .= area_of_parish[parish]
# end
#
# parish_3D[1,100,:]
#
# #### Parish Dictionary ####
#
# f_to_p_dict = Dict{Integer, Array}()
#
# for farm_uid in 1:size(cph_to_formatted, 1)
#
#   c_p = Array(@subset(results_df, :row_id .== farm_uid, :t .== 1)[:, 44:45])
#
#   parish_uid = parish_3D[(parish_3D[:, 1, 19] .== c_p[1]) .& (parish_3D[:, 1, 20] .== c_p[2]), 1, "parish_uid"][1]
#
#   parish_members_uid = Array(@subset(results_df, :t .== 1, :county .== c_p[1], :parish .== c_p[2])[:, :row_id])
#
#   parish_pos = 0
#   parish_members_pos = 0
#
#   f_to_p_dict[farm_uid] = [parish_uid, parish_pos, parish_members_uid, parish_members_pos]
# end

#### Movements ####

# movements_grouped = Array{NamedArray{Float64, 2}}(undef, 2173, 360)
#
# for time in 1:36
#   for farm in 1:100
#
#     move_subset_df = @where(record_of_movements, :off_row_id .== farm, :t .== time)
#     move_subset = NamedArray(zeros(size(move_subset_df, 1), 12))
#     println(move_subset_df)
#     move_subset[:, [1:3 ; 6:12]] .= move_subset_df
#     setnames!(move_subset,
#       [
#         "t",
#         "off_farm_uid",
#         "on_farm_uid",
#         "off_farm_pos",
#         "on_farm_pos",
#         "cS_gen",
#         "cE_gen",
#         "cI_gen",
#         "cS_StateMoved",
#         "cE_StateMoved",
#         "cI_StateMoved",
#         "total_moves_in_transaction",
#       ],
#       2,
#     )
#
#     movements_grouped[farm, time] = move_subset
#   end
#
#   move_subset_df = @where(record_of_movements, :off_row_id .== -1, :t .== time)
#   move_subset = NamedArray(zeros(size(move_subset_df, 1), 12))
#   move_subset[:, [1:3 ; 6:12]] .= move_subset_df
#   setnames!(move_subset,
#     [
#       "t",
#       "off_farm_uid",
#       "on_farm_uid",
#       "off_farm_pos",
#       "on_farm_pos",
#       "cS_gen",
#       "cE_gen",
#       "cI_gen",
#       "cS_StateMoved",
#       "cE_StateMoved",
#       "cI_StateMoved",
#       "total_moves_in_transaction",
#     ],
#     2,
#   )
#
#   movements_grouped[1273, time] = move_subset
# end
#
# NOTE: Doesn't work because there are many times when a farm doesn't move
# NOTE: Just keeping it as it is for now until I think of a better solution
# NOTE: but I think any solution may just add unnecessary complexity




################################
### Subset and Add Positions ###
################################

sum_of_infecteds = sum(res_3D[:, :, "cI_init"], dims = 2)[:, 1]

ids_oi = sortperm(Array{Int64}(sum_of_infecteds), rev = true)[1:100]

# sum_of_infecteds[ids_oi]

res_3D[ids_oi, :, "pos"] .= 1:length(ids_oi)
track_3D[ids_oi, :, "pos"] .= 1:length(ids_oi)
pers_3D[ids_oi, :, "pos"] .= 1:length(ids_oi)

res_3D_oi = res_3D[ids_oi, :, :]
track_3D_oi = track_3D[ids_oi, :, :]
pers_3D_oi = pers_3D[ids_oi, :, :]


#### PARISH ####

results_df_oi = @subset(results_df, in.(:row_id, [ids_oi]), :t .== 1)

first_h_of_p = combine(first, groupby(results_df_oi, [:county_idx, :parish_idx]))[:, :row_id]

function calc_area_of_parish()

  num_of_parishes = size(first_h_of_p, 1)

  area_of_parish = fill(-1, num_of_parishes)

  for i in 1:num_of_parishes

    county_idx = @subset(results_df_oi, :row_id .== first_h_of_p[i])[:, 47]
    parish_idx = @subset(results_df_oi, :row_id .== first_h_of_p[i])[:, 48]

    num_of_farms = size(@subset(results_df, :county_idx .== county_idx, :parish_idx .== parish_idx, :t .== 1), 1)

    area_of_parish[i] = num_of_farms * 1000
  end

  return(area_of_parish)
end

area_of_parish = calc_area_of_parish()

init_parish_3D = zeros(size(area_of_parish, 1), 360, 20)

parish_3D = NamedArray(init_parish_3D)

setdimnames!(parish_3D, ("Parish", "Time", "Data"))

setnames!(parish_3D, [
  "t",
  "parish_uid",
  "parish_pos",
  "pcS_init",
  "pcE_init",
  "pcI_init",
  "pcS_final",
  "pcE_final",
  "pcI_final",
  "pbS_init",
  "pbE_init",
  "pbI_init",
  "pbS_final",
  "pbE_final",
  "pbI_final",
  "remaining_pressure",
  "new_pressure",
  "scaling",
  "county",
  "parish",
], 3)

for parish in 1:size(area_of_parish, 1)
    parish_3D[parish, :, 2] = Array(@subset(results_df, :row_id .== first_h_of_p[parish]))[:, 45]
    parish_3D[parish, :, 3] .= parish
    parish_3D[parish, :, [1 ; 4:15 ; 19:20]] = Array(@subset(results_df, :row_id .== first_h_of_p[parish]))[:, [1 ; 32:45]]
    parish_3D[parish, :, 16:17] = Array(@subset(track_df, :row_id .== first_h_of_p[parish]))[:, 28:29]
    parish_3D[parish, :, 18] .= area_of_parish[parish]
end

parish_3D[1,100,:]

#### Parish Dictionary ####

f_to_p_dict = Dict{Integer, Array}()

for farm_uid in 1:size(ids_oi, 1)

  c_p = Array(@subset(results_df, :row_id .== ids_oi[farm_uid], :t .== 1)[:, 44:45])

  parish_uid, parish_pos = parish_3D[(parish_3D[:, 1, 19] .== c_p[1]) .& (parish_3D[:, 1, 20] .== c_p[2]), 1, 2:3]

  parish_members_uid = Array(res_3D_oi[(res_3D_oi[:, 1, 33] .== c_p[1]) .& (res_3D_oi[:, 1, 34] .== c_p[2]), 1, :][:, 2])

  parish_members_pos = Array(res_3D_oi[(res_3D_oi[:, 1, 33] .== c_p[1]) .& (res_3D_oi[:, 1, 34] .== c_p[2]), 1, :][:, 3])

  f_to_p_dict[ids_oi[farm_uid]] = [parish_uid, parish_pos, parish_members_uid, parish_members_pos]
end





##########################
### Combine into Array ###
##########################

combi_array = [res_3D_oi, track_3D_oi, pers_3D_oi, parish_3D]

f_to_p_dict


####################
### Save Results ###
####################

save("Restructure/Data/Set 1/res_3D_oi.jld2", "array", res_3D_oi)

save("Restructure/Data/Set 1/combi_array.jld2", "array", combi_array)
save("Restructure/Data/Set 1/f_to_p_dict.jld2", "dict", f_to_p_dict)
