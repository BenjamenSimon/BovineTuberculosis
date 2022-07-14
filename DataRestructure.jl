
include("LoadPackages.jl")

#################
### LOAD DATA ###
#################

results_df = CSV.read("Data/Set 1/results_df.csv", DataFrame)

track_df = CSV.read("Data/Set 1/track_df.csv", DataFrame)

pers_df = CSV.read("Data/Set 1/pers_df.csv", DataFrame)

record_of_movements = CSV.read("Data/Set 1/record_of_movements.csv", DataFrame)

cph_to_formatted = load("Data/cph_to_formatted_chesh.rds")






########################
### Create 3D Arrays ###
########################

######
### RESULTS
######

init_res_3D = Int64.(zeros(2172, 360, 34))

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


######
### TRACK
######

init_track_3D = Int64.(zeros(2172, 360, 28))

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


#########
### PERSISTENTS
#########

init_pers_3D = zeros(2172, 360, 7)

pers_3D = NamedArray(init_pers_3D)

setdimnames!(pers_3D, ("Farm", "Time", "Data"))

setnames!(pers_3D,
  [
    "t",
    "row_id",
    "pos",
    "c_exp_prob",
    "b_exp_prob",
    "c_inf_prob",
    "b_inf_prob",
  ],
  3,
)

for farm in 1:2172
    pers_3D[farm, :, [1:2 ; 4:7]] = Array(@subset(pers_df, :row_id .== farm))[:, [1:2 ; 5:8]]
end


#######
### PARISH
#######

first_h_of_p_ALL = Int64.(combine(first, groupby(results_df, [:county_idx, :parish_idx]))[:, :row_id])


function calc_area_of_parish_all()

  num_of_parishes = size(first_h_of_p_ALL, 1)

  area_of_parish = fill(-1, num_of_parishes)

  for i in 1:num_of_parishes

    county_idx = @subset(results_df, :row_id .== first_h_of_p_ALL[i])[1, 47]
    parish_idx = @subset(results_df, :row_id .== first_h_of_p_ALL[i])[1, 48]

    num_of_farms = size(@subset(results_df, :county_idx .== county_idx, :parish_idx .== parish_idx, :t .== 1), 1)

    area_of_parish[i] = num_of_farms * 1000
  end

  return(area_of_parish)
end

area_of_parish_all = calc_area_of_parish_all()


init_parish_3D = zeros(size(area_of_parish_all, 1), 360, 23)

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
  "p_env_prev",
  "p_env_cur",
  "county",
  "parish",
  "first_h_of_p_uid"
], 3)


for parish in 1:size(area_of_parish_all, 1)
    parish_3D[parish, :, 2] .= parish
    parish_3D[parish, :, [1 ; 4:15 ; 21:22]] = Array(@subset(results_df, :row_id .== first_h_of_p_ALL[parish]))[:, [1 ; 32:45]]
    parish_3D[parish, :, 16:17] = Array(@subset(track_df, :row_id .== first_h_of_p_ALL[parish]))[:, 28:29]
    parish_3D[parish, :, 18] .= area_of_parish_all[parish]
    parish_3D[parish, :, 19:20] = Array(@subset(pers_df, :row_id .== first_h_of_p_ALL[parish]))[:, [4,5]]
    parish_3D[parish, :, 23] .= first_h_of_p_ALL[parish]
end








################################
### Subset and Add Positions ###
################################

#############
### Identify top X parishes
#############

results_df_first_h_of_p = @subset(results_df, in.(:row_id, [first_h_of_p_ALL]))[:, [1:2 ; 32:37]]

sum_of_infecteds_parish = combine(groupby(results_df_first_h_of_p, [:row_id]), :pcI_init => sum)

uids_first_h_of_p_top = Int64.(sum_of_infecteds_parish[sortperm(Array{Int64}(sum_of_infecteds_parish[:, 2]), rev = true)[1:7], :][:, 1])


##################
#### Full Farm to Parish Dictionary
##################

f_to_p_dict_ALL = Dict{Integer, Array}()

for farm_uid in uids_first_h_of_p_top

  c_p = Array(@subset(results_df, :row_id .== farm_uid, :t .== 1)[:, 44:45])

  parish_uid = parish_3D[(parish_3D[:, 1, 21] .== c_p[1]) .& (parish_3D[:, 1, 22] .== c_p[2]), 1, 2][1]

  parish_members_uid = Array(res_3D[(res_3D[:, 1, 33] .== c_p[1]) .& (res_3D[:, 1, 34] .== c_p[2]), 1, :][:, 2])

  f_to_p_dict_ALL[farm_uid] = [Int64.(parish_uid), Int64.(parish_members_uid)]
end

f_to_p_dict_ALL


##############
#### Extract Farms of Interest
##############

uids_oi = Array[]
parish_uids_oi = Array[]

for i in 1:size(uids_first_h_of_p_top, 1)
  uids_oi = [uids_oi ; f_to_p_dict_ALL[uids_first_h_of_p_top[i]][2]]
  parish_uids_oi = [parish_uids_oi ; f_to_p_dict_ALL[uids_first_h_of_p_top[i]][1]]
end

uids_oi
parish_uids_oi


##########
### Subset arrays
##########

res_3D[uids_oi, :, "pos"] .= 1:length(uids_oi)
track_3D[uids_oi, :, "pos"] .= 1:length(uids_oi)
pers_3D[uids_oi, :, "pos"] .= 1:length(uids_oi)
parish_3D[parish_uids_oi, :, "parish_pos"] .= 1:length(parish_uids_oi)

res_3D_oi = res_3D[uids_oi, :, :]
track_3D_oi = track_3D[uids_oi, :, :]
pers_3D_oi = pers_3D[uids_oi, :, :]
parish_3D_oi = parish_3D[parish_uids_oi, :, :]


#################
### Update the Farm to Parish Dictionary
#################

f_to_p_dict = Dict{Integer, Array}()

for pos in 1:size(uids_oi, 1)

  c_p = Array(@subset(results_df, :row_id .== uids_oi[pos], :t .== 1)[:, 44:45])

  parish_uid, parish_pos = parish_3D[(parish_3D[:, 1, 21] .== c_p[1]) .& (parish_3D[:, 1, 22] .== c_p[2]), 1, 2:3]

  parish_members_uid = Array(res_3D_oi[(res_3D_oi[:, 1, 33] .== c_p[1]) .& (res_3D_oi[:, 1, 34] .== c_p[2]), 1, :][:, 2])

  parish_members_pos = Array(res_3D_oi[(res_3D_oi[:, 1, 33] .== c_p[1]) .& (res_3D_oi[:, 1, 34] .== c_p[2]), 1, :][:, 3])

  f_to_p_dict[pos] = [Int64.(parish_uid), Int64.(parish_pos), Int64.(parish_members_uid), Int64.(parish_members_pos), Int64.(uids_oi[pos])]
            #farm_pos                                                                                                         # farm_uid
end

f_to_p_dict


####################
### Create the farm_uid to position Dictionary
####################

ids_to_pos_dict = Dict{Int64, Int64}()

for pos in 1:size(uids_oi, 1)
    ids_to_pos_dict[res_3D_oi[pos, 1, 2]] = res_3D_oi[pos, 1, 3]
end

ids_to_pos_dict






#################
### MOVEMENTS ###
#################

#############
### Subset the Array of Movements
#############

record_of_movements

record_of_movements_temp = deepcopy(Int64.(record_of_movements))

for i in 1:size(record_of_movements_temp, 1)
  if ((record_of_movements_temp[i,2] in uids_oi) == false)
     record_of_movements_temp[i,2] = -1
     record_of_movements_temp[i,[4,5,6]] = [10800, 600, 600]
  end
  if ((record_of_movements_temp[i,3] in uids_oi) == false)
     record_of_movements_temp[i,3] = -1
  end
end

for t in 1:360
  indicies = findall((record_of_movements_temp[:,1] .== t) .& (record_of_movements_temp[:,2] .== -1) .& (record_of_movements_temp[:,3] .== -1))

  record_of_movements_temp[indicies,2] .= -99
  record_of_movements_temp[indicies,3] .= -99
end

record_of_movements_temp

record_of_movements_oi = record_of_movements_temp[(record_of_movements_temp[:,2] .> -10), :]


###########
### Movement dictionary
###########

dict_of_movements = Dict{Tuple{Int64, Int64}, Array}()

for pos in 1:size(uids_oi, 1)
  for t in 1:360
    dict_of_movements[(pos, t)] = findall((record_of_movements_oi[:,1] .== t) .& (record_of_movements_oi[:,2] .== uids_oi[pos]))
  end
end

dict_of_movements


dict_of_movements_out = Dict{Tuple{Int64, Int64, Int64}, Array}()

for t in 1:360
  for pos in 1:size(uids_oi, 1)
    dict_of_movements_out[(-1, t, pos)] = findall((record_of_movements_oi[:,1] .== t) .& (record_of_movements_oi[:,2] .== -1) .& (record_of_movements_oi[:,3] .== uids_oi[pos]))
  end
end

dict_of_movements_out


###########
### Update outside movements
###########

for t in 1:360
  for pos in 1:size(uids_oi, 1)

    rom = record_of_movements_oi[dict_of_movements_out[(-1, t, pos)], :]

    if (size(rom, 1) > 0)
      println("WAS: ", track_3D_oi[pos, t, 10:12], "AND NOW IS: ", sum.(eachcol(rom[:, 7:9])))

      track_3D_oi[pos, t, 10:12] = sum.(eachcol(rom[:, 7:9]))
    end
  end
end


##########################
### Combine into Array ###
##########################

combi_array = [res_3D_oi, track_3D_oi, pers_3D_oi, parish_3D_oi]

record_of_movements_oi

f_to_p_dict

ids_to_pos_dict

dict_of_movements

####################
### Save Results ###
####################

save("Data/Set 1/res_3D_oi.jld2", "array", res_3D_oi)

save("Data/Set 1/combi_array.jld2", "array", combi_array)

save("Data/Set 1/record_of_movements_oi.jld2", "array", record_of_movements_oi)


save("Data/Set 1/f_to_p_dict.jld2", "dict", f_to_p_dict)

save("Data/Set 1/ids_to_pos_dict.jld2", "dict", ids_to_pos_dict)

save("Data/Set 1/dict_of_movements.jld2", "dict", dict_of_movements)
