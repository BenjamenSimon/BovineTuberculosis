
include("LoadPackages.jl")

#################
### LOAD DATA ###
#################

real_data_results = load("Data/Set Real/real_data_results.jld2")["array"]

real_data_track = load("Data/Set Real/real_data_track.jld2")["array"]

real_data_pers = load("Data/Set Real/real_data_pers.jld2")["array"]

real_data_parish = load("Data/Set Real/real_data_parish.jld2")["array"]

df_chesh_movements = load("Data/Set Real/df_chesh_movements.jld2")["array"]

dict_of_movements = load("Data/Set Real/dict_of_movements.jld2")["array"]

f_to_p_dict = load("Data/Set Real/f_to_p_dict.jld2")["array"]

individual_movements_record = load("Data/Set Real/individual_movements_record.jld2")["array"]

sum(real_data_results[:,1,4:6])

###################################
### Subset and Change Positions ###
###################################

#############
### Choose Parishes
#############

parish_ids_oi = real_data_parish[:, 1, 16] .∈ ([139,354,336,367,349,360,345,347,357,335,344,348,363,404,358,333,359,351,332,342,343,356,327,329], )

first_h_of_p_oi = Array(real_data_parish[parish_ids_oi, 1, 17])

holding_ids_oi = Int64[]

for row_id in first_h_of_p_oi
    for entry in f_to_p_dict[ row_id ][3]
        push!(holding_ids_oi, entry)
    end
end

holding_ids_oi


## Results

real_data_results_subset = deepcopy(real_data_results[holding_ids_oi, :, :])

for pos in 1:307
    real_data_results_subset[pos, :, 3] .= fill(pos, 360)
end

## Track

real_data_track_subset = deepcopy(real_data_track[holding_ids_oi, :, :])

for pos in 1:307
    real_data_track_subset[pos, :, 3] .= fill(pos, 360)
end

## Persistents

real_data_pers_subset = deepcopy(real_data_pers[holding_ids_oi, :, :])

for pos in 1:307
    real_data_pers_subset[pos, :, 3] .= fill(pos, 360)
end

## Parishes

real_data_parish_subset = deepcopy(real_data_parish[parish_ids_oi, :, :])

for parish_pos in 1:24
    real_data_parish_subset[parish_pos, :, 3] .= fill(parish_pos, 360)
end


##########################
### Functional Objects ###
##########################

## Farm to Parish Dictionary

f_to_p_dict_subset = Dict{Integer, Array}()

for pos in 1:307

    row_id = real_data_results_subset[pos, 1, 2]

    farm_pos = real_data_results_subset[pos, 1, 3]

    parish_num = real_data_results_subset[pos, 1, 23]

    parish_uid, parish_pos = real_data_parish_subset[(real_data_parish_subset[:, 1, 16] .== parish_num), 1, 2:3]

    parish_members_uid = Array(real_data_results_subset[(real_data_results_subset[:, 1, 23] .== parish_num), 1, 2])

    parish_members_pos = Array(real_data_results_subset[(real_data_results_subset[:, 1, 23] .== parish_num), 1, 3])

    f_to_p_dict_subset[pos] = [Int64.(parish_uid), Int64.(parish_pos), Int64.(parish_members_uid), Int64.(parish_members_pos), Int64.(row_id), Int64.(farm_pos)]

end

f_to_p_dict_subset

struct Farm_Parish_info
  farm_UID::Int64
  farm_position::Int64
  parish_UID::Int64
  parish_position::Int64
  parish_members_UIDs::Vector{Int64}
  parish_members_positions::Vector{Int64}
end

f_to_p_structs = Array{Farm_Parish_info}(undef, length(f_to_p_dict))

for i in 1:length(f_to_p_dict_subset)

  f_to_p_structs[i] = Farm_Parish_info(f_to_p_dict_subset[i][5],
                                       f_to_p_dict_subset[i][6],
                                       f_to_p_dict_subset[i][1],
                                       f_to_p_dict_subset[i][2],
                                       f_to_p_dict_subset[i][3],
                                       f_to_p_dict_subset[i][4])
end

f_to_p_structs


## farm_uid to position Dictionary

ids_to_pos_dict = Dict{Int64, Int64}()

for pos in 1:307
    ids_to_pos_dict[real_data_results_subset[pos, 1, 2]] = real_data_results_subset[pos, 1, 3]
end

ids_to_pos_dict


#################
### Movements ###
#################

#############
### Subset the Array of Movements
#############

holding_ids_oi

individual_movements_record

record_of_movements_subset = deepcopy(Int64.(individual_movements_record))

for i in 1:size(record_of_movements_subset, 1)
  if ( (record_of_movements_subset[i,2] in holding_ids_oi) == false )
     record_of_movements_subset[i,2] = -1
  end
  if ((record_of_movements_subset[i,3] in holding_ids_oi) == false)
     record_of_movements_subset[i,3] = -1
  end
end

for t in 1:360
  # find all -1 to -1 movements and distinguish them
  indicies = findall((record_of_movements_subset[:,1] .== t) .& (record_of_movements_subset[:,2] .== -1) .& (record_of_movements_subset[:,3] .== -1))

  record_of_movements_subset[indicies,2] .= -99
  record_of_movements_subset[indicies,3] .= -99
end

record_of_movements_subset

record_of_movements_oi = record_of_movements_subset[(record_of_movements_subset[:,2] .> -10), :]

###########
### Movement dictionary
###########

dict_of_movements_oi = Dict{Tuple{Int64, Int64}, Vector{Int64}}()

for pos in 1:307
  for t in 1:360

    off_row_id = real_data_results_subset[pos, 1, 2]

    dict_of_movements_oi[(pos, t)] = findall((record_of_movements_oi[:,1] .== t) .& (record_of_movements_oi[:,2] .== off_row_id))
  end
end

dict_of_movements_oi

dict_of_movements_oi[(27, 1)]

###########
### Update movements from outside in tracking array
###########

dict_of_movements_out = Dict{Tuple{Int64, Int64}, Vector{Int64}}()

for pos in 1:307
  for t in 1:360

    on_row_id = real_data_results_subset[pos, 1, 2]

    dict_of_movements_out[(pos, t)] = findall((record_of_movements_oi[:,1] .== t) .& (record_of_movements_oi[:,2] .== -1) .& (record_of_movements_oi[:,3] .== on_row_id))
  end
end

for t in 1:360
  for pos in 1:307

    on_row_id = real_data_results_subset[pos, 1, 2]

    rom = record_of_movements_oi[dict_of_movements_out[(pos, t)], :]

    if (size(rom, 1) > 0)
      println(pos, "     ", t, "     ", on_row_id, "      ")
      println(rom)
      println("WAS: ", real_data_track_subset[pos, t, 10:12])
      println(" AND NOW IS: ", real_data_track_subset[pos, t, 10:12] + sum.(eachcol(rom[:, 7:9])))

      real_data_track_subset[pos, t, 10:12] += sum.(eachcol(rom[:, 7:9]))
      real_data_track_subset[pos, t, 4:6] -= sum.(eachcol(rom[:, 7:9]))
    end

  end
end

#############
### Subset the Array of Movements
#############

# holding_ids_oi
#
# df_chesh_movements
#
# record_of_movements_subset = deepcopy(Int64.(df_chesh_movements))
#
# for i in 1:size(record_of_movements_subset, 1)
#   if ( (record_of_movements_subset[i,11] in holding_ids_oi) == false )
#      record_of_movements_subset[i,[2,3,4,5,11]] = [-1,-1,-1,-1,-1]
#   end
#   if ((record_of_movements_subset[i,12] in holding_ids_oi) == false)
#      record_of_movements_subset[i,[6,7,8,9,12]] = [-1,-1,-1,-1,-1]
#   end
# end
#
# for t in 1:360
#   # find all -1 to -1 movements and distinguish them
#   indicies = findall((record_of_movements_subset[:,1] .== t) .& (record_of_movements_subset[:,2] .== -1) .& (record_of_movements_subset[:,6] .== -1))
#
#   record_of_movements_subset[indicies,2] .= -99
#   record_of_movements_subset[indicies,6] .= -99
# end
#
# record_of_movements_subset
#
# record_of_movements_subset = record_of_movements_subset[(record_of_movements_subset[:,2] .> -10), :]
#
# record_of_movements_oi = combine(groupby(record_of_movements_subset, [:week_no, :first_off_cph, :first_off_county, :first_off_parish, :first_off_holding, :last_on_cph, :last_on_county, :last_on_parish, :last_on_holding, :first_off_row_id, :last_on_row_id]), :n_moves_off => sum)
#
# record_of_movements_oi = select(record_of_movements_oi, :week_no, :first_off_cph, :first_off_county, :first_off_parish, :first_off_holding, :last_on_cph, :last_on_county, :last_on_parish, :last_on_holding, :n_moves_off_sum => :n_moves_off, :first_off_row_id, :last_on_row_id)

###########
### Movement dictionary
###########

# dict_of_movements_oi = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
#
# for off_row_id in unique(record_of_movements_oi.first_off_row_id)
#   for t in 1:360
#     dict_of_movements_oi[(off_row_id, t)] = findall((record_of_movements_oi[:,1] .== t) .& (record_of_movements_oi[:,11] .== off_row_id))
#   end
#   println("off_row_id = ", off_row_id)
# end


##########################
### Combine into Array ###
##########################

real_data_results_subset

real_data_track_subset

real_data_pers_subset

real_data_parish_subset


array1 = Array{Int64, 3}(real_data_results_subset)
array2 = Array{Int64, 3}(real_data_track_subset)
array3 = Array{Float64, 3}(real_data_pers_subset)
array4 = Array{Float64, 3}(real_data_parish_subset)


DATA_res_and_track = [array1, array2]
DATA_pers_and_parish = [array3, array4]



f_to_p_structs

ids_to_pos_dict

record_of_movements_oi

dict_of_movements_oi




####################
### Save Results ###
####################

save("Data/Set Real Subset/DATA_res_and_track.jld2", "array", DATA_res_and_track)
save("Data/Set Real Subset/DATA_pers_and_parish.jld2", "array", DATA_pers_and_parish)

save("Data/Set Real Subset/record_of_movements_oi.jld2", "array", record_of_movements_oi)


save("Data/Set Real Subset/f_to_p_structs.jld2", "struct", f_to_p_structs)

save("Data/Set Real Subset/ids_to_pos_dict.jld2", "dict", ids_to_pos_dict)

save("Data/Set Real Subset/dict_of_movements_oi_real.jld2", "dict", dict_of_movements_oi)
