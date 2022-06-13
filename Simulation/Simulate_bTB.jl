
include("../LoadPackages.jl")

##############################
### Load the Cheshire Data ###
##############################

cph_to_formatted = load("Data/cph_to_formatted_chesh.rds")

df_chesh_initial_states = load("Data/df_chesh_initial_states.rds")

df_chesh_birth_and_death = load("Data/df_chesh_birth_and_death.rds")

df_chesh_movements_final = load("Data/df_chesh_movements_final.rds")

##################################################
### Number of counties, parishes, and holdings ###
##################################################


n_cph = cph_to_formatted[:, [:n_county_total, :n_parish_in_county, :n_holding_in_county_parish, :county_idx, :parish_idx, :holding_idx]]


n_counties = unique(n_cph[:, 1])[1]

n_parishes = fill(0, n_counties)
for i in 1:n_counties
  county = @where(n_cph, :county_idx .== i)
  global n_parishes[i]  = (county[1, :n_parish_in_county])[1]
end

n_holdings = fill([], n_counties)
for i in 1:n_counties
  global n_holdings[i] = fill([], n_parishes[i])
  for j in 1:n_parishes[i]
    parish = @where(n_cph, :county_idx .== i, :parish_idx .== j)
    global n_holdings[i][j] = (parish[1, :n_holding_in_county_parish])[1]
  end
end

area_of_parish = fill([], n_counties)
for i in 1:n_counties
  global area_of_parish[i] = fill([], n_parishes[i])
  for j in 1:n_parishes[i]
    area = n_holdings[i][j] * 1000
    global area_of_parish[i][j] = area
  end
end

total_holdings = cph_to_formatted[:, :n_holding_total][1]

###################################
### Functions to create holding ###
###################################

# If you put an argument after ; it becomes a keyword argument
# If you then don't give it a default value, it must be named in the function call
# If you don't include the ; it makes it optional, and can be named

function create_holding_res(;Time::Int64, row_id::Int64,
                             cStates_init::Array{Int64, 1}, bStates_init::Array{Int64, 1},
                             cph::Array{Int64, 1}, cph_idx::Array{Int64, 1})

  holding_res = fill(0., Time, 49)


  # [:t, :row_id,
  # :cS_init, :cE_init, :cI_init,
  # :cS_Moves, :cE_Moves, :cI_Moves,
  # :cS_postM, :cE_postM, :cI_postM,
  # :cS_postEI, :cE_postEI, :cI_postEI,
  # :cS_postDet, :cE_postDet, :cI_postDet,
  # :cS_final, :cE_final, :cI_final,
  # :bS_init, :bE_init, :bI_init,
  # :bS_postEI, :bE_postEI, :bI_postEI,
  # :bS_final, :bE_final, :bI_final,
  # :move_res_init, :move_res_final,
  # :pcS_init, :pcE_init, :pcI_init,
  # :pcS_final, :pcE_final, :pcI_final,
  # :pbS_init, :pbE_init, :pbI_init,
  # :pbS_final, :pbE_final, :pbI_final,
  # :county, :parish, :holding,
  # :county_idx, :parish_idx, :holding_idx]

  holding_res[1, :] .= [1.; row_id;
                            cStates_init; fill(0., 15); bStates_init; fill(0., 6);
                            false; false;
                            fill(0., 12);
                            cph; cph_idx]

  return(holding_res)
end

# create_holding_res( Time = 360, row_id = 9,
#                     cStates_init = [100, 0, 19], bStates_init = [50, 0, 5],
#                     cph = [2, 3, 7], cph_idx = [1,2,4])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_holding_track(;Time::Int64, row_id::Int64,
                             next_test_week::Int64,
                             cph::Array{Int64, 1}, cph_idx::Array{Int64, 1})

   # [:t, :row_id,
   # :sus_on, :exp_on, :inf_on,
   # :sus_off, :exp_off, :inf_off,
   # :S_on_out, :E_on_out, :I_on_out,
   # :c_new_exp, :c_new_inf, :b_new_exp, :b_new_inf,
   # :test_occur, :next_test_week,
   # :E_detected, :I_detected,
   # :c_birth, :c_S_death, :c_E_death, :c_I_death,
   # :b_birth, :b_S_death, :b_E_death, :b_I_death,
   # :remaining_pressure, :new_pressure
   # :county, :parish, :holding,
   # :county_idx, :parish_idx, :holding_idx]

  holding_track = fill(0., Time, 35)
  holding_track[1,:] .= [1.; row_id; fill(0., 14); next_test_week; fill(0., 12); cph; cph_idx]

  return(holding_track)
end

# create_holding_track(Time = 360, row_id = 9, next_test_week = 5,
#                     cph = [2, 3, 7], cph_idx = [1,2,4])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_holding_pers(;Time::Int64, row_id::Int64,
                             cph::Array{Int64, 1}, cph_idx::Array{Int64, 1})

   # [:t, :row_id,
   # :p_env_prev, :p_env_cur,
   # :c_exp_prob, :b_exp_prob,
   # :c_inf_prob, :b_inf_prob,
   # :county, :parish, :holding,
   # :county_idx, :parish_idx, :holding_idx]

  holding_pers = fill(0., Time, 14)
  holding_pers[1,:] .= [1.; row_id; fill(0., 6); cph; cph_idx]

  return(holding_pers)
end

# create_holding_pers(Time = 360, row_id = 9,
#                     cph = [2, 3, 7], cph_idx = [1,2,4])


##################################
### Functions to create parish ###
##################################

function create_parish_res(;Time::Int64, initial_settings::DataFrame,
                            cp::Array{Int64, 1}, cp_idx::Array{Int64, 1})

  # Extract the dataframe for this parish
  parish_df = @where(initial_settings, :county .== cp[1], :parish .== cp[2])

  # Extract a list of holdings and their ids from that parish
  row_ids = unique(convert(Array{Int64}, parish_df[:, :row_id]))
  holdings = unique(convert(Array{Int64}, parish_df[:, :holding]))
  holdings_idx = unique(convert(Array{Int64}, parish_df[:, :holding_idx]))

  # Calculate how many holdings there are
  num_holdings = length(unique(holdings))

  # Extract the initial states of the cattle
  cN_init = convert(Array{Int64}, parish_df[:, :num_animals])
  cI_init = convert(Array{Int64}, parish_df[:, :estimated_inf])
  cS_init = cN_init - cI_init

  ##############
  ### Inputs ###
  ##############

  # Parish results
  parish_res = Array{Array{Int64, 2}}(undef, 1, num_holdings)

  for i = 1:num_holdings

    #println("row_ids", row_ids[i])

    cStates_i = [cS_init[i], 0, cI_init[i]]

    prop = cStates_i[3]/cN_init[i]
    prop = if isnan(prop) 0 else prop end

    bS_init = sample(10:100)
    bI_init = rand(Binomial(bS_init, prop))

    bStates_i = convert(Array{Int64}, [bS_init, 0., bI_init])

    # Time (::Int64)
    # row_id (::Int64)
    # cStates_init (Array = [cS, cE, cI])
    # bStates_init (Array = [bS, bE, bI])
    # cph (Array = [c, p, h])
    # cph_idx (Array = [c_idx, p_idx, h_idx])

    parish_res[1, i] = create_holding_res(
      Time = Time,
      row_id = row_ids[i],
      cStates_init = cStates_i,
      bStates_init = bStates_i,
      cph = [cp; holdings[i]],
      cph_idx = [cp_idx; holdings_idx[i]]
    )
  end

  return(parish_res)
end

# temp_parish_res = create_parish_res(Time = 360, initial_settings = df_chesh_initial_states,
#                     cp = [6, 1], cp_idx = [1, 1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_parish_track(;Time::Int64, initial_settings::DataFrame,
                            cp::Array{Int64, 1}, cp_idx::Array{Int64, 1})

  # Extract the dataframe for this parish
  parish_df = @where(initial_settings, :county .== cp[1], :parish .== cp[2])

  # Extract a list of holdings and their ids from that parish
  row_ids = unique(convert(Array{Int64}, parish_df[:, :row_id]))
  holdings = unique(convert(Array{Int64}, parish_df[:, :holding]))
  holdings_idx = unique(convert(Array{Int64}, parish_df[:, :holding_idx]))

  # Calculate how many holdings there are
  num_holdings = length(unique(holdings))

  # Extract the initial next test week
  test_weeks = convert(Array{Int64}, parish_df[:, :test_week])


  ##############
  ### Inputs ###
  ##############

  # Parish track
  parish_track = Array{Array{Int64, 2}}(undef, 1, num_holdings)

  # Parish
  for i = 1:num_holdings
    parish_track[1, i] = create_holding_track(Time = Time,
                                              row_id = row_ids[i],
                                              next_test_week = test_weeks[i],
                                              cph = [cp; holdings[i]],
                                              cph_idx = [cp_idx; holdings_idx[i]])
  end

  return(parish_track)
end

# temp_parish_track = create_parish_track(Time = 360, initial_settings = df_chesh_initial_states,
#                                         cp = [6, 1], cp_idx = [1, 1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_parish_pers(;Time::Int64, initial_settings::DataFrame,
                            cp::Array{Int64, 1}, cp_idx::Array{Int64, 1})

  # Extract the dataframe for this parish
  parish_df = @where(initial_settings, :county .== cp[1], :parish .== cp[2])

  # Extract a list of holdings and their ids from that parish
  row_ids = unique(convert(Array{Int64}, parish_df[:, :row_id]))
  holdings = unique(convert(Array{Int64}, parish_df[:, :holding]))
  holdings_idx = unique(convert(Array{Int64}, parish_df[:, :holding_idx]))

  # Calculate how many holdings there are
  num_holdings = length(unique(holdings))


  ##############
  ### Inputs ###
  ##############

  # Parish track
  parish_pers = Array{Array{Float64, 2}}(undef, 1, num_holdings)

  # Parish
  for i = 1:num_holdings
    parish_pers[1, i] = create_holding_pers(Time = Time,
                                              row_id = row_ids[i],
                                              cph = [cp; holdings[i]],
                                              cph_idx = [cp_idx; holdings_idx[i]])
  end

  return(parish_pers)
end

# temp_parish_pers = create_parish_pers(Time = 360, initial_settings = df_chesh_initial_states,
#                                       cp = [6, 1], cp_idx = [1, 1])


##################################
### Functions to create county ###
##################################

function create_county_res(;Time::Int64, initial_settings::DataFrame,
                            c::Int64, c_idx::Int64)

  # Extract the dataframe for this county
  county_df = @where(initial_settings, :county .== c)

  # Extract a list of parishes and their ids from that county
  parishes = unique(convert(Array{Int64}, county_df[:, :parish]))
  parishes_idx = unique(convert(Array{Int64}, county_df[:, :parish_idx]))

  # Calculate how many parishes there are
  num_parishes = length(unique(parishes))

  ##############
  ### Inputs ###
  ##############

  # County results
  county_res = Array{Array{Array{Int64,2}, 2}}(undef, 1, num_parishes)

  # Time (::Int64)
  # initial_settings (DataFrame)
  # cp (Array = [c, p])
  # cp_idx (Array = [c_idx, p_idx])

  for i = 1:num_parishes

    county_res[1, i] = create_parish_res(Time = Time,
                                          initial_settings = initial_settings,
                                          cp = [c; parishes[i]],
                                          cp_idx = [c_idx; parishes_idx[i]])
  end
  return(county_res)
end

# temp_county_res = create_county_res(Time = 360, initial_settings = df_chesh_initial_states,
#                     c = 6, c_idx = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_county_track(;Time::Int64, initial_settings::DataFrame,
                            c::Int64, c_idx::Int64)

  # Extract the dataframe for this county
  county_df = @where(initial_settings, :county .== c)

  # Extract a list of parishes and their ids from that county
  parishes = unique(convert(Array{Int64}, county_df[:, :parish]))
  parishes_idx = unique(convert(Array{Int64}, county_df[:, :parish_idx]))

  # Calculate how many parishes there are
  num_parishes = length(unique(parishes))

  ##############
  ### Inputs ###
  ##############

  # County track
  county_track = Array{Array{Array{Int64,2}, 2}}(undef, 1, num_parishes)

  # Time (::Int64)
  # initial_settings (DataFrame)
  # cp (Array = [c, p])
  # cp_idx (Array = [c_idx, p_idx])

  for i = 1:num_parishes
    county_track[1, i] = create_parish_track(Time = Time,
                                            initial_settings = initial_settings,
                                            cp = [c; parishes[i]],
                                            cp_idx = [c_idx; parishes_idx[i]])
  end

  return(county_track)
end

# temp_county_track = create_county_track(Time = 360, initial_settings = df_chesh_initial_states,
#                                         c = 6, c_idx = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_county_pers(;Time::Int64, initial_settings::DataFrame,
                            c::Int64, c_idx::Int64)

  # Extract the dataframe for this county
  county_df = @where(initial_settings, :county .== c)

  # Extract a list of parishes and their ids from that county
  parishes = unique(convert(Array{Int64}, county_df[:, :parish]))
  parishes_idx = unique(convert(Array{Int64}, county_df[:, :parish_idx]))

  # Calculate how many parishes there are
  num_parishes = length(unique(parishes))

  ##############
  ### Inputs ###
  ##############

  # County track
  county_pers = Array{Array{Array{Float64,2}, 2}}(undef, 1, num_parishes)

  # Time (::Int64)
  # initial_settings (DataFrame)
  # cp (Array = [c, p])
  # cp_idx (Array = [c_idx, p_idx])

  for i = 1:num_parishes
    county_pers[1, i] = create_parish_pers(Time = Time,
                                            initial_settings = initial_settings,
                                            cp = [c; parishes[i]],
                                            cp_idx = [c_idx; parishes_idx[i]])
  end

  return(county_pers)
end

# temp_county_pers = create_county_pers(Time = 360, initial_settings = df_chesh_initial_states,
#                                         c = 6, c_idx = 1)



###################################
### Functions to create country ###
###################################

function create_country_res(;Time::Int64, initial_settings::DataFrame)

  # Extract a list of counties and their ids
  counties = unique(convert(Array{Int64}, initial_settings[:, :county]))
  counties_idx = unique(convert(Array{Int64}, initial_settings[:, :county_idx]))

  # Calculate how many counties there are
  num_counties = length(unique(counties))

  ##############
  ### Inputs ###
  ##############

  # Country results
  country_res = Array{Array{Array{Array{Int64,2}, 2}, 2}}(undef, 1, num_counties)

  # Time (::Int64)
  # initial_settings (DataFrame)
  # c (::Int64)
  # c_idx (::Int64)

  for i = 1:num_counties
    country_res[1, i] = create_county_res(Time = Time,
                                          initial_settings = initial_settings,
                                          c = counties[i],
                                          c_idx = counties_idx[i])
  end
  return(country_res)
end

# temp_country_res = create_country_res(Time = 360, initial_settings = df_chesh_initial_states)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_country_track(;Time::Int64, initial_settings::DataFrame)

  # Extract a list of counties and their ids
  counties = unique(convert(Array{Int64}, initial_settings[:, :county]))
  counties_idx = unique(convert(Array{Int64}, initial_settings[:, :county_idx]))

  # Calculate how many counties there are
  num_counties = length(unique(counties))

  ##############
  ### Inputs ###
  ##############

  # Country results
  country_track = Array{Array{Array{Array{Int64,2}, 2}, 2}}(undef, 1, num_counties)

  for i = 1:num_counties
    country_track[1, i] = create_county_track(Time = Time,
                                              initial_settings = initial_settings,
                                              c = counties[i],
                                              c_idx = counties_idx[i])
  end
  return(country_track)
end

# temp_country_track = create_country_track(Time = 360, initial_settings = df_chesh_initial_states)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_country_pers(;Time::Int64, initial_settings::DataFrame)

  # Extract a list of counties and their ids
  counties = unique(convert(Array{Int64}, initial_settings[:, :county]))
  counties_idx = unique(convert(Array{Int64}, initial_settings[:, :county_idx]))

  # Calculate how many counties there are
  num_counties = length(unique(counties))

  ##############
  ### Inputs ###
  ##############

  # Country results
  country_pers = Array{Array{Array{Array{Float64,2}, 2}, 2}}(undef, 1, num_counties)

  for i = 1:num_counties
    country_pers[1, i] = create_county_pers(Time = Time,
                                              initial_settings = initial_settings,
                                              c = counties[i],
                                              c_idx = counties_idx[i])
  end
  return(country_pers)
end

# temp_country_pers = create_country_pers(Time = 360, initial_settings = df_chesh_initial_states)




############################################
### Extract the states at a given time t ###
############################################

function extract_states(;t, country_res, n_counties, n_parishes, n_holdings, total_holdings)

  n_cols = size(country_res[1][1][1], 2)

  states = fill(0., total_holdings, n_cols)

  #[:t, :row_id,
  # :cS_init, :cE_init, :cI_init,
  # :cS_postM, :cE_postM, :cI_postM,
  # :cS_postEI, :cE_postEI, :cI_postEI,
  # :cS_postDet, :cE_postDet, :cI_postDet,
  # :cS_final, :cE_final, :cI_final,
  # :bS_init, :bE_init, :bI_init,
  # :bS_postEI, :bE_postEI, :bI_postEI,
  # :bS_final, :bE_final, :bI_final,
  # :move_res_init, :move_res_final,
  # :county, :parish, :holding,
  # :county_idx, :parish_idx, :holding_idx]

  rowid = 1
  for h = 1:n_counties
    for i = 1:n_parishes[h]
      for j = 1:n_holdings[h][i]
        for k = 1:n_cols
          begin
            states[rowid, k] = country_res[h][i][j][t, k]
            states
          end # end of begin
        end # end of column loop
        rowid = rowid + 1
      end # end of holdings loop
    end # end of parishes loop
  end # end of counties loop
  return (states)
end

# temp_current_states = extract_states(t = 1, country_res = temp_country_res,
#                                      n_counties = n_counties, n_parishes = n_parishes,
#                                      n_holdings = n_holdings, total_holdings = total_holdings)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#


function extract_track(;t, country_track, n_counties, n_parishes, n_holdings, total_holdings)

  n_cols = size(country_track[1][1][1], 2)

  states = fill(0., total_holdings, n_cols)

  #[:t, :row_id,
  # :sus_on, :exp_on, :inf_on,
  # :sus_off, :exp_off, :inf_off,
  # :S_on_out, :E_on_out, :I_on_out,
  # :c_new_exp, :c_new_inf, :b_new_exp, :b_new_inf,
  # :test_occur, next_test_week,
  # :E_detected, :I_detected,
  # :c_birth, :c_S_death, :c_E_death, :c_I_death,
  # :b_birth, :b_S_death, :b_E_death, :b_I_death,
  # :county, :parish, :holding,
  # :county_idx, :parish_idx, :holding_idx]

  rowid = 1
  for h = 1:n_counties
    for i = 1:n_parishes[h]
      for j = 1:n_holdings[h][i]
        for k = 1:n_cols
          begin
            states[rowid, k] = country_track[h][i][j][t, k]
            states
          end # end of begin
        end # end of column loop
        rowid = rowid + 1
      end # end of holdings loop
    end # end of parishes loop
  end # end of counties loop
  return (states)
end

# temp_current_track = extract_track(t = 1, country_track = temp_country_track,
#                                      n_counties = n_counties, n_parishes = n_parishes,
#                                      n_holdings = n_holdings, total_holdings = total_holdings)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function extract_pers(;t, country_pers, n_counties, n_parishes, n_holdings, total_holdings)

  n_cols = size(country_pers[1][1][1], 2)

  states = fill(0., total_holdings, n_cols)

  #[:t, :row_id,
  # :p_env_prev, :p_env_cur,
  # c_exp_prob, :b_exp_prob,
  # c_inf_prob, :b_inf_prob,
  # :county, :parish, :holding,
  # :county_idx, :parish_idx, :holding_idx]

  rowid = 1
  for h = 1:n_counties
   for i = 1:n_parishes[h]
     for j = 1:n_holdings[h][i]
       for k = 1:n_cols
         begin
           states[rowid, k] = country_pers[h][i][j][t, k]
           states
         end # end of begin
       end # end of column loop
       rowid = rowid + 1
     end # end of holdings loop
   end # end of parishes loop
  end # end of counties loop
  return (states)
end

# temp_current_pers = extract_pers(t = 1, country_pers = temp_country_pers,
#                                     n_counties = n_counties, n_parishes = n_parishes,
#                                     n_holdings = n_holdings, total_holdings = total_holdings)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function extract_initial_state(;country_res, n_counties, n_parishes, n_holdings, total_holdings)

  temp_current_states = extract_states(t = 1, country_res = country_res,
                                       n_counties = n_counties, n_parishes = n_parishes,
                                       n_holdings = n_holdings, total_holdings = total_holdings)

  ###########################
  ### Set the initial, intermediate, and final states of animals to be identical
  ###########################

  # :cS_postM, :cE_postM, :cI_postM = :cS_init, :cE_init, :cI_init
  temp_current_states[:, 9:11] .= temp_current_states[:, 3:5]
  # :cS_postEI, :cE_postEI, :cI_postEI = :cS_init, :cE_init, :cI_init
  temp_current_states[:, 12:14] .= temp_current_states[:, 3:5]
  # :cS_postDet, :cE_postDet, :cI_postDet = :cS_init, :cE_init, :cI_init
  temp_current_states[:, 15:17] .= temp_current_states[:, 3:5]
  # :cS_final, :cE_final, :cI_final = :cS_init, :cE_init, :cI_init
  temp_current_states[:, 18:20] .= temp_current_states[:, 3:5]

  # :bS_postEI, :bE_postEI, :bI_postEI = :bS_init, :bE_init, :bI_init
  temp_current_states[:, 24:26] .= temp_current_states[:, 21:23]
  # :bS_final, :bE_final, :bI_final = :bS_init, :bE_init, :bI_init
  temp_current_states[:, 27:29] .= temp_current_states[:, 21:23]


  ###########################
  ### Calculate the total animals in the parishes
  ###########################

  for c in 1:n_counties
    for p in 1:n_parishes[c]

      # Find the members of each parish
      parish_mems = findall((temp_current_states[:,47] .== c) .& (temp_current_states[:,48] .== p))

      # Calcualte the total animals of each state in the parish
      # at the end of the time step
      pcS = sum(temp_current_states[parish_mems, 18])
      pcE = sum(temp_current_states[parish_mems, 19])
      pcI = sum(temp_current_states[parish_mems, 20])

      pbS = sum(temp_current_states[parish_mems, 27])
      pbE = sum(temp_current_states[parish_mems, 28])
      pbI = sum(temp_current_states[parish_mems, 29])

      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      temp_current_states[parish_mems, [32,35]] .= pcS
      temp_current_states[parish_mems, [33,36]] .= pcE
      temp_current_states[parish_mems, [34,37]] .= pcI

      temp_current_states[parish_mems, [38,41]] .= pbS
      temp_current_states[parish_mems, [39,42]] .= pbE
      temp_current_states[parish_mems, [40,43]] .= pbI

    end
  end

  return(temp_current_states)
end

# temp_initial_states = extract_initial_state(country_res = temp_country_res,
#                                      n_counties = n_counties, n_parishes = n_parishes,
#                                      n_holdings = n_holdings, total_holdings = total_holdings)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function extract_initial_track(;country_track, n_counties, n_parishes, n_holdings, total_holdings)

 temp_current_track = extract_track(t = 1, country_track = country_track,
                                      n_counties = n_counties, n_parishes = n_parishes,
                                      n_holdings = n_holdings, total_holdings = total_holdings)

 return(temp_current_track)
end

# temp_initial_track = extract_initial_track(country_track = temp_country_track,
#                                     n_counties = n_counties, n_parishes = n_parishes,
#                                     n_holdings = n_holdings, total_holdings = total_holdings)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function extract_initial_pers(;country_pers, n_counties, n_parishes, n_holdings, total_holdings)

 temp_current_pers = extract_pers(t = 1, country_pers = country_pers,
                                      n_counties = n_counties, n_parishes = n_parishes,
                                      n_holdings = n_holdings, total_holdings = total_holdings)

 return(temp_current_pers)
end

# temp_initial_pers = extract_initial_pers(country_pers = temp_country_pers,
#                                     n_counties = n_counties, n_parishes = n_parishes,
#                                     n_holdings = n_holdings, total_holdings = total_holdings)





#################################################################
### Functions for calculating and recording parish level data ###
#################################################################

function add_parish_final_totals(;current_state, n_counties, n_parishes)

  for c in 1:n_counties
    for p in 1:n_parishes[c]

      # Find the members of each parish
      parish_mems = findall((current_state[:,47] .== c) .& (current_state[:,48] .== p))

      # Calcualte the total animals of each state in the parish
      # at the end of the time step
      pcS_final = sum(current_state[parish_mems, 18])
      pcE_final = sum(current_state[parish_mems, 19])
      pcI_final = sum(current_state[parish_mems, 20])

      pbS_final = sum(current_state[parish_mems, 27])
      pbE_final = sum(current_state[parish_mems, 28])
      pbI_final = sum(current_state[parish_mems, 29])

      # Record for each holding in the parish, the total numbers
      # of animals in the parish
      current_state[parish_mems, 35] .= pcS_final
      current_state[parish_mems, 36] .= pcE_final
      current_state[parish_mems, 37] .= pcI_final

      current_state[parish_mems, 41] .= pbS_final
      current_state[parish_mems, 42] .= pbE_final
      current_state[parish_mems, 43] .= pbI_final

    end
  end

  return(current_state)
end

# temp_initial_states = add_parish_final_totals(current_state = temp_initial_states,
#                                               n_counties = n_counties, n_parishes = n_parishes)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#


function add_parish_env(;current_state, current_track, current_pers, n_counties, n_parishes, epi_pars, scaling)

  # Extract the decay rate
  ϵ = epi_pars[5]

  for c in 1:n_counties
    for p in 1:n_parishes[c]

      # Find the members of each parish
      parish_mems = findall((current_state[:,47] .== c) .& (current_state[:,48] .== p))

      # Extract the previous parish environmental level
      p_env_prev = current_pers[parish_mems[1], 3]

      # Extract the parish level states of all animals
      pI_init = sum(current_state[parish_mems[1], [34,40]])

      # Calculate new parish environmental level

      new_pressure = rand(Poisson(pI_init))

      remaining_pressure = rand(Binomial(round(Int, scaling[c][p] * p_env_prev), (1-ϵ)))

      parish_env = (remaining_pressure + new_pressure) / scaling[c][p]



      # Record for each holding in the parish, the parish level
      # environmental reservoir for the current timestep
      current_pers[parish_mems, 4] .= parish_env

      current_track[parish_mems, 28] .= remaining_pressure
      current_track[parish_mems, 29] .= new_pressure

    end
  end

  return(current_track, current_pers)
end

# temp_current_pers = add_parish_env(current_state = temp_initial_states,
#                                               current_pers = temp_current_pers,
#                                               n_counties = n_counties, n_parishes = n_parishes,
#                                               epi_pars = [0,0,0,0,0,0.05,0,0])


#######################################
### Functions to generate movements ###
#######################################

function rng_mvhyper(n, k)
  # n is a vector of size of pop for each group
  # k is the number of trials (total number moving)

  N = sum(n) # total pop size
  m = length(n) # number of groups
  n_otr = N - n[1] # number not in first group

  x = fill(0, m) # results

  x[1] = rand(Hypergeometric(n[1], n_otr, k))

  for i in 2:(m-1)
    n_otr = n_otr - n[i]
    k = k - x[i-1]
    x[i] = rand(Hypergeometric(n[i], n_otr, k))
  end

  x[m] = k - x[m-1]
  return(x)
end

#rng_mvhyper([150, 10, 20], 10)

function add_move_restrictions(moves_data, move_res)

  # Movement restrictions for farms moving animals off

  # moves_data = join(
  #   moves_data,
  #   move_res,
  #   on = [
  #     (:off_county => :county),
  #     (:off_parish => :parish),
  #     (:off_holding => :holding),
  #   ],
  #   kind = :left,
  # )

  moves_data = leftjoin(moves_data, move_res, on = [
                :off_county => :county,
                :off_parish => :parish,
                :off_holding => :holding],
                makeunique=false, indicator=nothing, validate=(false, false))

  rename!(moves_data, :move_res => :off_move_res)


  # Movement restrictions for farms moving animals on

  # moves_data = join(
  #   moves_data,
  #   move_res,
  #   on = [
  #     (:on_county => :county),
  #     (:on_parish => :parish),
  #     (:on_holding => :holding),
  #   ],
  #   kind = :left,
  # )

  moves_data = leftjoin(moves_data, move_res, on = [
                :on_county => :county,
                :on_parish => :parish,
                :on_holding => :holding],
                makeunique=false, indicator=nothing, validate=(false, false))

  rename!(moves_data, :move_res => :on_move_res)

  return(moves_data)
end


function generate_moves(;t, prev_states, hist_moves, unused_movements_MR, unused_movements_F)

  ##############################################
  ### Create an empty movement results table ###
  ##############################################

  #[:t, :row_id,
  # :cS_Moves, :cE_Moves, :cI_Moves, :cN_Moves,
  # :cMove_off, :cS_move_off, :cE_move_off, :cI_move_off,
  # :cS_move_on, :cE_move_on, :cI_move_on,
  # :cS_on_out, :cE_on_out, :cI_on_out
  # :move_res_init,
  # :county, :parish, :holding,
  # :country_idx, :parish_idx, :holding_idx]

  total_farms = size(prev_states, 1)

  movements = fill(0.0, total_farms, 23)


  #[:t, :row_id_off, :row_id_on, :cS_gen, :cE_gen, :cI_gen,
  # :cS_sent, :cE_sent, :cI_sent, cTotal_sent]
  record_of_movements = fill(-10., (total_farms*5), 10)
  rom_tr = 1
  # Would like
  # :off_county, :off_parish, :off_holding,
  # :off_country_idx, :off_parish_idx, :off_holding_idx
  # :on_county, :on_parish, :on_holding,
  # :on_country_idx, :on_parish_idx, :on_holding_idx]

  ##################################################
  ### Fill in elements with their current states ###
  ##################################################

  # :t, :row_id
  movements[:, 1:2] .= prev_states[:, 1:2]
  # :county, :parish, :holding, :country_idx, :parish_idx, :holding_idx
  movements[:, 18:23] .= prev_states[:, 44:49]
  # :cS_Moves, :cE_Moves, :cI_Moves <= :cS_final, :cE_final, :cI_final
  movements[:, 3:5] .= prev_states[:, 18:20]
  # :cN_Moves <= sum(:cS_final, :cE_final, :cI_final)
  movements[:, 6] .= sum.(eachrow(prev_states[:, 18:20]))
  # :move_res(trictions)_init <= :move_res_final
  movements[:, 17] .= prev_states[:, 31]


  # Make a seperate table for tracking movement restrictions
  move_res = DataFrame(
    county = movements[:, 18],
    parish = movements[:, 19],
    holding = movements[:, 20],
    move_res = movements[:, 17],
  )

  ############################################
  ### Extract the movement data for time t ###
  ############################################

  # This extracts the movement data at time t
  moves_data_raw = @where(hist_moves, :week_no .== t)

  # This then augments it with which farms are under movement restrictions
  moves_data_with_MR = add_move_restrictions(moves_data_raw, move_res)

  # Record the movements lost due to movement restrictions
  moves_data_MR_off = @where(moves_data_with_MR, isequal.(:off_move_res, 1))
  moves_data_MR_on = @where(moves_data_with_MR, isequal.(:on_move_res, 1))
  append!(unused_movements_MR, moves_data_MR_off)
  append!(unused_movements_MR, moves_data_MR_on)

  # Remove the movements that won't occur due to movement restrictions
  moves_data_wo_1MR = @where(moves_data_with_MR, isequal.(:off_move_res, 1) .== false)
  moves_data = @where(moves_data_wo_1MR, isequal.(:on_move_res, 1) .== false)

  ####################################################################
  ### Find all the animals coming into cheshire from farms outside ###
  ####################################################################

  # One row for each unique on_cph that recieved animals from outside of cheshire
  on_cphs_all = unique(moves_data[
    :,
    [:off_county, :on_county, :on_parish, :on_holding, :on_row_id],
  ])

  from_out_cphs = @where(on_cphs_all, :off_county .== -1)

  ############################################################
  ### Update farms with movements from outside of cheshire ###
  ############################################################

  ## For each of the farms inside of cheshire that recieved animals from
  ## outside of cheshire:

  for j = 1:size(from_out_cphs, 1)

    # Extract the row id of the j'th farm in cheshire that recieved animals from outside
    on_index = from_out_cphs[j, 5]

    # Extract the movements into that cheshire farm from farms outside of cheshire
    moves_from_out =
      @where(moves_data, isequal.(:on_row_id, on_index), :off_county .== -1)

    # Sum up the total number of animals moved onto that farm
    total_moved = sum(moves_from_out[:, :n_moves_off])

    # Randomly draw the states of the animals moved
    move_states = rng_mvhyper([90000, 5000, 5000], total_moved)

    # Update the total moves on:
    # :cS_move_on, :cE_move_on, :cI_move_on
    movements[on_index, 11:13] .= movements[on_index, 11:13] + move_states

    # Update the total moves on from outside of cheshire
    movements[on_index, 14:16] .= movements[on_index, 14:16] + move_states

    # Record Movement details
    record_of_movements[rom_tr, :] .= [t, -1., on_index, 90000., 5000., 5000.,
                                       move_states[1], move_states[2], move_states[3],
                                       total_moved]
    rom_tr += 1
  end


  ############################################################
  ### Find all the CPHs in cheshire that moved animals off ###
  ############################################################

  # One row for each cph that moved animals off
  off_cphs_all =
    unique(moves_data[:, [:off_county, :off_parish, :off_holding, :off_row_id]])

  # Remove cphs that aren't in cheshire
  off_cphs = @where(off_cphs_all, :off_county .!= -1)

  ###########################################################
  ### Update farms with movements from inside of cheshire ###
  ###########################################################

  # Due to the aggregation and some artefacts in the data set,
  # as well as the introduced stochasticity of the simulation,
  # some movements won't be valid as it will involve moving more
  # animals than currently exist on a farm.
  # The following code loops through all farms, records and skips the
  # ones that have insufficient animals, and then goes back to them on
  # the next loop.


  ## WHILE number of movements to be updated is > 0
  while (size(off_cphs, 1) > 0)

    # Temporary movement array to record the movements that didn't work in the
    # last loop
    new_off_cphs = DataFrame(
      off_county = Float64[],
      off_parish = Float64[],
      off_holding = Float64[],
      off_row_id = Int32[],
    )


    # FOR each of the farms that moved animals off
    for i = 1:size(off_cphs, 1)

      # Extract the row id of the farm moving animals off
      off_index = off_cphs[i, :off_row_id]

      # Extract the movement off data for that farm
      off_moves = @where(moves_data, isequal.(:off_row_id, off_index))

      # Calculate the total number of animals moved off the farm
      total_moved = sum(off_moves[:, :n_moves_off])

      # IF there are more animals moved than currently exist on the farm (:cN_init)
        # Push that movement data to the temporary movement array
        # and skip this movement for now
      if (total_moved > movements[off_index, 6])
        push!(new_off_cphs, off_cphs[i, :])
        continue
      end
      # ELSE

      # Record the total number of animals moved off the farm (:cMove_off)
      movements[off_index, 7] = total_moved

      # Randomly generate the states of those animals
      # n = (cS_init, cE_init, cI_init), N = cMove_off
      move_states = rng_mvhyper(movements[off_index, 3:5], total_moved)

      # Record the states moved off the farm
      # :cS_move_off, :cE_move_off, :cI_move_off
      movements[off_index, 8:10] .= move_states


      # FOR each farm that had animals moved onto it from this off-farm
       # IF it is in cheshire
        # IF the off farm moved animals to more than one farm
          # Extract the movements for this on-farm
          # Calculate the total movements
          # Generate the states
        # ELSE just calculate states
        # Update the states of the farm that recieved animals
        # Reduce the moved states by those that just got assigned
      # Look at the next on-farm

      # Extract all the row ids of all the farms that recieved animals from this farm
      on_indices = off_moves[:, :on_row_id]


      for idx in on_indices
        # IF the farm is in Cheshire
        if (isequal(idx, missing) == false)

          # IF there is more than one farm that recieved animals from this farm
          if (size(off_moves, 1) > 1)
            # Extract the movements for this on-farm
            moves_on = @where(off_moves, isequal.(:on_row_id, idx))
            # Calculate the total movements
            total_moved = sum(moves_on[:, :n_moves_off])
            # Generate the states
            states_on = rng_mvhyper(move_states, total_moved)
          else # ELSE there is only one farm
            states_on = move_states
          end

          # Update the states of the farm that recieved animals
          # :cS_move_on, :cE_move_on, :cI_move_on
          movements[idx, 11:13] = movements[idx, 11:13] + states_on


          # Record Movement details
          record_of_movements[rom_tr, :] .= [t, off_index, idx,
                                             move_states[1], move_states[2], move_states[3],
                                             states_on[1], states_on[2], states_on[3],
                                             total_moved]
          rom_tr += 1

          # Update the running total of the moved states
          # Reduce the moved states by those that just got assigned
          move_states = move_states - states_on

        end #end of if in cheshire

      end #end of for farms that recieved animals

      if sum(move_states) > 0
        record_of_movements[rom_tr, :] .= [t, off_index, -1.,
                                         move_states[1], move_states[2], move_states[3],
                                         move_states[1], move_states[2], move_states[3],
                                         sum(move_states)]
        rom_tr += 1
      end

    end #end of for off_cphs

    # Calculate how many moves could not be processed
    not_worked = size(new_off_cphs, 1)

    # If it is the same number as the last loop
    # (as in no new moves got processed)
    # Record the moves that didn't work and end the loop
    if (not_worked == size(off_cphs, 1))
      println("Stuck in a loop. ")
      println("t = ", t, " and ", not_worked, " farms did not work. ")

      for Fidx in new_off_cphs[:, :off_row_id]
        F_moves = @where(moves_data, isequal.(:off_row_id, Fidx))
        append!(unused_movements_F, F_moves)
      end

      break
    end

    # IF there were no moves left that didn't work
    # OR some more farms were able to be processed this loop

    # Update the list of farms that still need to processed
    off_cphs = deepcopy(new_off_cphs)

    new_off_cph_idxs = new_off_cphs[:, :off_row_id]

    # These farms did not work because they relied on other movements
    # in the same week to be processed first.
    # Now we update the "initial" states to reflect that these movements
    # have occured, but only of the farms that still need to be processed.
    movements[new_off_cph_idxs, 3:5] .=
      (prev_states[new_off_cph_idxs, 18:20] -
                                          movements[new_off_cph_idxs, 8:10] +
                                          movements[new_off_cph_idxs, 11:13])



    movements[new_off_cph_idxs, 6] = sum(movements[new_off_cph_idxs, 3:5], dims = 2)

  end # end of while

  record_of_movements = record_of_movements[findall(record_of_movements[:,1] .> 0), :]

  return ([movements, unused_movements_MR, unused_movements_F], record_of_movements)
end

# movements = generate_moves(t = 1, prev_states = temp_initial_states,
#                             hist_moves = df_chesh_movements_final,
#                             unused_movements_MR = unused_moves_MR,
#                             unused_movements_F = unused_moves_F)
# unused_moves_MR_t = movements[2]
# unused_moves_F_t = movements[3]
# movements = movements[1]


#########################################
### Function to generate timestep t+1 ###
#########################################

function update_states(;t, row_id, epi_pars::Array{Float64, 1},
                        prev_state, prev_track, prev_pers,
                        movements_t, births_and_deaths,
                        unused_death_t)

  ############################
  ### Establish Parameters ###
  ############################

  β_c = epi_pars[1]
  β_b = epi_pars[2]
  γ = epi_pars[3]

  F = epi_pars[4]
  ϵ = epi_pars[5]

  ρ = epi_pars[6]
  ρ_E = epi_pars[7]

  θ_bb = epi_pars[8]
  θ_bd = epi_pars[9]


  #####################
  ### Find the farm ###
  #####################

  holdings = 1:total_holdings

  holding_i_prev_state = deepcopy(prev_state[row_id, :])
  holding_i_prev_track = deepcopy(prev_track[row_id, :])
  holding_i_prev_pers = deepcopy(prev_pers[row_id, :])

  movements_i = movements_t[row_id, :]

  c_birth_and_death_df = @where(births_and_deaths, :row_id .== row_id, :week_no .== t)

  ##############################
  ### Initialise the records ###
  ##############################

  new_state = deepcopy(holding_i_prev_state)
  new_state[Not([1,2,44,45,46,47,48,49])] .= 0

  new_track = deepcopy(holding_i_prev_track)
  new_track[Not([1,2,30,31,32,33,34,35])] .= 0

  new_pers = deepcopy(holding_i_prev_pers)
  new_pers[Not([1,2,9,10,11,12,13,14])] .= 0

  ##############################
  ### Extract initial states ###
  ##############################

  cS_init = holding_i_prev_state[18] #cS_final from prev timestep
  cE_init = holding_i_prev_state[19]
  cI_init = holding_i_prev_state[20]
  cN_init = cS_init + cE_init + cI_init

  bS_init = holding_i_prev_state[27] #bS_final from prev timestep
  bE_init = holding_i_prev_state[28]
  bI_init = holding_i_prev_state[29]
  bN_init = bS_init + bE_init + bI_init

  pcS_init = holding_i_prev_state[35] #pcS_final from prev timestep
  pcE_init = holding_i_prev_state[36]
  pcI_init = holding_i_prev_state[37]

  pbS_init = holding_i_prev_state[41] #pbS_final from prev timestep
  pbE_init = holding_i_prev_state[42]
  pbI_init = holding_i_prev_state[43]

  next_test_week_prev = holding_i_prev_track[17]
  move_res_init = holding_i_prev_state[31] #move_res_final from prev timestep

  p_env_prev = holding_i_prev_pers[4] #p_env_cur from prev timestep

  ######################################################
  ### Calculate exposure and infection probabilities ###
  ######################################################

  if cN_init > 0
    c_exp_prob = 1 - exp( - cI_init/cN_init * β_c - F * p_env_prev)
  else
    c_exp_prob = 0
  end

  if bN_init > 0
    b_exp_prob = 1 - exp( - bI_init/bN_init * β_b - F * p_env_prev)
  else
    b_exp_prob = 0
  end

  c_inf_prob = 1 - exp( - γ)
  b_inf_prob = 1 - exp( - γ)

  ########################
  ### Update movements ###
  ########################

  cS_Moves = movements_i[3] # cS when the movements of that farm occured
  cE_Moves = movements_i[4]
  cI_Moves = movements_i[5]

  cS_move_off = movements_i[8]
  cE_move_off = movements_i[9]
  cI_move_off = movements_i[10]

  cS_move_on = movements_i[11]
  cE_move_on = movements_i[12]
  cI_move_on = movements_i[13]

  cS_postM = cS_init + cS_move_on - cS_move_off
  cE_postM = cE_init + cE_move_on - cE_move_off
  cI_postM = cI_init + cI_move_on - cI_move_off
  cN_postM = cS_postM + cE_postM + cI_postM

  cS_on_out = movements_i[14]
  cE_on_out = movements_i[15]
  cI_on_out = movements_i[16]


  #####################################################
  ### Draw new exposures and infectious individuals ###
  #####################################################

  c_new_exp = rand(Binomial(convert(Int64, cS_postM), c_exp_prob))
  b_new_exp = rand(Binomial(convert(Int64, bS_init), b_exp_prob))

  c_new_inf = rand(Binomial(convert(Int64, cE_postM), c_inf_prob))
  b_new_inf = rand(Binomial(convert(Int64, bE_init), b_inf_prob))

  cS_postEI = cS_postM - c_new_exp
  cE_postEI = cE_postM + c_new_exp - c_new_inf
  cI_postEI = cI_postM + c_new_inf
  cN_postEI = cS_postEI + cE_postEI + cI_postEI

  bS_postEI = bS_init - b_new_exp
  bE_postEI = bE_init + b_new_exp - b_new_inf
  bI_postEI = bI_init + b_new_inf
  bN_postEI = bS_postEI + bE_postEI + bI_postEI

  #######################
  ### Draw detections ###
  #######################

  E_detec = 0
  I_detec = 0

  next_test_week = deepcopy(next_test_week_prev)
  move_res_final = deepcopy(move_res_init)

  test_occur = 0

  if t == next_test_week_prev
    E_detec = rand(Binomial(convert(Int64, cE_postEI), (ρ * ρ_E)))
    I_detec = rand(Binomial(convert(Int64, cI_postEI), ρ))
    test_occur = 1

    if (I_detec + E_detec) > 0
      next_test_week = next_test_week_prev + 4
      move_res_final = 1 #true
    else
      next_test_week = next_test_week_prev + 24
      move_res_final = 0 #false
    end
  end

  cS_postDet = cS_postEI
  cE_postDet = cE_postEI - E_detec
  cI_postDet = cI_postEI - I_detec
  cN_postDet = cS_postDet + cE_postDet + cI_postDet

  ##############################
  ### Draw births and deaths ###
  ##############################

  # c_birth = 0
  if size(c_birth_and_death_df, 1) == 0
    c_birth = 0
  else
    c_birth = c_birth_and_death_df[1, :n_births]
  end

  # c_death = 0
  if size(c_birth_and_death_df, 1) == 0
    c_death = 0
  else
    c_death = c_birth_and_death_df[1, :n_deaths]
  end

  #println(county_idx, " ", parish_idx, " ", holding_idx)
  #println("c_S_new ", c_S_new, " c_E_new ", c_E_new, " c_I_new ", c_I_new, " c_death ", c_death)

  if (cN_postDet >= c_death)
    c_S_death, c_E_death, c_I_death = rng_mvhyper([cS_postDet, cE_postDet, cI_postDet], c_death)
  else
    c_S_death = 0
    c_E_death = 0
    c_I_death = 0
    # println("county_idx = ", county_idx, " parish_idx = ", parish_idx, " holding_idx = ", holding_idx, ", ")
    println("Row_id = ", row_id, " lost ", c_death, " death(s). ")
    append!(unused_death_t, c_birth_and_death_df)
  end


  b_birth = rand(Poisson(bN_postEI * θ_bb))

  b_S_death = rand(Binomial(convert(Int64, bS_postEI), θ_bd))
  b_E_death = rand(Binomial(convert(Int64, bE_postEI), θ_bd))
  b_I_death = rand(Binomial(convert(Int64, bI_postEI), θ_bd))

  cS_final = cS_postDet - c_S_death + c_birth
  cE_final = cE_postDet - c_E_death
  cI_final = cI_postDet - c_I_death

  bS_final = bS_postEI - b_S_death + b_birth
  bE_final = bE_postEI - b_E_death
  bI_final = bI_postEI - b_I_death

  #########################
  ### Record new states ###
  #########################

  new_state[1] = t

  new_state[Not([1,2,44,45,46,47,48,49])] .= [cS_init, cE_init, cI_init,
                                              cS_Moves, cE_Moves, cI_Moves,
                                              cS_postM, cE_postM, cI_postM,
                                              cS_postEI, cE_postEI, cI_postEI,
                                              cS_postDet, cE_postDet, cI_postDet,
                                              cS_final, cE_final, cI_final,
                                              bS_init, bE_init, bI_init,
                                              bS_postEI, bE_postEI, bI_postEI,
                                              bS_final, bE_final, bI_final,
                                              move_res_init, move_res_final,
                                              pcS_init, pcE_init, pcI_init,
                                              -100, -100, -100,
                                              pbS_init, pbE_init, pbI_init,
                                              -100, -100, -100]

  ###################################
  ### Record tracking information ###
  ###################################

  new_track[1] = t

  new_track[Not([1,2,30,31,32,33,34,35])] .= [cS_move_on, cE_move_on, cI_move_on,
                                              cS_move_off, cE_move_off, cI_move_off,
                                              cS_on_out, cE_on_out, cI_on_out,
                                              c_new_exp, c_new_inf, b_new_exp, b_new_inf,
                                              test_occur, next_test_week,
                                              E_detec, I_detec,
                                              c_birth, c_S_death, c_E_death, c_I_death,
                                              b_birth, b_S_death, b_E_death, b_I_death,
                                              -1000, -1000]

  ###################################
  ### Record tracking information ###
  ###################################

  new_pers[1] = t

  new_pers[Not([1,2,9,10,11,12,13,14])] .= [p_env_prev, -1000,
                                             c_exp_prob, b_exp_prob,
                                             c_inf_prob, b_inf_prob]

  return([new_state, new_track, new_pers, unused_death_t])
end

# update = update_states(t = 1, row_id = 1, epi_pars = epi_pars,
#                         prev_state = temp_initial_states,
#                         prev_track = temp_initial_track,
#                         prev_pers = temp_initial_pers,
#                         movements_t = movements,
#                         births_and_deaths = df_chesh_birth_and_death,
#                         unused_death_t = unused_deaths)
