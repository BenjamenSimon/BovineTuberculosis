
include("../LoadPackages.jl")


######################
#### Load in Data ####
######################

cph_to_formatted_raw = load("Data/cph_to_formatted_chesh.rds")

df_chesh_initial_states_raw = load("Data/RealDataInit/df_initial_conditions.rds")

df_chesh_movements_raw = load("Data/RealDataInit/df_movements_aggregated_farm_pairs_week.rds")

df_chesh_tests_raw = load("Data/RealDataInit/df_tests_aggregated_cph_week_counts.rds")

df_chesh_births_raw = load("Data/RealDataInit/df_births_aggregated_farm_week.rds")

df_chesh_deaths_raw = load("Data/RealDataInit/df_deaths_aggregated_farm_week.rds")


#######################################
#### Small initial conditions edit ####
#######################################

rows_to_update_1 = findall(df_chesh_initial_states_raw.cph .== (62640016))
rows_to_update_2 = findall(df_chesh_initial_states_raw.cph .== (61510129))

df_chesh_initial_states_raw[rows_to_update_1, [:num_not_infected, :num_infected]] .= [570.0 4.0]

df_chesh_initial_states_raw[rows_to_update_2, [:num_not_infected, :num_infected]] .= [133.0 2.0]


#############################
#### Preprocess the Data ####
#############################

###~~~ IDs ~~~###

cph_to_formatted = cph_to_formatted_raw[:, [:cph, :county, :parish, :holding, :n_holding_in_county_parish, :county_idx, :parish_idx, :row_id]]

###~~~ Initial states ~~~###

df_chesh_initial_states = leftjoin(df_chesh_initial_states_raw[:, [1,2,3,4,7,8]], cph_to_formatted[:, [1,8]], on = [:cph => :cph],
                                                        makeunique=false, indicator=nothing, validate=(false, false))

###~~~ Tests ~~~###

df_chesh_tests_raw_subset = @subset(df_chesh_tests_raw, :week_no .< 361)

df_chesh_tests = leftjoin(df_chesh_tests_raw_subset[:, [1,3,4,5]], cph_to_formatted[:, [1,2,3,4,8]], on = [:cph => :cph],
                                                        makeunique=false, indicator=nothing, validate=(false, false))

###~~~ Births ~~~###

df_chesh_births_raw_subset = @subset(df_chesh_births_raw, :week_no .< 361)

df_chesh_births = leftjoin(df_chesh_births_raw_subset[:, [1,2,3,4,5,7]], cph_to_formatted[:, [1, 8]], on = [:on_cph => :cph],
                                                            makeunique=false, indicator=nothing, validate=(false, false))

###~~~ Deaths ~~~###

df_chesh_deaths_raw_subset = @subset(df_chesh_deaths_raw, :week_no .< 361)

df_chesh_deaths = leftjoin(df_chesh_deaths_raw_subset[:, [1,2,3,4,5,7]], cph_to_formatted[:, [1, 8]], on = [:off_cph => :cph],
                                                            makeunique=false, indicator=nothing, validate=(false, false))

###~~~ Movements ~~~###

df_chesh_movements_raw_subset = @subset(df_chesh_movements_raw, :week_no .< 361)

df_chesh_movements_raw_2 = leftjoin(df_chesh_movements_raw[:, [1,2,3,4,5,7,8,9,10,12]], cph_to_formatted[:, [1, 8]], on = [:first_off_cph => :cph],
                                                            makeunique=false, indicator=nothing, validate=(false, false))

rename!(df_chesh_movements_raw_2, :row_id => :first_off_row_id)

df_chesh_movements = leftjoin(df_chesh_movements_raw_2, cph_to_formatted[:, [1, 8]], on = [:last_on_cph => :cph],
                                                            makeunique=false, indicator=nothing, validate=(false, false))

rename!(df_chesh_movements, :row_id => :last_on_row_id)

df_chesh_movements.first_off_row_id = coalesce.(df_chesh_movements.first_off_row_id, -1.0)
df_chesh_movements.last_on_row_id = coalesce.(df_chesh_movements.last_on_row_id, -1.0)


########################
#### Results arrays ####
########################

###~~~ Results (States) ~~~###

temp_real_data_res = Int64.(zeros(2172, 360, 23))

real_data_results = NamedArray(temp_real_data_res)

setdimnames!(real_data_results, ("Farm", "Time", "Data"))

setnames!(real_data_results, [
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
  "county",
  "parish"
], 3)

for farm in 1:2172
    real_data_results[farm, :, [1,2,3]] .= hcat(collect(1:360), fill(Float64(farm), 360), fill(Float64(farm), 360))
end

###~~~ Track (Events) ~~~###

temp_real_data_track = Int64.(zeros(2172, 360, 24))

real_data_track = NamedArray(temp_real_data_track)

setdimnames!(real_data_track, ("Farm", "Time", "Data"))

setnames!(real_data_track,
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
    "number_tested",
    "number_slaughtered",
    "E_detected",
    "I_detected",
    "c_birth",
    "c_S_death",
    "c_E_death",
    "c_I_death",
    "number_movements_off",
    "number_deaths"
  ],
  3,
)

for farm in 1:2172
    real_data_track[farm, :, [1,2,3]] .= hcat(collect(1:360), fill(Float64(farm), 360), fill(Float64(farm), 360))
end


###~~~ Persistents (Probs) ~~~###

temp_real_data_pers = zeros(2172, 360, 5)

real_data_pers = NamedArray(temp_real_data_pers)

setdimnames!(real_data_pers, ("Farm", "Time", "Data"))

setnames!(real_data_pers,
  [
    "t",
    "row_id",
    "pos",
    "c_exp_prob",
    "c_inf_prob",
  ],
  3,
)

for farm in 1:2172
    real_data_track[farm, :, [1,2,3]] .= hcat(collect(1:360), fill(Float64(farm), 360), fill(Float64(farm), 360))
end


###~~~ Parish (States) ~~~###

first_h_of_p_ALL = Int64.(combine(first, groupby(cph_to_formatted, [:county, :parish]))[:, :row_id])

function calc_area_of_parish_all()

  num_of_parishes = size(first_h_of_p_ALL, 1)

  area_of_parish = fill(-1, num_of_parishes)

  for i in 1:num_of_parishes

    county_idx = @subset(cph_to_formatted, :row_id .== first_h_of_p_ALL[i])[1, 6]
    parish_idx = @subset(cph_to_formatted, :row_id .== first_h_of_p_ALL[i])[1, 7]

    num_of_farms = size(@subset(cph_to_formatted, :county_idx .== county_idx, :parish_idx .== parish_idx), 1)

    area_of_parish[i] = num_of_farms * 1000
  end

  return(area_of_parish)
end

area_of_parish_all = calc_area_of_parish_all()

temp_real_data_parish = zeros(size(area_of_parish_all, 1), 360, 17)

real_data_parish = NamedArray(temp_real_data_parish)

setdimnames!(real_data_parish, ("Parish", "Time", "Data"))

setnames!(real_data_parish, [
  "t",
  "parish_uid",
  "parish_pos",
  "pcS_init",
  "pcE_init",
  "pcI_init",
  "pcS_final",
  "pcE_final",
  "pcI_final",
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
    real_data_parish[parish, :, 1] = 1:360
    real_data_parish[parish, :, 2] .= parish
    real_data_parish[parish, :, 3] .= parish
    real_data_parish[parish, :, 12] .= area_of_parish_all[parish]
    real_data_parish[parish, :, [15,16]] .= Array(@subset(cph_to_formatted, :row_id .== first_h_of_p_ALL[parish]))[:, [2,3]]
    real_data_parish[parish, :, 17] .= first_h_of_p_ALL[parish]
end

########################################
#### Fill in the initial conditions ####
########################################

for row_id in 1:size(real_data_results, 1)
  if row_id in unique(df_chesh_initial_states[:, :row_id])
    real_data_results[row_id, :, [4,6]] .= Array(@subset(df_chesh_initial_states, :row_id .== row_id))[:, [5,6]]
  end
  real_data_results[row_id, :, [22,23]] .= Array(@subset(cph_to_formatted, :row_id .== row_id))[:, [2,3]]
end


###################################
#### Fill in the events matrix ####
###################################

###~~~ Testing ~~~###

for row_id in 1:size(real_data_results, 1)
  times_when_tests_happened = Array(@subset(df_chesh_tests, :row_id .== row_id))[:,2]
  for t in 1:360
    if t in times_when_tests_happened
      real_data_track[row_id, t, [15,16]] = Array(@subset(df_chesh_tests, :row_id .== row_id, :week_no .== t))[:, [3,4]]
    end
  end
end

###~~~ Births ~~~###

for row_id in 1:size(real_data_results, 1)
  times_when_births_happened = Array(@subset(df_chesh_births, :row_id .== row_id))[:,1]
  for t in 1:360
    if t in times_when_births_happened
      real_data_track[row_id, t, [19]] .= Array(@subset(df_chesh_births, :row_id .== row_id, :week_no .== t))[:, [6]]
    end
  end
end

###~~~ Deaths ~~~###

for row_id in 1:size(real_data_results, 1)
  times_when_deaths_happened = Array(@subset(df_chesh_deaths, :row_id .== row_id))[:,1]
  for t in 1:360
    if t in times_when_deaths_happened
      real_data_track[row_id, t, [24]] .= Array(@subset(df_chesh_deaths, :row_id .== row_id, :week_no .== t))[:, [6]]
    end
  end
end


############################
#### Functional objects ####
############################

f_to_p_dict = Dict{Integer, Array}()

for row_id in 1:2172

  parish_num = Array(real_data_results[(real_data_results[:, 1, 2] .== row_id), 1, 23])

  parish_uid, parish_pos = Array(real_data_parish[(real_data_parish[:, 1, 16] .== parish_num), 1, 2:3])

  parish_members_uid = Array(real_data_results[(real_data_results[:, 1, 23] .== parish_num), 1, 2])

  parish_members_pos = Array(real_data_results[(real_data_results[:, 1, 23] .== parish_num), 1, 3])

  f_to_p_dict[row_id] = [Int64.(parish_uid), Int64.(parish_pos), Int64.(parish_members_uid), Int64.(parish_members_pos)]
                                                                                                       # farm_uid
end

f_to_p_dict

#####################################
#### Simulating events functions ####
#####################################

###~~~ Initialise week ~~~###

function initialise_week(t, real_data_results)

  for row_id in 1:2172
    real_data_results[row_id, t, 4:6] = copy(real_data_results[row_id, (t-1), 19:21])
  end

  return(real_data_results)
end


###~~~ Movements ~~~###

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

# dict_of_movements = Dict{Tuple{Int64, Int64}, Vector{Int64}}()
#
# for off_row_id in unique(df_chesh_movements.first_off_row_id)
#   for t in 1:360
#     dict_of_movements[(off_row_id, t)] = findall((df_chesh_movements[:,1] .== t) .& (df_chesh_movements[:,11] .== off_row_id))
#   end
#   println("off_row_id = ", off_row_id)
# end

# save("dict_of_movements_real.jld2", "dict", dict_of_movements)

dict_of_movements = load("dict_of_movements_real.jld2")["dict"]

function create_temp_movement_states(t, real_data_results)

  temp_movements_states = deepcopy(real_data_results[:, t, 4:6])

  return(temp_movements_states)
end

function process_movements_from_outside(t, real_data_results, real_data_track, dict_of_movements, df_chesh_movements, temp_movements_states)

  movements_t_outside = df_chesh_movements[ dict_of_movements[(-1, t)] , : ]

  for move in 1:size(movements_t_outside, 1)

        on_row_id = movements_t_outside[move, :last_on_row_id]

        S_E_I_on = rng_mvhyper([10800, 600, 600], movements_t_outside[move, :n_moves_off])

        temp_movements_states[on_row_id, 1:3] .+= S_E_I_on
        real_data_results[on_row_id, t, 10:12] .+= S_E_I_on

        real_data_track[on_row_id, t, 10:12] .+= S_E_I_on

  end

  return(real_data_results, real_data_track)
end

function process_movements_off(t, real_data_results, real_data_track, dict_of_movements, df_chesh_movements, failed_moves, temp_movements_states)

  recheck_ids = Int64[]

  ids_to_check = unique(df_chesh_movements.first_off_row_id)[unique(df_chesh_movements.first_off_row_id) .> 0]

  while size(ids_to_check, 1) > 0

    println("SIZE OF IDs TO CHECK = ", size(ids_to_check, 1))

    for off_row_id in ids_to_check
      movements_t = df_chesh_movements[ dict_of_movements[(off_row_id, t)] , : ]

      if size(movements_t, 1) > 0
        if size(movements_t, 1) == 1

          if ( sum(movements_t[:, :n_moves_off]) > sum(temp_movements_states[off_row_id, 1:3]) )
            push!(recheck_ids, copy(off_row_id))
            continue
          end

          on_row_id = movements_t[1, :last_on_row_id]

          S_E_I_on = rng_mvhyper(temp_movements_states[off_row_id, 1:3], movements_t[1, :n_moves_off])

          real_data_results[off_row_id, t, 7:9] .+= temp_movements_states[off_row_id, 1:3]
          real_data_results[off_row_id, t, 10:12] .+= (temp_movements_states[off_row_id, 1:3] - S_E_I_on)

          real_data_track[off_row_id, t, 7:9] .+= S_E_I_on

          if on_row_id > -1.0
            temp_movements_states[on_row_id, 1:3] .+= S_E_I_on
            real_data_results[on_row_id, t, 10:12] .+= S_E_I_on

            real_data_track[on_row_id, t, 4:6] .+= S_E_I_on
          end

        else # size(movements_t, 1) > 1

          if ( sum(movements_t[:, :n_moves_off]) > sum(temp_movements_states[off_row_id, 1:3]) )
            push!(recheck_ids, copy(off_row_id))
            continue
          end

          total_S_E_I_off = rng_mvhyper(temp_movements_states[off_row_id, 1:3], sum(movements_t[:, :n_moves_off]))

          real_data_results[off_row_id, t, 7:9] .+= temp_movements_states[off_row_id, 1:3]

          real_data_results[off_row_id, t, 10:12] .+= (temp_movements_states[off_row_id, 1:3] - total_S_E_I_off)

          real_data_track[off_row_id, t, 7:9] .+= total_S_E_I_off


          for move in 1:size(movements_t, 1)

            on_row_id = movements_t[move, :last_on_row_id]

            S_E_I_on = rng_mvhyper(total_S_E_I_off, movements_t[move, :n_moves_off])

            total_S_E_I_off -= S_E_I_on

            if on_row_id > -1.0
              temp_movements_states[on_row_id, 1:3] .+= S_E_I_on
              real_data_results[on_row_id, t, 10:12] .+= S_E_I_on

              real_data_track[on_row_id, t, 4:6] .+= S_E_I_on
            end
          end

        end
      end
    end

    if (size(recheck_ids, 1) == size(ids_to_check, 1))
      println("Stuck in a loop. ")
      println("t = ", t, " and ", size(recheck_ids, 1), " farms did not work. ")

      for off_row_id in ids_to_check
        push!(failed_moves, copy([off_row_id, t]))
      end

      break
    end


    ids_to_check = copy(recheck_ids)

    recheck_ids = Int64[]

  end #while

  return(real_data_results, real_data_track, failed_moves)
end


###~~~ Exposures and Infections ~~~###

function calc_exp_prop(States_init, p_env_prev, β, F)

  I = States_init[3]
  N = sum(States_init)

  if N > 0
    exp_prob = 1.0 - exp(- (I/N) * β - F * p_env_prev)
  else
    exp_prob = 0.0
  end

  return(exp_prob)
end

function calc_inf_prob(γ)

  inf_prob = 1.0 - exp(-γ)

  return(inf_prob)
end

function generate_infection_process(t, real_data_results, real_data_track, real_data_pers, real_data_parish, f_to_p_dict, β, F, γ)

  for row_id in 1:2172

    states = real_data_results[row_id, t, [10,11,12]]

    p_env_prev = real_data_parish[f_to_p_dict[(row_id)][1], t, 13]

    real_data_pers[row_id, t, [4]] .= calc_exp_prop(states, p_env_prev, β, F)

    real_data_pers[row_id, t, [5]] .= calc_inf_prob(γ)

    new_E = rand(Binomial(convert(Int64, states[1]), real_data_pers[row_id, t, 4]))

    new_I = rand(Binomial(convert(Int64, states[2]), real_data_pers[row_id, t, 5]))


    real_data_track[row_id, t, [13, 14]] .+= [new_E, new_I]

    real_data_results[row_id, t, [14, 15]] .+= [new_E, new_I]

    real_data_results[row_id, t, [13, 14]] .-= [new_E, new_I]

  end

  return(real_data_results, real_data_track, real_data_pers)
end


###~~~ Testing ~~~###

function process_testing(t, real_data_results, real_data_track)

  for row_id in 1:2172

    real_data_results[row_id, t, [16,17,18]] .= copy(real_data_results[row_id, t, [13,14,15]])

    num_tests, num_slaughtered = real_data_track[row_id, t, [15,16]]

    if num_tests > 0

      states = real_data_results[row_id, t, [13,14,15]]

      if sum(states) < num_tests
        println("Testing: Not enough animals to test. t = ", t, ", row_id = ", row_id)
      end

      if sum(states[2:3]) < num_slaughtered
        println("Testing: Not enough infected animals to slaughter. t = ", t, ", row_id = ", row_id)
        continue
      end

      positive_tests = rng_mvhyper(states[2:3], num_slaughtered)

      real_data_track[row_id, t, [17, 18]] .+= positive_tests

      real_data_results[row_id, t, [17,18]] .-= positive_tests

    end

  end

  return(real_data_results, real_data_track)
end


###~~~ Births & Deaths ~~~###

function process_births(t, real_data_results, real_data_track)

  for row_id in 1:2172

    real_data_results[row_id, t, [19,20,21]] .= copy(real_data_results[row_id, t, [16,17,18]])

    real_data_results[row_id, t, [19]] .+= real_data_track[row_id, t, [19]]

  end

  return(real_data_results, real_data_track)
end

function process_deaths(t, real_data_results, real_data_track)

  for row_id in 1:2172

    num_deaths = real_data_track[row_id, t, 24]

    states = real_data_results[row_id, t, [19,20,21]] #post births

    if sum(states) < num_deaths
      println("Deaths: Not enough animals to remove. t = ", t, ", row_id = ", row_id)
      continue
    end

    death_states = rng_mvhyper(states, num_deaths)

    real_data_track[row_id, t, [20,21,22]] .+= death_states

    real_data_results[row_id, t, [19,20,21]] .-= death_states

  end

  return(real_data_results, real_data_track)
end


###~~~ Parish ~~~###

function update_parish(t, real_data_results, real_data_track, real_data_pers, real_data_parish, f_to_p_dict, ϵ)

  for parish_row_id in 1:311

    parish_members_row_ids = f_to_p_dict[ real_data_parish[parish_row_id, 1, 17] ][3]

    real_data_parish[parish_row_id, t, 4:6] = sum(Array(real_data_results[parish_members_row_ids, t, 4:6]), dims = 1)

    real_data_parish[parish_row_id, t, 7:9] = sum(Array(real_data_results[parish_members_row_ids, t, 19:21]), dims = 1)

    # Calculate new parish environmental level

    new_pressure = rand(Poisson(real_data_parish[parish_row_id, t, 6]))
                                                  #pI

    remaining_pressure = rand(Binomial(round(Int, prod(real_data_parish[parish_row_id, t, 12:13])), (1-ϵ)))
                                                                # scaling * p_env_prev

    real_data_parish[parish_row_id, t, 14] = (remaining_pressure + new_pressure) / real_data_parish[parish_row_id, t, 12]
    # p_env_cur

    real_data_parish[parish_row_id, t, 10:11] = [remaining_pressure, new_pressure]

    if t < 360
      real_data_parish[parish_row_id, (t+1), 13] = copy(real_data_parish[parish_row_id, t, 14])
      # p_env_prev(t+1)
    end

  end

  return(real_data_results, real_data_track, real_data_pers, real_data_parish)
end



###################
#### Simulator ####
###################

real_data_results_TEST = deepcopy(real_data_results)
real_data_track_TEST = deepcopy(real_data_track)
real_data_pers_TEST = deepcopy(real_data_pers)
real_data_parish_TEST = deepcopy(real_data_parish)

failed_moves = []

for week in 1:360

  println("----------------------------------- WEEK = ", week, " -----------------------------------")

  println("----------------------------------- INITIALISING -----------------------------------")

  if week > 1
    real_data_results_TEST = initialise_week(week, real_data_results_TEST)
  end

  println("----------------------------------- CREATING TEMP MOVEMENTS -----------------------------------")

  temp_movement_states = create_temp_movement_states(week, real_data_results_TEST)

  println("----------------------------------- PROCESSING MOVEMENTS FROM OUTSIDE -----------------------------------")

  real_data_results_TEST, real_data_track_TEST = process_movements_from_outside(week, real_data_results_TEST, real_data_track_TEST, dict_of_movements, df_chesh_movements, temp_movement_states)

  println("----------------------------------- PROCESSING MOVEMENTS OFF -----------------------------------")

  real_data_results_TEST, real_data_track_TEST, failed_moves = process_movements_off(week, real_data_results_TEST, real_data_track_TEST, dict_of_movements, df_chesh_movements, failed_moves, temp_movement_states)

  println("----------------------------------- GENERATING INFECTION PROCESS -----------------------------------")

  real_data_results_TEST, real_data_track_TEST, real_data_pers_TEST = generate_infection_process(week, real_data_results_TEST, real_data_track_TEST, real_data_pers_TEST, real_data_parish_TEST, f_to_p_dict, 0.0004, 0.0004, 0.015)

  println("----------------------------------- PROCESSING TESTING -----------------------------------")

  real_data_results_TEST, real_data_track_TEST = process_testing(week, real_data_results_TEST, real_data_track_TEST)

  println("----------------------------------- PROCESSING BIRTHS -----------------------------------")

  real_data_results_TEST, real_data_track_TEST = process_births(week, real_data_results_TEST, real_data_track_TEST)

  println("----------------------------------- PROCESSING DEATHS -----------------------------------")

  real_data_results_TEST, real_data_track_TEST = process_deaths(week, real_data_results_TEST, real_data_track_TEST)

  println("----------------------------------- UPDATING PARISHES -----------------------------------")

  real_data_results_TEST, real_data_track_TEST, real_data_pers_TEST, real_data_parish_TEST = update_parish(week, real_data_results_TEST, real_data_track_TEST, real_data_pers_TEST, real_data_parish_TEST, f_to_p_dict, 0.05)

end
