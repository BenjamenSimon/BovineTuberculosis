
include("Simulate_bTB.jl")

# xtraM is assumed to mean more movements have been processed by looping

# TODO:
# - Remove the globals


###################
### The Process ###
###################

## The Parameters ##

β_c = 0.002
β_b = 0.004
γ = 0.015

F = 0.004
ϵ = 0.05

ρ = 0.75
ρ_E = 0.2

θ_bb = 0.25/52
θ_bd = 0.25/52

epi_pars = [β_c, β_b, γ, F, ϵ, ρ, ρ_E, θ_bb, θ_bd]

## Initialise the country ##

Random.seed!(3)

Time = 360

country_res = create_country_res(Time = Time, initial_settings = df_chesh_initial_states)
country_track = create_country_track(Time = Time, initial_settings = df_chesh_initial_states)
country_pers = create_country_pers(Time = Time, initial_settings = df_chesh_initial_states)


## Missing events ##

function create_unused_X_templates()

  unused_moves_template = DataFrame(
    off_cph = Float64[],
    off_county = Float64[],
    off_parish = Float64[],
    off_holding = Float64[],
    on_cph = Float64[],
    on_county = Float64[],
    on_parish = Float64[],
    on_holding = Float64[],
    week_no = Float64[],
    n_moves_off = Int32[],
    off_county_idx = Union{Missing, Float64}[],
    off_parish_idx = Union{Missing, Float64}[],
    off_holding_idx = Union{Missing, Float64}[],
    off_row_id = Union{Missing, Int32}[],
    on_county_idx = Union{Missing, Float64}[],
    on_parish_idx = Union{Missing, Float64}[],
    on_holding_idx = Union{Missing, Float64}[],
    on_row_id = Union{Missing, Int32}[],
    off_move_res = Union{Missing, Float64}[],
    on_move_res = Union{Missing, Float64}[],
  )

  unused_deaths_template = DataFrame(
    cph = Int32[],
    week_no = Float64[],
    n_births = Float64[],
    n_deaths = Float64[],
    county = Float64[],
    parish = Float64[],
    holding = Float64[],
    county_idx = Float64[],
    parish_idx = Float64[],
    holding_idx = Float64[],
    row_id = Int32[],
  )

  return([unused_moves_template, unused_deaths_template])
end

unused_moves_template, unused_deaths_template = create_unused_X_templates()

## Functional objects ##

unused_moves_MR = deepcopy(unused_moves_template) # Movement Restrictions
unused_moves_F = deepcopy(unused_moves_template) # Failed Movements

unused_deaths = deepcopy(unused_deaths_template)

record_the_movements = DataFrame(
  t = Float64[],
  off_row_id = Float64[],
  on_row_id = Float64[],
  cS_gen = Float64[],
  cE_gen = Float64[],
  cI_gen = Float64[],
  cS_StateMoved = Float64[],
  cE_StateMoved = Float64[],
  cI_StateMoved = Float64[],
  total_moves_to_farm = Float64[]
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### NEW ONE ###

@time begin

  ## Extract previous timestep states ##
  prev_states_t = extract_initial_state(country_res = country_res,
                                       n_counties = n_counties, n_parishes = n_parishes,
                                       n_holdings = n_holdings, total_holdings = total_holdings)

  prev_track_t = extract_initial_track(country_track = country_track,
                                       n_counties = n_counties, n_parishes = n_parishes,
                                       n_holdings = n_holdings, total_holdings = total_holdings)

  prev_pers_t = extract_initial_pers(country_pers = country_pers,
                                       n_counties = n_counties, n_parishes = n_parishes,
                                       n_holdings = n_holdings, total_holdings = total_holdings)

  t = 1

  while t <= Time

    println("t = ", t)

    global prev_states_t
    global prev_track_t
    global prev_pers_t

    ## Create record array ##
    unused_moves_MR_t = deepcopy(unused_moves_template)
    unused_moves_F_t = deepcopy(unused_moves_template)

    ## Generate new movements ##
    gen_moves_t, temp_record_of_move = generate_moves(t = t, prev_states = prev_states_t,
                                          hist_moves = df_chesh_movements_final,
                                          unused_movements_MR = unused_moves_MR_t,
                                          unused_movements_F = unused_moves_F_t)
    movements_t = gen_moves_t[1]

    ## Record unused movements ##
    append!(unused_moves_MR, gen_moves_t[2])
    append!(unused_moves_F, gen_moves_t[3])
    ## Record successful movements ##
    temp_record_of_move = DataFrame(temp_record_of_move, [:t, :off_row_id, :on_row_id, :cS_gen, :cE_gen, :cI_gen,
                       :cS_StateMoved, :cE_StateMoved, :cI_StateMoved, :total_moves_to_farm])
    append!(record_the_movements, temp_record_of_move)

    ## Calculate updates to states ##

    for r in 1:total_holdings
      ## Create record array ##
      unused_deaths_t = deepcopy(unused_deaths_template)

      ## Calculate update for holding with row_id 'r'
      update_t = update_states(t = t, row_id = r, epi_pars = epi_pars,
                                prev_state = prev_states_t,
                                prev_track = prev_track_t,
                                prev_pers = prev_pers_t,
                                movements_t = movements_t,
                                births_and_deaths = df_chesh_birth_and_death,
                                unused_death_t = unused_deaths_t)

      prev_states_t[r, :] .= update_t[1]
      prev_track_t[r, :] .= update_t[2]
      prev_pers_t[r, :] .= update_t[3]

      append!(unused_deaths, update_t[4])
    end

    prev_states_t = add_parish_final_totals(current_state = prev_states_t,
                                             n_counties = n_counties, n_parishes = n_parishes)

    prev_track_t, prev_pers_t = add_parish_env(current_state = prev_states_t,
                                  current_track = prev_track_t,
                                  current_pers = prev_pers_t,
                                  n_counties = n_counties, n_parishes = n_parishes,
                                  epi_pars = epi_pars,
                                  scaling = area_of_parish)

    ## Record updates to states ##

    row_id_count = 1

    for c in 1:n_counties
      for p in 1:n_parishes[c]
        for h in 1:n_holdings[c][p]

          country_res[c][p][h][t, :] .= prev_states_t[row_id_count, :]
          country_track[c][p][h][t, :] .= prev_track_t[row_id_count, :]
          country_pers[c][p][h][t, :] .= prev_pers_t[row_id_count, :]

          row_id_count = row_id_count+1

        end
      end
    end

    global t = t+1

  end
end


###################
### The Results ###
###################

## Inspect the final results ##

end_states = extract_states(t = Time, country_res = country_res,
                                     n_counties = n_counties, n_parishes = n_parishes,
                                     n_holdings = n_holdings, total_holdings = total_holdings)

end_track = extract_track(t = Time, country_track = country_track,
                                    n_counties = n_counties, n_parishes = n_parishes,
                                    n_holdings = n_holdings, total_holdings = total_holdings)

end_pers = extract_pers(t = Time, country_pers = country_pers,
                                    n_counties = n_counties, n_parishes = n_parishes,
                                    n_holdings = n_holdings, total_holdings = total_holdings)

unused_moves_MR

unused_moves_F

unused_deaths

######################
#### Save results ####
######################

results_chesh = fill(0., (2172*Time), 49)
track_chesh = fill(0., (2172*Time), 35)
pers_chesh = fill(0., (2172*Time), 14)

for i in 1:Time
  states_t = extract_states(t = i, country_res = country_res,
                                       n_counties = n_counties, n_parishes = n_parishes,
                                       n_holdings = n_holdings, total_holdings = total_holdings)

  track_t = extract_track(t = i, country_track = country_track,
                                     n_counties = n_counties, n_parishes = n_parishes,
                                     n_holdings = n_holdings, total_holdings = total_holdings)

  pers_t = extract_pers(t = i, country_pers = country_pers,
                                   n_counties = n_counties, n_parishes = n_parishes,
                                   n_holdings = n_holdings, total_holdings = total_holdings)

  start = (2172*(i-1))+1
  fin = (2172*i)

  results_chesh[start:fin, :] = states_t
  track_chesh[start:fin, :] = track_t
  pers_chesh[start:fin, :] = pers_t
end



results_df = DataFrame(results_chesh, :auto)
# Add column names to dataframe
rename!(results_df,   [:t, :row_id,
                        :cS_init, :cE_init, :cI_init,
                        :cS_Moves, :cE_Moves, :cI_Moves,
                        :cS_postM, :cE_postM, :cI_postM,
                        :cS_postEI, :cE_postEI, :cI_postEI,
                        :cS_postDet, :cE_postDet, :cI_postDet,
                        :cS_final, :cE_final, :cI_final,
                        :bS_init, :bE_init, :bI_init,
                        :bS_postEI, :bE_postEI, :bI_postEI,
                        :bS_final, :bE_final, :bI_final,
                        :move_res_init, :move_res_final,
                        :pcS_init, :pcE_init, :pcI_init,
                        :pcS_final, :pcE_final, :pcI_final,
                        :pbS_init, :pbE_init, :pbI_init,
                        :pbS_final, :pbE_final, :pbI_final,
                        :county, :parish, :holding,
                        :county_idx, :parish_idx, :holding_idx])

track_df = DataFrame(track_chesh, :auto)
# Add column names to dataframe
rename!(track_df, [:t, :row_id,
                    :sus_on, :exp_on, :inf_on,
                    :sus_off, :exp_off, :inf_off,
                    :S_on_out, :E_on_out, :I_on_out,
                    :c_new_exp, :c_new_inf, :b_new_exp, :b_new_inf,
                    :test_occur, :next_test_week,
                    :E_detected, :I_detected,
                    :c_birth, :c_S_death, :c_E_death, :c_I_death,
                    :b_birth, :b_S_death, :b_E_death, :b_I_death,
                    :remaining_pressure, :new_pressure,
                    :county, :parish, :holding,
                    :county_idx, :parish_idx, :holding_idx])

pers_df = DataFrame(pers_chesh, :auto)
# Add column names to dataframe
rename!(pers_df, [:t, :row_id,
                      :p_env_prev,
                      :p_env_cur,
                      :c_exp_prob, :b_exp_prob,
                      :c_inf_prob, :b_inf_prob,
                      :county, :parish, :holding,
                      :county_idx, :parish_idx, :holding_idx])



record_the_movements

##### SAVE ######

CSV.write("Data/Set 1/results_df.csv", results_df, header=true)

CSV.write("Data/Set 1/track_df.csv", track_df, header=true)

CSV.write("Data/Set 1/pers_df.csv", pers_df, header=true)

CSV.write("Data/Set 1/record_of_movements.csv", record_the_movements, header=true)


CSV.write("Data/Set 1/unused_moves_MR.csv", unused_moves_MR, header=true)

CSV.write("Data/Set 1/unused_moves_F.csv", unused_moves_F, header=true)

CSV.write("Data/Set 1/unused_deaths.csv", unused_deaths, header=true)
