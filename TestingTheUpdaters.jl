
include("Likelihood.jl")

#########################################
### Constructing the likelihood array ###
#########################################

llh_array = zeros(100, 360, 13)

p_env_llh_array = zeros(77, 360, 2)

scope = [1, 360, 1:100, 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array,
                                                            record_of_movements, epi_params_true, dict_of_movements)

############################
### Testing the Updaters ###
############################

# be careful to only input a valid proposed move
# ie. moving 1 when there are 0 is invalid

### update_data_Move_SE

combi_array_prime, scope_new, valid = update_data_Move_SE(combi_array, 1, 336, -1, 1, f_to_p_dict) #combi_array_cur, pos_id, t, Δ, num_SE_moved
combi_array_prime, scope_new, valid = update_data_Move_SE(combi_array, 1, 336, 1, 1, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_Move_SE(combi_array, 1, 336, -20, 1, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_Move_SE(combi_array, 1, 336, 20, 1, f_to_p_dict)


###update_data_Move_EI

combi_array_prime, scope_new, valid = update_data_Move_EI(combi_array, 1, 281, -1, 1, epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_Move_EI(combi_array, 1, 281, 1, 1, epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_Move_EI(combi_array, 1, 281, -20, 1, epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_Move_EI(combi_array, 1, 281, 20, 1, epi_params_true, f_to_p_dict)


### update_data_AddRem_SE

combi_array_prime, scope_new, valid = update_data_AddRem_SE(combi_array, 1, 336, 1, f_to_p_dict) #combi_array_cur, pos_id, t, Δ
combi_array_prime, scope_new, valid = update_data_AddRem_SE(combi_array, 1, 336, -1, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_SE(combi_array, 1, 336, 20, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_SE(combi_array, 1, 336, -20, f_to_p_dict)


### update_data_AddRem_EI

combi_array_prime, scope_new, valid = update_data_AddRem_EI(combi_array, 1, 281, -1, epi_params_true, f_to_p_dict) # combi_array_cur, position, t, Δ, epi_params, f_to_p_dict
combi_array_prime, scope_new, valid = update_data_AddRem_EI(combi_array, 1, 281, 1, epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_EI(combi_array, 1, 281, -20, epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_EI(combi_array, 1, 281, 20, epi_params_true, f_to_p_dict)


### update_data_AddRem_Det

combi_array_prime, scope_new, valid = update_data_AddRem_Det(combi_array, 1, 250, [2,0], epi_params_true, f_to_p_dict) #combi_array_cur, position, t, Δs, epi_params, f_to_p_dict
combi_array_prime, scope_new, valid = update_data_AddRem_Det(combi_array, 1, 250, [0,2], epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_Det(combi_array, 1, 250, [1,1], epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_Det(combi_array, 1, 250, [2,2], epi_params_true, f_to_p_dict)


### update_data_AddRem_Deaths

combi_array_prime, scope_new, valid = update_data_AddRem_Deaths(combi_array, 4, 71, [2,1,0], epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_Deaths(combi_array, 4, 71, [1,0,2], epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_Deaths(combi_array, 4, 71, [0,2,1], epi_params_true, f_to_p_dict)
combi_array_prime, scope_new, valid = update_data_AddRem_Deaths(combi_array, 4, 71, [1,1,1], epi_params_true, f_to_p_dict) # combi_array_cur, position, t, Δs, epi_params, f_to_p_dict


### generate_new_movement

generate_new_movement(combi_array, 4, 169, epi_params_true, record_of_movements, dict_of_movements, f_to_p_dict, ids_to_pos_dict) #combi_array_cur, position, t, movement_record, dict_of_movements, f_to_p_dict, ids_to_pos_dict


### update_data_AddRem_penv

update_data_AddRem_penv(combi_array, 1, 300, [31, 6], epi_params_true, f_to_p_dict, ids_to_pos_dict) # combi_array_cur, p_position, t, Δs, epi_params, f_to_p_dict, ids_to_pos_dict



#################################
### Benchmarking the Updaters ###
#################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

@benchmark begin
  update_data_Move_SE(combi_array, 1, 26, 50, 1)
end



using ProfileView

@profview begin

end
