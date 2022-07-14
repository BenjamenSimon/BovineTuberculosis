
include("Likelihood.jl")
include("DataUpdaters.jl")

#########################################
### Constructing the likelihood array ###
#########################################

llh_array = zeros(size(combi_array[1], 1), 360, 13)

p_env_llh_array = zeros(size(combi_array[4], 1), 360, 2)

scope = [1, 360, 1:size(combi_array[1], 1), 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array, record_of_movements,
                                                          epi_params_true, dict_of_movements, f_to_p_dict)

############################
### Testing the Updaters ###
############################

### propose_Move_SE

combi_array_prime, log_q_ratio, scope, Move_SE_track = propose_Move_SE(combi_array, f_to_p_dict)


### propose_Move_EI

combi_array_prime, log_q_ratio, scope, Move_EI_track = propose_Move_EI(combi_array, epi_params_true, f_to_p_dict)


### propose_AddRem_SE

combi_array_prime, log_q_ratio, scope, AddRem_SE_track = propose_AddRem_SE(combi_array, f_to_p_dict)


### propose_AddRem_EI

combi_array_prime, log_q_ratio, scope, AddRem_EI_track = propose_AddRem_EI(combi_array, epi_params_true, f_to_p_dict)


### propose_AddRem_Det

combi_array_prime, log_q_ratio, scope, AddRem_Det_track = propose_AddRem_Det(combi_array, epi_params_true, f_to_p_dict)


### propose_AddRem_Deaths

combi_array_prime, log_q_ratio, scope, AddRem_Deaths_track = propose_AddRem_Deaths(combi_array, epi_params_true, f_to_p_dict)

# Investigate as produces no valid update


### generate_new_movement and propose_AddRem_Movements

generate_new_movement(combi_array, 102, 110, epi_params_true, record_of_movements, dict_of_movements, f_to_p_dict, ids_to_pos_dict)


combi_array_prime, log_q_ratio, scope, AddRem_Movements_track = propose_AddRem_Movements(combi_array, epi_params_true, record_of_movements, dict_of_movements, f_to_p_dict, ids_to_pos_dict)


### update_data_AddRem_penv

combi_array_prime, log_q_ratio, scope, AddRem_penv_track = propose_AddRem_penv(combi_array, epi_params_true, f_to_p_dict)





#################################
### Benchmarking the Updaters ###
#################################

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000

@benchmark begin

end



using ProfileView

@profview begin

end
