
include("Likelihood.jl")

#########################################
### Constructing the likelihood array ###
#########################################

llh_array = zeros(100, 360, 13)

p_env_llh_array = zeros(77, 360, 2)

scope = [1, 360, 1:100, 1:13]

llh_array_cur, p_env_llh_array_cur = update_llh_array_ALL(scope, llh_array, p_env_llh_array, combi_array,
                                                            record_of_movements, epi_params_true)

############################
### Testing the Updaters ###
############################

combi_array_prime, scope_new, valid = update_data_Move_SE(combi_array, 1, 26, 50, 1)




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
