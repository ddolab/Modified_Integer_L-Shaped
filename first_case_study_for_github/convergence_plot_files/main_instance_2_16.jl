#create the directory at this level
rng_seed=4
#= #save results (this doesn't work because this part of the code is run every problem instance.)
currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
dirstring="rng$(rng_seed)_"*currenttime
mkdir(dirstring)
current_dir=pwd()
cd(current_dir*"/"*dirstring)
using Base.Filesystem =#

#convert from input parameter to algorithm, demand sites, processing sites, scenarios, time periods
input_parameter=16
if input_parameter <= 9
    algorithm=1
elseif input_parameter <= 18
    algorithm=2
elseif input_parameter <= 27
    algorithm=3
end
#= num_demand_sites = 50
num_processing_sites = 10
num_scenarios = 125
num_time_periods = 24 =#
modified_input_parameter=input_parameter % 9 #now lies from 0 to 8
if modified_input_parameter <= 2
    num_processing_sites=5
    num_demand_sites=25
elseif modified_input_parameter <= 5
    num_processing_sites=10
    num_demand_sites=50
else
    num_processing_sites=15
    num_demand_sites=75
end
modified_input_parameter=modified_input_parameter % 3 #now lies from 0 to 2
if modified_input_parameter <= 0
    num_scenarios=32
elseif modified_input_parameter <= 1
    num_scenarios=64
else
    num_scenarios=96
end
num_time_periods=24
include("problem_parameters_new_costs.jl")
if algorithm == 1
    include("Gurobi.jl")
elseif algorithm == 2
    include("L_shaped_SolveMPOnce.jl")
else
    include("modified_L_shaped_SolveMPOnce.jl")
end
