#create the directory at this level

#= #save results (this doesn't work because this part of the code is run every problem instance.)
currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
dirstring="rng$(rng_seed)_"*currenttime
mkdir(dirstring)
current_dir=pwd()
cd(current_dir*"/"*dirstring)
using Base.Filesystem =#

#convert from input parameter to algorithm, demand sites, processing sites, scenarios, time periods
input_parameter=parse(Int64,ARGS[1])
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
if modified_input_parameter <= 2 #24 periods
    days_per_scheduling_horizon = 1
    dt = 1  # length of time period (h)
elseif modified_input_parameter <= 5 #36 periods
    days_per_scheduling_horizon = 3
    dt = 2  # length of time period (h)
else #48 periods
    days_per_scheduling_horizon = 2
    dt = 1  # length of time period (h)
end
modified_input_parameter=modified_input_parameter % 3 #now lies from 0 to 2
if modified_input_parameter <= 0
    num_scenarios=8
elseif modified_input_parameter <= 1
    num_scenarios=16
else
    num_scenarios=24
end

rng_seed=5
include("problem_input.jl")
if algorithm == 1
    include("Gurobi.jl")
elseif algorithm == 2
    include("L_shaped.jl")
else
    include("modified_L_shaped.jl")
end
