#Gurobi problem 1


using Distributed
rmprocs(workers())
#addprocs(16)
@show workers()
@everywhere using JuMP, Gurobi
using XLSX,CSV,DataFrames
using Random
using Dates





#read in data
##########JuMP specific code ###############

# Create a new model
model=Model(Gurobi.Optimizer)
#first stage 
@variable(model, y[[0],m=M,k=K[1]], Bin) #binary for small gasifiers and turbines
@constraint(model, y .<= 1)
@constraint(model, -y.<= 0)
#symmetry breaking constraints
@constraint(model, [m=M,k=2:length(K[1])], y[0,m,k] <= y[0,m,k-1])
#second stage
#z is the auxiliary investment variable
@variable(model, z[J,M,S], Int)
@constraint(model,[j=[0],m=M,s=S], z[j,m,s]<=sum(y[j,m,k] for k in K[1]))
@constraint(model,[j=[0],m=M,s=S], -z[j,m,s]<=-(sum(y[j,m,k] for k in K[1])))
@constraint(model,[j=1:num_processing_sites,m=M,s=S], z[j,m,s]<=0)
@constraint(model,[j=1:num_processing_sites,m=M,s=S], -z[j,m,s]<=0)
@variable(model, x[I,J,T,S])
@constraint(model, x .<= 1)
@constraint(model, -x .<= 0)
@variable(model, q[I,T,S]) #fraction of unserved demand for site i
@constraint(model, q .<= 1)
@constraint(model, -q .<= 0)
@constraint(model, [i=I,t=T,s=S], sum(x[i,j,t,s] for j=J) +q[i,t,s] <= 1)
@constraint(model, [i=I,t=T,s=S], -sum(x[i,j,t,s] for j=J) -q[i,t,s] <= -1)
@variable(model, v[J,M,T,S], Int)
@constraint(model, -v[J,M,T,S] .<= 0)
@variable(model, w[J,J,M,T,S], Int)
@constraint(model, -w[J,J,M,T,S] .<= 0)
@constraint(model, [j=J,m=M,t=T,s=S], v[j,m,t,s] <= z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t))
@constraint(model, [j=J,m=M,t=T,s=S], -v[j,m,t,s] <= -(z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t)))
@constraint(model, [j=J,t=T,s=S], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= sum(v[j,m,t,s]*u[m] for m in M)) #u[m] is the production rate of module type m
@constraint(model, [j=J,t=T,s=S], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= umax)


@objective(model, Max, -1*(sum(p[s]*(sum(sum(c[i,j,t]*d[i,t,s]*x[i,j,t,s] for i in I) + sum(h[j,jprime,m,t]*w[j,jprime,m,t,s] for jprime in J, m in M) for j in J, t in T) + sum(penalty_cost[i,t]*q[i,t,s]*d[i,t,s] for i in I, t in T)) for s in S)
+sum(g[m]*y[0,m,k] for m in M for k in K[1])))

#set the MIPGap to 0.1%
set_optimizer_attribute(model, "MIPGap", 0.005)
set_optimizer_attribute(model, "TimeLimit", 3600*3)
optimize!(model)

#save results
currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
dirstring="$(num_demand_sites)_$(num_processing_sites)_$(num_scenarios)_$(num_time_periods)_"*currenttime
using Base.Filesystem
#= # Get the full path of the current file
full_path = @__FILE__
# Get just the file name
file_name = basename(full_path)
open(file_name[1:end-3]*"_solution.txt", "w") do io
    for x in all_variables(model)
        println(io, value(x))
    end
end =#

XLSX.openxlsx(dirstring*"results.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "new_sheet")
    sheet["A1"] = "num_demand_sites"
    sheet["A2"] = num_demand_sites
    sheet["B1"] = "num_processing_sites"
    sheet["B2"] = num_processing_sites
    sheet["C1"] = "num_scenarios"
    sheet["C2"] = num_scenarios
    sheet["D1"] = "num_time_periods"
    sheet["D2"] = num_time_periods
    sheet["E1"] = "objective_value"
    sheet["E2"] = objective_value(model)
    sheet["F1"] = "objective_bound"
    sheet["F2"] = objective_bound(model)
    sheet["G1"] = "MIPGap"
    sheet["G2"] = abs(objective_bound(model)-objective_value(model))/abs(objective_value(model))
    sheet["H1"] = "termination_status"
    sheet["H2"] = string(termination_status(model))
    sheet["J1"] = "time"
    sheet["J2"] = solve_time(model)
    sheet["K1"] = "algorithm"
    sheet["K2"] = "Gurobi"
    
 
end
