#Gurobi problem 1


using Distributed
rmprocs(workers())
#addprocs(16)
@show workers()
@everywhere using JuMP, Gurobi






#read in data
##########JuMP specific code ###############
# Create model.
modelMP = direct_model(Gurobi.Optimizer())


## first stage

@variables modelMP begin
  V[i=1:I] >= 0  # capital cost
  lambda[i=1:I, l=1:1+L_i[i]] >= 0  # coefficient for piecewise-linear approximation
  w[i=1:I, l=1:1+L_i[i]], Bin  # 1 if capacity lies in inteval
  x[i=1:I], Bin  # 1 if process built
  x_[j=1:J], Bin  # 1 if resource inventory built
  #new variables 
  u[i=1:I, k=1:K], Bin
  u_[j=1:J, k=1:K], Bin
  Cmaster[i=1:I] >= 0
  C_master[j=1:J] >= 0
end
X = vcat(x...,x_...,u...,u_...) #vector of all first stage state variables
#Extra Network design constraints:
#Extra Network design constraints:
@constraints modelMP begin
    conHD[i=[9,10]], x[i] <= sum(x[i] for i = [15,23])
    conMS[i=[23]], x[i] <= sum(x[i] for i = [8])
    conP[i=[5,10,14]], x[i] <= sum(x[i] for i = [7,13,18,23,2,4,15])
    conSG[i=[2,3,16,19,20,22]], x[i] <= sum(x[i] for i = [1])
    conH[i=[6,17]], x[i] <= sum(x[i] for i = [5,3,19])
    ConO[i=[12]], x[i] <= sum(x[i] for i = [11])
    conMethane[i=[18]], x[i] <= sum(x[i] for i = [17])
    conMethanol[i=[12]], x[i] <= sum(x[i] for i = [20,6])
    conTE[i=[12]], x[i] <= sum(x[i] for i = [21,22,20,19,16])
    #= prodHD[i=[23,15]], x[i] <= sum(x[i] for i = [9,10])
    prodMS[i=[8]], x[i] <= sum(x[i] for i = [23])
    prodSG[i=[1]], x[i] <= sum(x[i] for i = [2,3,16,19,20,22])
    prodH[i=[5,3,19]], x[i] <= sum(x[i] for i = [6,17]) #this is sorta wrong because we could just discharge H for free/no reward
    prodO[i=[11]], x[i] <= sum(x[i] for i = [12])
    prodMethane[i=[17]], x[i] <= sum(x[i] for i = [18])
    prodMethanol[i=[20,6]], x[i] <= sum(x[i] for i = [12])
    prodTE[i=[21,22,20,19,16]], x[i] <= sum(x[i] for i = [12]) #also sorta wrong for same reason =#
    conStorage, x_ .== 1
end
# Network design constraints:
@constraints modelMP begin
  ndconstr3[i=1:I], sum(u[i, k] for k = 1:K) == x[i]
  ndconstr4[j=1:J], sum(u_[j, k] for k = 1:K) == x_[j]
  ndconstr7[i=1:I], Cmaster[i] == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)
  ndconstr8[j=1:J], C_master[j] == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)
end

# Piecewise-linear cost functions:
@constraints modelMP begin
  lcconstr1[i=1:I], Cmaster[i] == sum(lambda[i, l] * (Cl[i, l-1] - Cl[i, l]) + Cl[i, l] * w[i, l] for l = 2:1+L_i[i])
  lcconstr2[i=1:I], V[i] == sum(lambda[i, l] * (Vl[i, l-1] - Vl[i, l]) + Vl[i, l] * w[i, l] for l = 2:1+L_i[i])
  lcconstr3[i=1:I, l=2:1+L_i[i]], lambda[i, l] <= w[i, l]
  lcconstr4[i=1:I], sum(w[i, l] for l = 2:1+L_i[i]) == x[i]
end

## second stage
# Variables:
@variables modelMP begin
  B[j=1:J, h=1:H, t=1:θ_max+T_h[h], s=1:S] >= 0  # resource input
  C[i=1:I, s=1:S] >= 0  # process capacity
  C_[j=1:J, s=1:S] >= 0  # inventory capacity
  P[i=1:I, h=1:H, t=1:θ_max+T_h[h], s=1:S] >= 0  # production or consumption of reference resource
  P_[i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h], s=1:S] >= 0  # production or consumption in operating mode
  Q[s=1:S, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # inventory
  Q_[s=1:S, j=1:J, h=1:H]  # inventory carried over to the next season
  Sale[s=1:S, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # sales or waste
  y[s=1:S, i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if mode selected
  z[s=1:S, i=1:I, m1=1:M, m2=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if switching from mode m1 to mode m2 at t
  Dslack[s=1:S, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # slack for unmet demand
end

# Network design constraints:
@constraints modelMP begin
  ndconstr1[s=1:S, i=1:I], C[i, s] <= Cmax[i] * x[i]
  ndconstr2[s=1:S, j=1:J], C_[j, s] <= C_max[j] * x_[j]
  ndconstr5[s=1:S, i=1:I], C[i, s] == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)
  ndconstr6[s=1:S, j=1:J], C_[j, s] == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)
end

# Mass balance constraitns:
@constraints modelMP begin
  mbconstr1[s=1:S, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] == (1 - epsilon[j, h]) * Q[s, j, h, t-1] + sum(rho[i, j, h] * P[i, h, t, s] for i = 1:I) + B[j, h, t, s] - Sale[s, j, h, t]
  mbconstr2[s=1:S, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] <= eta[i, h, t, s] * C[i, s]
  mbconstr3[s=1:S, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] <= C_[j, s]
  mbconstr4[s=1:S, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], B[j, h, t, s] <= Bmax[j, h, t]
  mbconstr5[s=1:S, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J_[j] == 1], Sale[s, j, h, t] + Dslack[s, j, h, t] >= D[j, h, t]
  #mbconstr5[j=1:J,h=1:H,t=θ_max+1:θ_max+T_h[h]; J_[j]==1], Sale[s, j,h,t] >= D[j,h,t]
  mbconstr6[s=1:S, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J0[j] == 1], Sale[s, j, h, t] == 0
end

# Mode-based operation:
@constraints modelMP begin
  moconstr1[s=1:S, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(y[s, i, m, h, t] for m = 1:M if M_i[i, m] == 1) == x[i]
  moconstr2[s=1:S, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] == sum(P_[i, m, h, t, s] for m = 1:M if M_i[i, m] == 1)
  moconstr3[s=1:S, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], CCmin[i, m] * y[s, i, m, h, t] <= P_[i, m, h, t, s]
  moconstr4[s=1:S, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], P_[i, m, h, t, s] <= CCmax[i, m] * y[s, i, m, h, t]
end

# Transition constraints:
@constraints modelMP begin
  trconstr1[s=1:S, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], -Delta[i, m] - bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1]) <= P_[i, m, h, t, s] - P_[i, m, h, t-1, s]
  trconstr2[s=1:S, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], P_[i, m, h, t, s] - P_[i, m, h, t-1, s] <= Delta[i, m] + bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1])
  trconstr3[s=1:S, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], sum(z[s, i, m1, m, h, t-1] for m1 = 1:M if TR_i[i, m1, m] == 1) - sum(z[s, i, m, m2, h, t-1] for m2 = 1:M if TR_i[i, m, m2] == 1) == y[s, i, m, h, t] - y[s, i, m, h, t-1]
  trconstr4[s=1:S, i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; TR_i[i, m1, m2] == 1], y[s, i, m2, h, t] >= sum(z[s, i, m1, m2, h, t-t_prime] for t_prime = 1:theta[i, m1, m2])
  trconstr5[s=1:S, i=1:I, m1=1:M, m2=1:M, m3=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; SQ_i[i, m1, m2, m3] == 1], z[s, i, m1, m2, h, t-theta_[i, m1, m2, m3]] == z[s, i, m2, m3, h, t]
end

# Continuity equations:
@constraints modelMP begin
  ctconstr1[s=1:S, i=1:I, m=1:M, h=1:H; M_i[i, m] == 1], y[s, i, m, h, θ_max] == y[s, i, m, h, θ_max+T_h[h]]
  ctconstr2[s=1:S, i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t] == z[s, i, m1, m2, h, t+T_h[h]]
  ctconstr3[s=1:S, i=1:I, m=1:M, h=1:H-1; M_i[i, m] == 1], y[s, i, m, h, θ_max+T_h[h]] == y[s, i, m, h+1, θ_max]
  ctconstr4[s=1:S, i=1:I, m1=1:M, m2=1:M, h=1:H-1, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t+T_h[h]] == z[s, i, m1, m2, h+1, t]
  ctconstr5[s=1:S, j=1:J, h=1:H], Q_[s, j, h] == Q[s, j, h, θ_max+T_h[h]] - Q[s, j, h, θ_max]
  ctconstr6[s=1:S, j=1:J, h=1:H-1], Q[s, j, h, θ_max] + n[h] * Q_[s, j, h] == Q[s, j, h+1, θ_max]
  ctconstr7[s=1:S, j=1:J], Q[s, j, H, θ_max] + n[H] * Q_[s, j, H] >= Q[s, j, 1, θ_max]
end



# Objective function:
@objective(modelMP, Max, -(sum(V[i] for i = 1:I) + sum(alpha[j] * x_[j] + beta[j] * C_master[j] for j = 1:J)
                         + sum(p_s* n[h] * (sum(delta[i, m] * y[s, i, m, h, t] + gamma[i, m] * P_[i, m, h, t, s] for i = 1:I, m = 1:M if M_i[i, m] == 1)
                                       + sum(phi[j] * B[j, h, t, s] for j = 1:J) + sum(psi[j] * Sale[s, j, h, t] for j = 1:J)) for h = 1:H, t = θ_max+1:θ_max+T_h[h], s = 1:S)
                         + sum(1E4 * Dslack[s, j, h, t] for j = 1:J, h = 1:H, t = θ_max+1:θ_max+T_h[h], s = 1:S))/1e3)

#set the MIPGap to 0.1%
set_optimizer_attribute(modelMP, "MIPGap", 0.005)
set_optimizer_attribute(modelMP, "TimeLimit", 3600*4)
set_optimizer_attribute(modelMP, "ScaleFlag", 2)
set_optimizer_attribute(modelMP, "NumericFocus", 3)
set_optimizer_attribute(modelMP, "Method", 2)
optimize!(modelMP)

#save results
currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
dirstring="$(S)_$(hours_per_scheduling_horizon/dt)_$(rng_seed)_$(algorithm)_"*currenttime
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
  sheet["A1"] = "num_time_periods"
  sheet["A2"] = hours_per_scheduling_horizon/dt
  sheet["B1"] = "rng_seed"
  sheet["B2"] = rng_seed
  sheet["C1"] = "num_scenarios"
  sheet["C2"] = num_scenarios
  sheet["D1"] = "hours_per_scheduling_horizon"
  sheet["D2"] = hours_per_scheduling_horizon
  sheet["E1"] = "objective_value"
  try
  sheet["E2"] = objective_value(modelMP)
  catch
  sheet["E2"] = "no_sol"
  end
  sheet["F1"] = "objective_bound"
  try
  sheet["F2"] = objective_bound(modelMP)
  catch
  sheet["F2"] = "no_sol"
  end
  sheet["G1"] = "MIPGap"
  try
      sheet["G2"] = abs(objective_bound(modelMP)-objective_value(modelMP))/abs(objective_value(modelMP))
      catch
      sheet["G2"] = "no_sol"
      end
  sheet["H1"] = "termination_status"
  sheet["H2"] = string(termination_status(modelMP))
  sheet["J1"] = "time"
  sheet["J2"] = solve_time(RMP)
  sheet["K1"] = "algorithm"
  sheet["K2"] = "L_shaped"
  

end
