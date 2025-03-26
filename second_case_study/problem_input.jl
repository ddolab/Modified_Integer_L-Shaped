
# -----------------------------------------------------------------------------
# Preprocess data.
# -----------------------------------------------------------------------------

# If not passed on from main.jl:
# using DataArrays, DataFrames, JLD
# case = "Case1"
# instance = "a"

# -----------------------------------------------------------------------------

# Processes: 1-22 = P1-P22 including P8a, 23 = P8b
# Resources: 1 = Biomass (kg), 2 = Waste (kg), 3 = Water (kg),
#            4 = Solar (P7: kWh power at 1 kWh/m^2/day with 1 m^2 PV panel area;
#                       P8: kWh power at 1 kWh/m^2/day with 1 m^2 mirror area;
#                       P11: kg oil at 1 kWh/m^2/day with one pond, 1000 m^2 area),
#            5 = Wind (kWh at 10 m/s wind velocity), 6 = Syngas (kg),
#            7 = CO2 (kg), 8 = Hydrogen (kg), 9 = Oxygen (kg),
#            10 = Gasoline (kg), 11 = Methanol (kg), 12 = Diesel (kg),
#            13 = Methane (kg), 14 = Oil (kg), 15 = Glycerol (kg),
#            16 = Molten_Salt (kWh), 17 = Heat_Duty (kWh),
#            18 = Elevated_Water (m^3 at 100 m height), 19 = Thermal (kWh),
#            20 = Power (kWh)
#            Note: [res] = unit for resource, [refres] = unit for reference resource
# Modes: 1 = off, 2 = startup, 3 = on, 4 = shutdown, >=5 = additional on modes
using DataFrames, JLD, CSV, Statistics
using XLSX,CSV,DataFrames
using Random
using Dates
case = "Case1"
instance = "b"
model = "MP"  # MP = multiscale model, AP = aggregate model
# Import data from files.
Cmax_dat = CSV.File(string("data/", case, "_input_Cmax.csv")) |> DataFrame # ([refres]/h)
C_max_dat = CSV.File(string("data/", case, "_input_C_max.csv")) |> DataFrame # ([res])
rho_dat = CSV.File(string("data/", case, "_input_rho.csv")) |> DataFrame # ([res]/[refres])
theta_dat = CSV.File(string("data/", case, "_input_theta.csv")) |> DataFrame # (h)
Cl_dat = CSV.File(string("data/", case, "_input_Cl.csv")) |> DataFrame # ([refres]/h)
Vl_dat = CSV.File(string("data/", case, "_input_Vl.csv")) |> DataFrame # (EUR)
D_gasoline = CSV.File(string("data/", case, "_input_D_gasoline.csv")) |> DataFrame # (kg/h)
D_diesel = CSV.File(string("data/", case, "_input_D_diesel.csv")) |> DataFrame # (kg/h)
D_power = CSV.File(string("data/", case, "_input_D_power.csv")) |> DataFrame # (kW)
solar = CSV.File(string("data/", case, "_input_solar.csv")) |> DataFrame # (kJ/m^2/h)
wind = CSV.File(string("data/", case, "_input_wind.csv")) |> DataFrame # (km/h)
water = CSV.File(string("data/", case, "_input_water.csv")) |> DataFrame # (kg/h)
alpha_dat = CSV.File(string("data/", case, "_input_alpha.csv")) |> DataFrame # (EUR)
beta_dat = CSV.File(string("data/", case, "_input_beta.csv")) |> DataFrame # (EUR/[refres])
delta_dat = CSV.File(string("data/", case, "_input_delta.csv")) |> DataFrame # (EUR/h)
gamma_dat = CSV.File(string("data/", case, "_input_gamma.csv")) |> DataFrame # (EUR/[refres])

Cmax_dat[13,2]=1e6
Cl_dat[42,3]=1e6
Cmax_dat[18,2]=1e6
Cl_dat[56,3]=1e6

input_data = "multiscale"  # "multiscale" or "detailed"
Random.seed!(rng_seed)
# Set parameters depending on selected level of input data granularity.
H = 1 # number of seasons
if input_data == "multiscale"
  days_per_season = 90
  hours_per_scheduling_horizon = 24 * days_per_scheduling_horizon

  dT = hours_per_scheduling_horizon  # maximum length of scheduling horizon (h)
  dT_h = [hours_per_scheduling_horizon for i in 1:H]  # scheduling horizon length for each season (h)
  n = [90 * 24 / hours_per_scheduling_horizon for i in 1:H] # number of repetitions of cyclic schedule
  @assert n[1] == floor(n[1])
  n = Int.(n)
end

# Replace NA values in data frames.
solar.solar[ismissing.(solar.solar)] .= 0

# Index sets:
I = maximum(Cmax_dat[:, 1])  # number of processes
J = maximum(C_max_dat[:, 1])  # number of resources
L = maximum(Cl_dat[:, 2]) - 1  # maximum number of line segments
M = 7 # maximum(theta_dat[:, 2])  # maximum number of operating modes
R = 1  # number of regions
K = 1 # number of sizes of each investment
S = num_scenarios # number of scenarios
p_s = 1/S  # probability of scenarios

# Sets for time representation:

if dt > minimum(theta_dat[:, 4]/6)
  println("Warning: Δt is larger than smallest transition time. Transition times smaller than Δt will be set to Δt.")
  for i = 1:size(theta_dat, 1)
    if theta_dat[i, 4] < dt
      theta_dat[i, 4] = dt
    end
  end
end
T = round(Int64, dT / dt)  # maximum number of time periods (does not include time periods <=0)
dK = maximum(theta_dat[:, 4]/6)  # maximum possible stay time (h)
θ_max = round(Int64, dK / dt)  # maximum number of stay time periods
T_h = round.(Int64, dT_h / dt)  # number of time periods in set T_h[h]
lifetime = 15  # lifetime of facilities (year)

# Set of resources with demand:
J_ = zeros(Int64, J)
J_[10] = 1  # gasoline
J_[12] = 1  # diesel
J_[20] = 1  # power

# Set of resources that must not be discharged:
J0 = zeros(Int64, J)
J0[1:6] .= 1
J0[17] = 1

# Sets of operating modes:
M_i = zeros(Int64, I, M)
for i = 1:size(theta_dat, 1)
  M_i[theta_dat[i, 1], theta_dat[i, 2]] = 1
end
for i in 1:I
  if M_i[i,3]==1
    M_i[i,5:7] .=1
  end
end
# Wind and solar power generation only in on mode.
M_i[4, 3] = 1
M_i[7, 3] = 1
M_i[8, 3] = 1

# Sets of possible transitions:
TR_i = zeros(Int64, I, M, M)
for i = 1:size(theta_dat, 1)
  TR_i[theta_dat[i, 1], theta_dat[i, 2], theta_dat[i, 3]] = 1
end
for i in 1:I
  if M_i[i,3]==1
    TR_i[i,3,5] =1
    TR_i[i,5,3] =1
    TR_i[i,3,6] =1
    TR_i[i,6,3] =1
    TR_i[i,3,7] =1
    TR_i[i,7,3] =1
    TR_i[i,5,6] =1
    TR_i[i,6,5] =1
    TR_i[i,5,7] =1
    TR_i[i,7,5] =1
    TR_i[i,6,7] =1
    TR_i[i,7,6] =1
  end
end

# Sets of fixed sequences:
SQ_i = zeros(Int64, I, M, M, M)
SQ_i[1:I, 1, 2, 3] .= 1  # startup
SQ_i[1:I, 3, 4, 1] .= 1  # shutdown
# Manually add shutdown sequences for processes with more than one on mode.

# Sets of line segments:
L_i = ones(Int64, I)
L_i[3] = 5
L_i[9] = 5
L_i[11] = 5
L_i[12] = 5
L_i[16] = 5
L_i[19] = 5

# Parameter Cmax[r,i], process capacities ([refres]):
# Assume 2% of area in Almeria (8774 km^2) can be used for PV panels, and 2% for algae production.
Cmax = zeros(Float64, I)
for i = 1:size(Cmax_dat, 1)
  Cmax[Cmax_dat[i, 1]] = dt * Cmax_dat[i, 2]
end
# Manually add region-dependent capacities if necessary, e.g. due to availability of land.

# Parameter C_max[r,j], resource inventory capacities ([res]):
C_max = zeros(Float64, J)
for i = 1:size(C_max_dat, 1)
  C_max[C_max_dat[i, 1]] = C_max_dat[i, 2]
end
Cnameplate = [[  1][k] * Cmax[i] for i = 1:I, k = 1:K]
C_nameplate = [[ 1][k] * C_max[j] for j = 1:J, k = 1:K]
# Manually add region-dependent capacities if necessary.

# Parameter epsilon[j,h], losses:
epsilon = zeros(Float64, J, H)
# Assume 2%/h loss for molten salt used for concentrated solar power.
epsilon[16, 1:H] .= dt * 0.02

# Parameter rho[i,j,h,t], conversion factor w.r.t. the reference resource ([res]/[refres]):
rho = zeros(Float64, I, J, H)
for i = 1:size(rho_dat, 1)
  rho[rho_dat[i, 1], rho_dat[i, 2], 1:H] .= rho_dat[i, 3]
end

# Parameter Bmax[r,j,h,t], resource availability, e.g. wind and solar ([res]):
# CO2 availability may also be set to be nonzero if we want to assume that CO2 discharged from other plants can be used.
Bmax = zeros(Float64, J, H, θ_max + T)
Bmax[1, 1:H, θ_max+1:θ_max+T] .= dt * 12000  # biomass  ***
Bmax[2, 1:H, θ_max+1:θ_max+T] .= dt * 150000  # waste
# Water:
start = 0
for h = 1:H
  global start
  for t = 1:T_h[h]
    Bmax[3, h, θ_max+t] = sum(sum(water[start+(seas-1)*dT_h[h]+(t-1)*dt+i, 2] for i = 1:dt) for seas = 1:n[h]) / n[h]
  end
  start = start + n[h] * dT_h[h]
end
Bmax[4, 1:H, θ_max+1:θ_max+T] .= 1E6*dt  # solar
Bmax[5, 1:H, θ_max+1:θ_max+T] .= 125000*dt  # wind
Bmax[18, :, :] = 0.01 * Bmax[3, :, :]  # elevated water (the part that can be readily used for hydropower, i.e. does not have to be pumped up first)
Bmax[7, 1:H, θ_max+1:θ_max+T] .= 1E6 *dt # CO2 ***CHECK***
# Aggregate Bmax to BBmax for (AP).
BBmax = zeros(Float64, J, H)
for h = 1:H
  BBmax[1:J, h] = sum(Bmax[1:J, h, t] for t = θ_max+1:θ_max+T_h[h])
end

# Parameter D[j,h,t], demand ([res]): ***CHECK***
D = zeros(Float64, J, 4, θ_max + T)
dratioGasoline = 0.1*0.5  # fraction of total gasoline demand to be met
dratioDiesel = 0.1*0.5  # fraction of total diesel demand to be met
dratioPower = 0.5*0.5  # fraction of total power demand to be met, 64 is feasible
# dratioGasoline = 0.2  # fraction of total gasoline demand to be met
# dratioDiesel = 0.1  # fraction of total diesel demand to be met
# dratioPower = 0.8  # fraction of total power demand to be met
# Gasoline:
for i = 1:size(D_gasoline, 1)
  D[10, D_gasoline[i, 1], θ_max+1:θ_max+T] .= dratioGasoline * dt * D_gasoline[i, 2] * (0.75 + 0.5 * rand())
end
# Diesel:
for i = 1:size(D_diesel, 1)
  D[12, D_diesel[i, 1], θ_max+1:θ_max+T] .= dratioDiesel * dt * D_diesel[i, 2] * (0.75 + 0.5 * rand())
end
# Power (average over each season's time periods):
start = 0
for h = 1:H
  global start
  for t = 1:T_h[h]
    D[20, h, θ_max+t] = sum(sum(dratioPower * D_power[start+(seas-1)*dT_h[h]+(t-1)*dt+i, 2] for i = 1:dt) for seas = 1:n[h]) / n[h] * (0.75 + 0.5 * rand())
  end
  start = start + n[h] * dT_h[h]
end
# Aggregate D to DD for (AP).
DD = zeros(Float64, J, H)
for h = 1:H
  DD[1:J, h] = sum(D[1:J, h, t] for t = θ_max+1:θ_max+T_h[h])
end

# Parameter eta[r,i,h,t], normalized capacity:
eta = ones(Float64, I, H, θ_max + T, S+1)
# Convert sun radiation (s in kWh/m^2/day) into eta. Note the data is given in kJ/m^2/h.
# For P7, eta=s. For P8, eta=s. For P11, eta=s.
# Average values for each season depending on the definition of the time horizons.
for s in 1:S
  for h = 1:H
    #pick a random day from the season currently considered to be the starting point of the sample we take.
    starting_day = (90 * (h - 1)) + rand(1:(days_per_season-days_per_scheduling_horizon+1))
    starting_hour = (starting_day - 1) * 24 + 1
    ending_hour = starting_hour + hours_per_scheduling_horizon - 1
    eta[7, h, θ_max+1:θ_max+T_h[h], s] = solar[starting_hour:dt:ending_hour, 2] * 24 / 3600
    eta[8, h, θ_max+1:θ_max+T_h[h], s] = solar[starting_hour:dt:ending_hour, 2] * 24 / 3600
    eta[11, h, θ_max+1:θ_max+T_h[h], s] = solar[starting_hour:dt:ending_hour, 2] * 24 / 3600
  end
  # Convert wind velocity (v in m/s) into eta.
  # For v<=5m/s, Power=0. For v>5, Power=(v-5)*250.
  # Use v=10 as reference. For v<=5m/s, eta=0. For v>5, eta=(v-5)/(10-5)=v/5-1.
  for h = 1:H
    #pick a random day from the season currently considered to be the starting point of the sample we take.
    starting_day = (90 * (h - 1)) + rand(1:(days_per_season-days_per_scheduling_horizon+1))
    starting_hour = (starting_day - 1) * 24 + 1
    ending_hour = starting_hour + hours_per_scheduling_horizon - 1
    # if sum(sum(wind[start+(seas-1)*dT_h[h]+(t-1)*dt+i,2]*1000/3600 for i=1:dt) for seas=1:n[h])/n[h] < 5  # *** used to be <5 ***
    #= if sum(sum(wind[start+(seas-1)*dT_h[h]+(t-1)*dt+i, 2] for i = 1:dt) for seas = 1:n[h]) / n[h] < 5
      eta[4, h, θ_max+t] = 0
    else =#
    eta[4, h, θ_max+1:θ_max+T_h[h], s] = (wind[starting_hour:dt:ending_hour, 2] .- 5) ./ (10 - 5)
    #eta[4, h, θ_max+t] = (sum(sum(wind[start+(seas-1)*dT_h[h]+(t-1)*dt+i, 2] for i = 1:dt) for seas = 1:n[h]) / n[h] - 5) / (8 - 5)
  end
end
eta[eta.<0] .= 0
eta_avg=[mean(eta[i,h,t,s] for s=1:S) for i=1:I,h=1:H,t=1:θ_max+T_h[1]]
#eta[:,:,:,S+1]=[mean(eta[i,h,t,s] for s=2:S) for i=1:I,h=1:H,t=1:θ_max+T_h[1]]

# Parameter theta[i,m1,m2], minimum stay time (in number of time periods):
theta = zeros(Int64, I, M, M)
for i = 1:size(theta_dat, 1)
  theta[theta_dat[i, 1], theta_dat[i, 2], theta_dat[i, 3]] = ceil.(Int64, theta_dat[i, 4] / dt/6)
end
theta[:, 2,3] .= ceil(4/dt)
theta[:, 3, 5:7] .= ceil(4/dt)
theta[:, 5:7, 3] .= ceil(4/dt)
theta[:, 5,6:7] .= ceil(4/dt)
theta[:, 6:7,5] .= ceil(4/dt)
theta[:, 6,7] .= ceil(4/dt)
theta[:, 7,6] .= ceil(4/dt)

# Parameter theta_[i,m1,m2,m3], fixed stay time (in number of time periods):
theta_ = zeros(Int64, I, M, M, M)
theta_[1:I, 1, 2, 3] = theta[1:I, 1, 2]
theta_[1:I, 3, 4, 1] = theta[1:I, 3, 4]
# Manullay add shutdown sequences for processes with more than one on mode.

# Parameter Cl[i,l], maximum capacity in interval ([refres]):
Cl = zeros(Float64, I, 1 + L)
for i = 1:size(Cl_dat, 1)
  Cl[Cl_dat[i, 1], Cl_dat[i, 2]] = dt * Cl_dat[i, 3]
end

# Parameter Vl[i,l], maximum capital cost in interval (EUR):
Vl = zeros(Float64, I, 1 + L)
for i = 1:size(Vl_dat, 1)
  Vl[Vl_dat[i, 1], Vl_dat[i, 2]] = Vl_dat[i, 3] / lifetime
end
for i in eachindex(Vl[:, 1])
  Vl[i, :] = Vl[i, :] * (0.75 + 0.5 * rand())
end

# Parameter CCmin[i,m], mode-dependent lower bound on production ([refres]):
# For wind and solar, it is zero. For chemical processes, assume 1% of maximum capacity of first line segment in on mode, zero otherwise.
CCmin = zeros(Float64, I, M)
CCmin[1:I, 3] = 0.03 * Cmax[1:I]


# Parameter CCmax[i,m], mode-dependent upper bound on production ([refres]):
# Use full capacity in on mode, zero otherwise.
CCmax = zeros(Float64, I, M)
CCmax[1:I, 3] = Cmax[1:I]*0.25
CCmin[1:I, 5] = CCmax[1:I, 3]
CCmax[1:I, 5] = Cmax[1:I]*0.5
CCmin[1:I, 6] = CCmax[1:I, 5]
CCmax[1:I, 6] = Cmax[1:I]*0.75
CCmin[1:I, 7] = CCmax[1:I, 6]
CCmax[1:I, 7] = Cmax[1:I]

CCmin[4, 3] = 0
CCmin[7, 3] = 0
CCmin[8, 3] = 0
CCmin[11, 3] = 0
CCmax[4, 3] = Cmax[4]
CCmax[7, 3] = Cmax[7]
CCmax[8, 3] = Cmax[8]
CCmax[11, 3] = Cmax[11]
# Parameter Delta[i,m], maximum rate change in the same mode ([res]):
# Full flexibility for wind and solar. For chemical processes, assume 30% of the entire range.
Delta = zeros(Float64, I, M)
Delta = 0.2 *dt * (CCmax - CCmin)
Delta[4, 3] = Cmax[4]
Delta[7, 3] = Cmax[7]
Delta[8, 3] = Cmax[8]
Delta[11, 3] = Cmax[11]

# Big-M parameter:
bigM = 1E6

# For approximating cost for storage tanks, we linearize the correlation 6839.8*V^0.65 ($, volume V in m^3).
# The linearized correlation is then 79000 + 113*V = 79000 + 113/rho*m (where rho here is the density and m is the mass).
# Parameter alpha[j], fixed capital cost for storage capacity (EUR):
alpha = zeros(Float64, J)
for i = 1:size(alpha_dat, 1)
  alpha[alpha_dat[i, 1]] = alpha_dat[i, 2] / lifetime * (0.75 + 0.5 * rand())
end
# Parameter beta[j], unit capital cost for storage capacity (EUR/[refres]):
beta = zeros(Float64, J)
for i = 1:size(beta_dat, 1)
  beta[beta_dat[i, 1]] = beta_dat[i, 2] / lifetime * (0.75 + 0.5 * rand())
end

# Parameter delta[i,m], fixed operating cost (EUR):
delta = zeros(Float64, I, M)
for i = 1:size(delta_dat, 1)
  delta[delta_dat[i, 1], delta_dat[i, 2]] = dt * delta_dat[i, 3] * (0.75 + 0.5 * rand())
end
Delta[:, 5:7] .= Delta[:, 3]
# Prameter gamma[i,m], unit operating cost (EUR/[refres])
gamma = zeros(Float64, I, M)
for i = 1:size(gamma_dat, 1)
  gamma[gamma_dat[i, 1], gamma_dat[i, 2]] = gamma_dat[i, 3] * (0.75 + 0.5 * rand())
end
gamma[:, 5:7] .= gamma[:, 3]
# Parameter phi[j], purchasing cost [EUR/[refres]):
phi = zeros(Float64, J)
phi[1] = 0.05 * (0.75 + 0.5 * rand()) # biomass
phi[7] = 1E-3 * (0.75 + 0.5 * rand()) # CO2

# Parameter psi[j], discharging cost (EUR/[refres]):
psi = zeros(Float64, J)
psi[7] = 0.05 * (0.75 + 0.5 * rand()) # CO2
#psi[[10,12,20]] .= -1E-2
# -----------------------------------------------------------------------------
