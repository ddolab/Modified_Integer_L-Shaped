using Random
using Distributions

#= num_demand_sites = 75
num_processing_sites = 15
num_scenarios = 96
num_time_periods = 24 =#
umax=200

I=1:num_demand_sites
J=0:num_processing_sites
S=1:num_scenarios
M = 1:3
T=1:num_time_periods
K=Dict()
K[1]=1:4*num_processing_sites 
K[2]=1:2*num_processing_sites
K[3]=1:num_processing_sites

Random.seed!(rng_seed)




vmax=Dict()
for j in J
    vmax[j,1]=4
    vmax[j,2]=1
end

u= Dict() #upper bound on production rate
u[1]=50
u[2]=100
u[3]=200

function calculate_distances(num_demand_sites)
    
    # Step 1: Generate num_demand_sites points with random coordinates in a 100x100 unit plane
    points = [(rand(Uniform(0,100)), rand(Uniform(0,100))) for _ in 1:num_demand_sites]

    # Step 2: Create a data structure (dictionary) to store distances between points
    distances = Dict()

    # Step 3: Calculate distances between all pairs of points
    for i in 1:length(points)
        for j in i+1:length(points)
            # Calculate Euclidean distance
            dist = sqrt((points[i][1] - points[j][1])^2 + (points[i][2] - points[j][2])^2)
            # Store the distance in the dictionary
            distances[(i, j)] = dist
            distances[(j, i)] = dist
        end
    end
    
    return distances, points
end
distances, points = calculate_distances(num_demand_sites)

# Select num_processing_sites random unique points from the points array
function select_random_points(points, num_processing_sites)
    return sample(points, num_processing_sites, replace=false)
end
production_centers = Dict()
production_centers[0] = (50, 50)  # Add production center with key 0 at location (50, 50)
for (index, point) in enumerate(select_random_points(points, num_processing_sites))
    production_centers[index] = point
end
# Calculate the distance from each point to each production center
function calculate_distances_to_production_centers(production_centers, points)
    distances_to_production_centers = Dict()
    for i in 1:length(points)
        for j in J
            dist = sqrt((production_centers[j][1] - points[i][1])^2 + (production_centers[j][2] - points[i][2])^2)
            distances_to_production_centers[(i, j)] = dist
        end
    end
    return distances_to_production_centers
end
distances_to_production_centers = calculate_distances_to_production_centers(production_centers, points)

#calculate the distance between each pair of production centers
distances_between_production_centers = Dict()
for i in J
    for j in J
        dist = sqrt((production_centers[i][1] - production_centers[j][1])^2 + (production_centers[i][2] - production_centers[j][2])^2)
        distances_between_production_centers[(i, j)] = dist
        distances_between_production_centers[(j, i)] = dist
    end
end


c = Dict()
product_transport_cost_per_distance = rand(Uniform(5,20))
operating_cost_per_unit_production=rand(Uniform(1,2))
for i in I
    for j in J
        for t in T
            c[(i, j, t)] = operating_cost_per_unit_production + distances_to_production_centers[(i, j)] * product_transport_cost_per_distance
        end
    end
end
demandRV=Uniform(0,num_processing_sites*umax/num_demand_sites)
d = Dict((i, t, s) => rand(demandRV) for i in I, t in T, s in S)
small_module_relocation_cost_per_distance = rand(Uniform(5,20))
h = Dict((j, jprime, m, t) =>  small_module_relocation_cost_per_distance * distances_between_production_centers[j, jprime] for j in J, jprime in J, m in [1], t in T)
relocation_cost_module_size_factor = rand(Uniform(2.5,3))
merge!(h, Dict((j, jprime, m, t) => h[j,jprime,m-1,t]*relocation_cost_module_size_factor for j in J, jprime in J, m in [2], t in T))
merge!(h, Dict((j, jprime, m, t) => h[j,jprime,m-1,t]*relocation_cost_module_size_factor for j in J, jprime in J, m in [3], t in T))
penalty_cost = Dict((i, t) => 1e6 for i in I, t in T)
g_small=rand(Uniform(1e5,3e5))
capital_cost_capacity_exponent=rand(Uniform(0.8,0.95))
g= Dict()
g[1]= g_small
g[2]= g_small*2^capital_cost_capacity_exponent#100 vs 50 production capacity
g[3]= g_small*4^capital_cost_capacity_exponent#50 vs 200 production capacity

p=[1/length(S) for s in S]