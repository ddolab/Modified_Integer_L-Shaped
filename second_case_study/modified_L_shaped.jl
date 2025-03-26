

using Distributed
rmprocs(workers())
addprocs(8)
@show workers()
@everywhere using JuMP, Gurobi
@everywhere using NLPModelsJuMP, NLPModels







##########JuMP specific code ###############



# Declare variables


# see https://github.com/jump-dev/JuMP.jl/blob/master/src/lp_sensitivity2.jl#L302-L412
# Converts an LP model into a standard form. This is used to find F_vec, G_vec, and d_vec automatically for the given problem
# While it might be possible to do this for any LP, it currently only works if all constraints are <= constraints

#calculate hi, Fi, Gi, and di for all i (all scenarios/subproblems)
num_scenarios = S
S=1:num_scenarios
S_avg=1:1 #scenarios for master problem
S_sp=1:num_scenarios #scenarios for subproblems
G_vec = Vector{Any}(undef, length(S))
F_vec = Vector{Any}(undef, length(S))
h_vec = Vector{Any}(undef, length(S))
d_vec = Vector{Any}(undef, length(S))
#set up model for subproblem of scenario i (same Fi,Gi,di as subproblem i)
@everywhere const env = Gurobi.Env(;output_flag=0)
@everywhere Gurobi.GRBsetintparam(env, "OutputFlag", 0)#supress output of dual separation problems
@everywhere function setup_DSP_1(s,problem_data)
    global env
    (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)=problem_data
    # see https://github.com/jump-dev/JuMP.jl/blob/master/src/lp_sensitivity2.jl#L302-L412
    # Converts an LP model into a standard form. This is used to find F_vec, G_vec, and d_vec automatically for the given problem
    # While it might be possible to do this for any LP, it currently only works if all constraints are <= constraints

    #set up model that includes constraints for scenario i (same Fi,Gi,di as subproblem i)
    modelMP = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
    ## first stage state vars
    #= @variables modelMP begin
        x[i=1:I], Bin  # 1 if process built
        x_[j=1:J], Bin  # 1 if resource inventory built
        #new variables 
        u[i=1:I, k=1:K], Bin
        u_[j=1:J, k=1:K], Bin
    end =#
    @variables modelMP begin
        x[i=1:I] # 1 if process built
        x_[j=1:J]  # 1 if resource inventory built
        #new variables 
        u[i=1:I, k=1:K]
        u_[j=1:J, k=1:K]
    end
    X = vcat(x...,x_...,u...,u_...) #vector of all first stage state variables

    #= @variables modelMP begin
        w[i=1:I, l=1:1+L_i[i]], Bin  # 1 if capacity lies in inteval
        V[i=1:I] >= 0  # capital cost
        lambda[i=1:I, l=1:1+L_i[i]] >= 0  # coefficient for piecewise-linear approximation
        Cmaster[i=1:I] >= 0
        C_master[j=1:J] >= 0
    end
    Z = vcat(w...,V...,lambda...,Cmaster...,C_master...) #vector of all first stage stage variables =#
        
        
        ## second stage
        # Variables:
    @variables modelMP begin
        B[j=1:J, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # resource input
        C[i=1:I, s=[s]] >= 0  # process capacity
        C_[j=1:J, s=[s]] >= 0  # inventory capacity
        P[i=1:I, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption of reference resource
        P_[i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption in operating mode
        Q[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # inventory
        Q_[s=[s], j=1:J, h=1:H]  # inventory carried over to the next season
        Sale[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # sales or waste
        y[s=[s], i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if mode selected
        z[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if switching from mode m1 to mode m2 at t
        Dslack[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # slack for unmet demand
    end
        
        # Network design constraints:
        @constraints modelMP begin
        ndconstr1[s=[s], i=1:I], C[i, s]/1e3 <= Cmax[i] * x[i]/1e3
        ndconstr2[s=[s], j=1:J], C_[j, s]/1e3 <= C_max[j] * x_[j]/1e3
        ndconstr5[s=[s], i=1:I], C[i, s]/1e3 == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)/1e3
        ndconstr6[s=[s], j=1:J], C_[j, s]/1e3 == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)/1e3
        end
        
        # Mass balance constraitns:
        @constraints modelMP begin
        mbconstr1[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] == (1 - epsilon[j, h]) * Q[s, j, h, t-1] + sum(rho[i, j, h] * P[i, h, t, s] for i = 1:I) + B[j, h, t, s] - Sale[s, j, h, t]
        mbconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] <= eta[i, h, t, s] * C[i, s]
        mbconstr3[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] <= C_[j, s]
        mbconstr4[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], B[j, h, t, s] <= Bmax[j, h, t]
        mbconstr5[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J_[j] == 1], Sale[s, j, h, t] + Dslack[s, j, h, t] >= D[j, h, t]
        #mbconstr5[j=1:J,h=1:H,t=θ_max+1:θ_max+T_h[h]; J_[j]==1], Sale[s, j,h,t] >= D[j,h,t]
        mbconstr6[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J0[j] == 1], Sale[s, j, h, t] == 0
        end
        
        # Mode-based operation:
        @constraints modelMP begin
        moconstr1[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(y[s, i, m, h, t] for m = 1:M if M_i[i, m] == 1) == x[i]
        moconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] == sum(P_[i, m, h, t, s] for m = 1:M if M_i[i, m] == 1)
        moconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(CCmin[i, m] * y[s, i, m, h, t]) <= 1e-3*(P_[i, m, h, t, s])
        moconstr4[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s]) <= 1e-3*(CCmax[i, m] * y[s, i, m, h, t])
        end
        
        # Transition constraints:
        @constraints modelMP begin
        trconstr1[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(-Delta[i, m] - bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1])) <= 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s])
        trconstr2[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s]) <= 1e-3*(Delta[i, m] + bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1]))
        trconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], sum(z[s, i, m1, m, h, t-1] for m1 = 1:M if TR_i[i, m1, m] == 1) - sum(z[s, i, m, m2, h, t-1] for m2 = 1:M if TR_i[i, m, m2] == 1) == y[s, i, m, h, t] - y[s, i, m, h, t-1]
        trconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; TR_i[i, m1, m2] == 1], y[s, i, m2, h, t] >= sum(z[s, i, m1, m2, h, t-t_prime] for t_prime = 1:theta[i, m1, m2])
        trconstr5[s=[s], i=1:I, m1=1:M, m2=1:M, m3=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; SQ_i[i, m1, m2, m3] == 1], z[s, i, m1, m2, h, t-theta_[i, m1, m2, m3]] == z[s, i, m2, m3, h, t]
        end
        #= @constraints modelMP begin
            [s=[s], h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(z[s,i,m1,m2,h,t] for i=1:I, m1=1:M, m2=1:M) <= 1
          end =#
        
        # Continuity equations:
        @constraints modelMP begin
        ctconstr1[s=[s], i=1:I, m=1:M, h=1:H; M_i[i, m] == 1], y[s, i, m, h, θ_max] == y[s, i, m, h, θ_max+T_h[h]]
        ctconstr2[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t] == z[s, i, m1, m2, h, t+T_h[h]]
        ctconstr3[s=[s], i=1:I, m=1:M, h=1:H-1; M_i[i, m] == 1], y[s, i, m, h, θ_max+T_h[h]] == y[s, i, m, h+1, θ_max]
        ctconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H-1, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t+T_h[h]] == z[s, i, m1, m2, h+1, t]
        ctconstr5[s=[s], j=1:J, h=1:H], Q_[s, j, h] == Q[s, j, h, θ_max+T_h[h]] - Q[s, j, h, θ_max]
        ctconstr6[s=[s], j=1:J, h=1:H-1], Q[s, j, h, θ_max] + n[h] * Q_[s, j, h] == Q[s, j, h+1, θ_max]
        ctconstr7[s=[s], j=1:J], Q[s, j, H, θ_max] + n[H] * Q_[s, j, H] >= Q[s, j, 1, θ_max]
        end
        
        
        
        # Objective function:
        @objective(modelMP, Max, -p_s*(sum(n[h] * (sum(delta[i, m] * y[s, i, m, h, t] + gamma[i, m] * P_[i, m, h, t, s] for i = 1:I, m = 1:M if M_i[i, m] == 1)
                                                + sum(phi[j] * B[j, h, t, s] for j = 1:J) + sum(psi[j] * Sale[s, j, h, t] for j = 1:J)) for h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s])
                                + sum(1E4 * Dslack[s, j, h, t] for j = 1:J, h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s]))/1e3)
        
                
    
    
    model=modelMP
    vars = JuMP.all_variables(model)
    JuMP.relax_integrality(model) # Relax the integrality of the binary variables
    for variable in JuMP.all_variables(model)
        # Check if the variable has a lower bound before accessing it
        if JuMP.has_lower_bound(variable)
        lb = JuMP.lower_bound(variable)
        JuMP.@constraint(model, variable >= lb)
        JuMP.delete_lower_bound(variable) # Correctly remove lower bound
        end
    
        # Check if the variable has an upper bound before accessing it
        if JuMP.has_upper_bound(variable)
        ub = JuMP.upper_bound(variable)
        JuMP.@constraint(model, variable <= ub)
        JuMP.delete_upper_bound(variable) # Correctly remove upper bound
        end
    end
    nlp = MathOptNLPModel(model) # NLPModelsJuMP enters here
    x = zeros(nlp.meta.nvar)
    c = grad(nlp, x) # = c = [1, 2, 3, -2, -4, -8]
    A = jac(nlp, x) # = A = 5x6 12 entries
    constraints = cons(nlp, x) # = g = zeros(5)
    nlp.meta.lcon, nlp.meta.ucon # l <= Ax + g <= u = ([0,0,-Inf,-Inf,-Inf], [0, 0, 1, 1, 1])
    # constraint indexes for each situation
    # [1, 2], [], [3, 4, 5], []
    nlp.meta.jfix, nlp.meta.jlow, nlp.meta.jupp, nlp.meta.jrng
    x_vector = filter(var -> var in X, vars)
    @assert X==x_vector
    y_vector = filter(var -> !(var in X), vars)
    Afix = A[nlp.meta.jfix, :]
    Alow = A[nlp.meta.jlow, :]
    Aupp = A[nlp.meta.jupp, :]
    ub = nlp.meta.ucon[nlp.meta.jupp]
    lb = nlp.meta.lcon[nlp.meta.jlow]
    b = nlp.meta.ucon[nlp.meta.jfix]
    #= The following is another statement of the constraints
    Aupp*vars<=ub (or -Aupp*vars>=-ub)
    Alow*vars>=lb (or -Alow*vars<=-lb)
    Afix*vars==b

    equivalently,
    [Aupp;   * vars <= [ub;
    -Alow;]               -lb;]
    Afix*vars==b
    =#
    Aineq = vcat(Aupp, -Alow,Afix,-Afix)
    bineq = vcat(ub, -lb,b,-b)
    F_vec_s = Aineq[:, findall(var -> var in x_vector, vars)]
    G_vec_s = Aineq[:, findall(var -> var in y_vector, vars)]
    d_vec_s=bineq
    part_of_c_for_y = c[findall(var -> var in y_vector, vars)]
    h_vec_s = part_of_c_for_y
    return F_vec_s,G_vec_s,d_vec_s,h_vec_s
end
results=let problem_data= (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)
    wp = CachingPool(workers())
    pmap(s->setup_DSP_1(s,problem_data),wp,S; retry_delays= ExponentialBackOff(n = 3))
end
F_vec=[results[i][1] for i in eachindex(results)]
G_vec=[results[i][2] for i in eachindex(results)]
d_vec=[results[i][3] for i in eachindex(results)]
h_vec=[results[i][4] for i in eachindex(results)]





UB_model=Vector{Any}(undef,num_scenarios)

#Calculate U (upper bound on subproblem cost)
function getUB()
    U = Vector{Any}(undef, length(S))

    U=let problem_data= (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)
        wp = CachingPool(workers())
    pmap(s->do_UB_optimization(s,problem_data),wp,S; retry_delays= ExponentialBackOff(n = 3))
    #continue here, make this a distributed function that solves the model. SEE SHARED ARRAYS LOOKS PROMISING. share only the model variables, solve them in workers, then get results undistributed
    end
    return U
end

@everywhere function do_UB_optimization(s,problem_data)
    (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)=problem_data
    #print(UB_model_s)

    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    function setup_UB_optimization(s)
        modelMP = direct_model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
        ## first stage state vars
            @variables modelMP begin
                x[i=1:I], Bin  # 1 if process built
                x_[j=1:J], Bin  # 1 if resource inventory built
                #new variables 
                u[i=1:I, k=1:K], Bin
                u_[j=1:J, k=1:K], Bin
            end
            #= @variables modelMP begin
                x[i=1:I] # 1 if process built
                x_[j=1:J]  # 1 if resource inventory built
                #new variables 
                u[i=1:I, k=1:K]
                u_[j=1:J, k=1:K]
            end =#
            X = vcat(x...,x_...,u...,u_...) #vector of all first stage state variables

            @variables modelMP begin
                w[i=1:I, l=1:1+L_i[i]], Bin  # 1 if capacity lies in inteval
                V[i=1:I] >= 0  # capital cost
                lambda[i=1:I, l=1:1+L_i[i]] >= 0  # coefficient for piecewise-linear approximation
                Cmaster[i=1:I] >= 0
                C_master[j=1:J] >= 0
            end
            Z = vcat(w...,V...,lambda...,Cmaster...,C_master...) #vector of all first stage stage variables
          
          
          ## second stage
          # Variables:
        @variables modelMP begin
            B[j=1:J, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # resource input
            C[i=1:I, s=[s]] >= 0  # process capacity
            C_[j=1:J, s=[s]] >= 0  # inventory capacity
            P[i=1:I, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption of reference resource
            P_[i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption in operating mode
            Q[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # inventory
            Q_[s=[s], j=1:J, h=1:H]  # inventory carried over to the next season
            Sale[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # sales or waste
            y[s=[s], i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if mode selected
            z[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if switching from mode m1 to mode m2 at t
            Dslack[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # slack for unmet demand
        end
          
          # Network design constraints:
          @constraints modelMP begin
            ndconstr1[s=[s], i=1:I], C[i, s]/1e3 <= Cmax[i] * x[i]/1e3
            ndconstr2[s=[s], j=1:J], C_[j, s]/1e3 <= C_max[j] * x_[j]/1e3
            ndconstr5[s=[s], i=1:I], C[i, s]/1e3 == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)/1e3
            ndconstr6[s=[s], j=1:J], C_[j, s]/1e3 == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)/1e3
          end
          
          # Mass balance constraitns:
          @constraints modelMP begin
            mbconstr1[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] == (1 - epsilon[j, h]) * Q[s, j, h, t-1] + sum(rho[i, j, h] * P[i, h, t, s] for i = 1:I) + B[j, h, t, s] - Sale[s, j, h, t]
            mbconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] <= eta[i, h, t, s] * C[i, s]
            mbconstr3[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] <= C_[j, s]
            mbconstr4[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], B[j, h, t, s] <= Bmax[j, h, t]
            mbconstr5[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J_[j] == 1], Sale[s, j, h, t] + Dslack[s, j, h, t] >= D[j, h, t]
            #mbconstr5[j=1:J,h=1:H,t=θ_max+1:θ_max+T_h[h]; J_[j]==1], Sale[s, j,h,t] >= D[j,h,t]
            mbconstr6[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J0[j] == 1], Sale[s, j, h, t] == 0
          end
          
          # Mode-based operation:
          @constraints modelMP begin
            moconstr1[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(y[s, i, m, h, t] for m = 1:M if M_i[i, m] == 1) == x[i]
            moconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] == sum(P_[i, m, h, t, s] for m = 1:M if M_i[i, m] == 1)
            moconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(CCmin[i, m] * y[s, i, m, h, t]) <= 1e-3*(P_[i, m, h, t, s])
            moconstr4[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s]) <= 1e-3*(CCmax[i, m] * y[s, i, m, h, t])
          end
          
          # Transition constraints:
          @constraints modelMP begin
            trconstr1[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(-Delta[i, m] - bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1])) <= 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s])
            trconstr2[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s]) <= 1e-3*(Delta[i, m] + bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1]))
            trconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], sum(z[s, i, m1, m, h, t-1] for m1 = 1:M if TR_i[i, m1, m] == 1) - sum(z[s, i, m, m2, h, t-1] for m2 = 1:M if TR_i[i, m, m2] == 1) == y[s, i, m, h, t] - y[s, i, m, h, t-1]
            trconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; TR_i[i, m1, m2] == 1], y[s, i, m2, h, t] >= sum(z[s, i, m1, m2, h, t-t_prime] for t_prime = 1:theta[i, m1, m2])
            trconstr5[s=[s], i=1:I, m1=1:M, m2=1:M, m3=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; SQ_i[i, m1, m2, m3] == 1], z[s, i, m1, m2, h, t-theta_[i, m1, m2, m3]] == z[s, i, m2, m3, h, t]
          end
          #= @constraints modelMP begin
            [s=[s], h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(z[s,i,m1,m2,h,t] for i=1:I, m1=1:M, m2=1:M) <= 1
          end =#
          
          # Continuity equations:
          @constraints modelMP begin
            ctconstr1[s=[s], i=1:I, m=1:M, h=1:H; M_i[i, m] == 1], y[s, i, m, h, θ_max] == y[s, i, m, h, θ_max+T_h[h]]
            ctconstr2[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t] == z[s, i, m1, m2, h, t+T_h[h]]
            ctconstr3[s=[s], i=1:I, m=1:M, h=1:H-1; M_i[i, m] == 1], y[s, i, m, h, θ_max+T_h[h]] == y[s, i, m, h+1, θ_max]
            ctconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H-1, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t+T_h[h]] == z[s, i, m1, m2, h+1, t]
            ctconstr5[s=[s], j=1:J, h=1:H], Q_[s, j, h] == Q[s, j, h, θ_max+T_h[h]] - Q[s, j, h, θ_max]
            ctconstr6[s=[s], j=1:J, h=1:H-1], Q[s, j, h, θ_max] + n[h] * Q_[s, j, h] == Q[s, j, h+1, θ_max]
            ctconstr7[s=[s], j=1:J], Q[s, j, H, θ_max] + n[H] * Q_[s, j, H] >= Q[s, j, 1, θ_max]
          end
          
          
          
          # Objective function:
          @objective(modelMP, Max, -p_s*(sum(n[h] * (sum(delta[i, m] * y[s, i, m, h, t] + gamma[i, m] * P_[i, m, h, t, s] for i = 1:I, m = 1:M if M_i[i, m] == 1)
                                                 + sum(phi[j] * B[j, h, t, s] for j = 1:J) + sum(psi[j] * Sale[s, j, h, t] for j = 1:J)) for h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s])
                                   + sum(1E4 * Dslack[s, j, h, t] for j = 1:J, h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s]))/1e3)
          
        return modelMP
    end


    global env
    UB_model_s=setup_UB_optimization(s)
    set_optimizer_attribute(UB_model_s, "MIPGap", 0.001)
    set_optimizer_attribute(UB_model_s, "TimeLimit", 250)
    set_optimizer_attribute(UB_model_s, "Threads",8)
    optimize!(UB_model_s)

    #= #Memory related code
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end =#

    return objective_bound(UB_model_s)
end
@time U=getUB()












RMP = direct_model( Gurobi.Optimizer())
#first stage 
## first stage state vars
@variables RMP begin
    x[i=1:I], Bin  # 1 if process built
    x_[j=1:J], Bin  # 1 if resource inventory built
    #new variables 
    u[i=1:I, k=1:K], Bin
    u_[j=1:J, k=1:K], Bin
end
#= @variables modelMP begin
    x[i=1:I] # 1 if process built
    x_[j=1:J]  # 1 if resource inventory built
    #new variables 
    u[i=1:I, k=1:K]
    u_[j=1:J, k=1:K]
end =#
X = vcat(x...,x_...,u...,u_...) #vector of all first stage state variables

@variables RMP begin
    w[i=1:I, l=1:1+L_i[i]], Bin  # 1 if capacity lies in inteval
    V[i=1:I] >= 0  # capital cost
    lambda[i=1:I, l=1:1+L_i[i]] >= 0  # coefficient for piecewise-linear approximation
    Cmaster[i=1:I] >= 0
    C_master[j=1:J] >= 0
end
Z = vcat(w...,V...,lambda...,Cmaster...,C_master...) #vector of all first stage stage variables

# Network design constraints:
@constraints RMP begin
    ndconstr3[i=1:I], sum(u[i, k] for k = 1:K) == x[i]
    ndconstr4[j=1:J], sum(u_[j, k] for k = 1:K) == x_[j]
    ndconstr7[i=1:I], Cmaster[i]/1e3 == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)/1e3
    ndconstr8[j=1:J], C_master[j]/1e3 == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)/1e3
end

#Extra Network design constraints:
@constraints RMP begin
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

# Piecewise-linear cost functions:
@constraints RMP begin
    lcconstr1[i=1:I], Cmaster[i]/1e3 == sum(lambda[i, l] * (Cl[i, l-1] - Cl[i, l]) + Cl[i, l] * w[i, l] for l = 2:1+L_i[i])/1e3
    lcconstr2[i=1:I], V[i]/1e3 == sum(lambda[i, l] * (Vl[i, l-1] - Vl[i, l]) + Vl[i, l] * w[i, l] for l = 2:1+L_i[i])/1e3
    lcconstr3[i=1:I, l=2:1+L_i[i]], lambda[i, l] <= w[i, l]
    lcconstr4[i=1:I], sum(w[i, l] for l = 2:1+L_i[i]) == x[i]
end
#= 
#second stage for mp scenarios
## second stage
# Variables:
@variables RMP begin
    B[j=1:J, h=1:H, t=1:θ_max+T_h[h], s=S_avg] >= 0  # resource input
    C[i=1:I, s=S_avg] >= 0  # process capacity
    C_[j=1:J, s=S_avg] >= 0  # inventory capacity
    P[i=1:I, h=1:H, t=1:θ_max+T_h[h], s=S_avg] >= 0  # production or consumption of reference resource
    P_[i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h], s=S_avg] >= 0  # production or consumption in operating mode
    Q[s=S_avg, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # inventory
    Q_[s=S_avg, j=1:J, h=1:H]  # inventory carried over to the next season
    Sale[s=S_avg, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # sales or waste
    0 <= y[s=S_avg, i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h]] <=1  # 1 if mode selected
    0<= z[s=S_avg, i=1:I, m1=1:M, m2=1:M, h=1:H, t=1:θ_max+T_h[h]] <=1  # 1 if switching from mode m1 to mode m2 at t
    Dslack[s=S_avg, j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # slack for unmet demand (penalty applied)
end
Y_avg = [vcat(B[:,:,:,s]...,C[:,s]...,C_[:,s]...,P[:,:,:,s]...,P_[:,:,:,:,s]...,Q[s,:,:,:]...,Q_[s,:,:]...,Sale[s,:,:,:]...,y[s,:,:,:,:]...,z[s,:,:,:,:,:]...,Dslack[s,:,:,:]...) for s in S_avg]
# Network design constraints:
@constraints RMP begin
    ndconstr1[s=S_avg, i=1:I], C[i, s]/1e3 <= Cmax[i] * x[i]/1e3
    ndconstr2[s=S_avg, j=1:J], C_[j, s]/1e3 <= C_max[j] * x_[j]/1e3
    ndconstr5[s=S_avg, i=1:I], C[i, s]/1e3 == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)/1e3
    ndconstr6[s=S_avg, j=1:J], C_[j, s]/1e3 == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)/1e3
end
    
    # Mass balance constraitns:
@constraints RMP begin
    mbconstr1[s=S_avg, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] == (1 - epsilon[j, h]) * Q[s, j, h, t-1] + sum(rho[i, j, h] * P[i, h, t, s] for i = 1:I) + B[j, h, t, s] - Sale[s, j, h, t]
    mbconstr2[s=S_avg, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] <= eta_avg[i, h, t] * C[i, s]
    mbconstr3[s=S_avg, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] <= C_[j, s]
    mbconstr4[s=S_avg, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], B[j, h, t, s] <= Bmax[j, h, t]
    mbconstr5[s=S_avg, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J_[j] == 1], Sale[s, j, h, t] + Dslack[s, j, h, t] >= D[j, h, t]
    #mbconstr5[j=1:J,h=1:H,t=θ_max+1:θ_max+T_h[h]; J_[j]==1], Sale[s, j,h,t] >= D[j,h,t]
    mbconstr6[s=S_avg, j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J0[j] == 1], Sale[s, j, h, t] == 0
end
    
    # Mode-based operation:
@constraints RMP begin
    moconstr1[s=S_avg, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(y[s, i, m, h, t] for m = 1:M if M_i[i, m] == 1) == x[i]
    moconstr2[s=S_avg, i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] == sum(P_[i, m, h, t, s] for m = 1:M if M_i[i, m] == 1)
    moconstr3[s=S_avg, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(CCmin[i, m] * y[s, i, m, h, t]) <= 1e-3*(P_[i, m, h, t, s])
    moconstr4[s=S_avg, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s]) <= 1e-3*(CCmax[i, m] * y[s, i, m, h, t])
end
    
    # Transition constraints:
@constraints RMP begin
    trconstr1[s=S_avg, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(-Delta[i, m] - bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1])) <= 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s])
    trconstr2[s=S_avg, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s]) <= 1e-3*(Delta[i, m] + bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1]))
    trconstr3[s=S_avg, i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], sum(z[s, i, m1, m, h, t-1] for m1 = 1:M if TR_i[i, m1, m] == 1) - sum(z[s, i, m, m2, h, t-1] for m2 = 1:M if TR_i[i, m, m2] == 1) == y[s, i, m, h, t] - y[s, i, m, h, t-1]
    trconstr4[s=S_avg, i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; TR_i[i, m1, m2] == 1], y[s, i, m2, h, t] >= sum(z[s, i, m1, m2, h, t-t_prime] for t_prime = 1:theta[i, m1, m2])
    trconstr5[s=S_avg, i=1:I, m1=1:M, m2=1:M, m3=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; SQ_i[i, m1, m2, m3] == 1], z[s, i, m1, m2, h, t-theta_[i, m1, m2, m3]] == z[s, i, m2, m3, h, t]
end
#= @constraints RMP begin
    [s=S_avg, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(z[s,i,m1,m2,h,t] for i=1:I, m1=1:M, m2=1:M) <= 1
  end =#

    # Continuity equations:
@constraints RMP begin
    ctconstr1[s=S_avg, i=1:I, m=1:M, h=1:H; M_i[i, m] == 1], y[s, i, m, h, θ_max] == y[s, i, m, h, θ_max+T_h[h]]
    ctconstr2[s=S_avg, i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t] == z[s, i, m1, m2, h, t+T_h[h]]
    ctconstr3[s=S_avg, i=1:I, m=1:M, h=1:H-1; M_i[i, m] == 1], y[s, i, m, h, θ_max+T_h[h]] == y[s, i, m, h+1, θ_max]
    ctconstr4[s=S_avg, i=1:I, m1=1:M, m2=1:M, h=1:H-1, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t+T_h[h]] == z[s, i, m1, m2, h+1, t]
    ctconstr5[s=S_avg, j=1:J, h=1:H], Q_[s, j, h] == Q[s, j, h, θ_max+T_h[h]] - Q[s, j, h, θ_max]
    ctconstr6[s=S_avg, j=1:J, h=1:H-1], Q[s, j, h, θ_max] + n[h] * Q_[s, j, h] == Q[s, j, h+1, θ_max]
    ctconstr7[s=S_avg, j=1:J], Q[s, j, H, θ_max] + n[H] * Q_[s, j, H] >= Q[s, j, 1, θ_max]
end

 =#
    

@variable(RMP, η[s=S_sp] .<= U[s]) #upper bound placed on ηi so problem is always bounded
#= @constraint(RMP, sum(η[i] for i in S_sp) <= -(sum(n[h] * (sum(delta[i, m] * y[s, i, m, h, t] + gamma[i, m] * P_[i, m, h, t, s] for i = 1:I, m = 1:M if M_i[i, m] == 1)
+ sum(phi[j] * B[j, h, t, s] for j = 1:J) + sum(psi[j] * Sale[s, j, h, t] for j = 1:J)) for h = 1:H, t = θ_max+1:θ_max+T_h[h], s = S_avg)
+ sum(1E4 * Dslack[s, j, h, t] for j = 1:J, h = 1:H, t = θ_max+1:θ_max+T_h[h], s = S_avg))) =#

@objective(RMP, Max, sum(η[i] for i in S_sp)  - (sum(V[i] for i = 1:I) + sum(alpha[j] * x_[j] + beta[j] * C_master[j] for j = 1:J))/1e3)
RMP_MIPgap=0.005
RMP_time_limit=3600*4










#To see how callbacks work with JuMP and Gurobi.jl: 
# https://github.com/jump-dev/Gurobi.jl#callbacks
 
cb_calls = Cint[]
Normal_cuts_added = Vector{Any}(undef, 0)
No_good_cuts_added = Vector{Any}(undef, 0)
DSP_times = Vector{Any}(undef, 0)
ISP_times = Vector{Any}(undef, 0)
ISP_s_times = Vector{Any}(undef, 0)
X_list_int_sol = Vector{Any}(undef, 0)
X_list_dsp_sol = Vector{Any}(undef, 0)
num_of_solutions_for_ISP=Vector{Any}(undef, 0)
feasible_mp_sol_list=Vector{Any}(undef, 0)
objval_list_dsp_sol=Vector{Any}(undef, 0)
time_limits=Vector{Any}(undef, 0)
ISP_solutions=Vector{Any}(nothing, length(S))
savedVBases=Vector{Any}(nothing, length(S))
savedCBases=Vector{Any}(nothing, length(S))
objval = Vector{Any}(undef, length(S))
DSP_objval = Dict()
bestbound = Vector{Any}(undef, length(S))
incumbent_objval=-Inf
just_improved_incumbent=false

ISP_optgaps=[0.10,0.01,0.00]
time_lim_initial=600
function my_callback_function(cb_data, cb_where::Cint)
    # You can reference variables outside the function as normal
    global DSP_times, ISP_times, Normal_cuts_added, No_good_cuts_added, X, η, incumbent_objval,num_of_solutions_for_ISP, ISP_solutions, time_limits, ISP_s_times, objval, DSP_objval, just_improved_incumbent,time_lim_initial
    push!(cb_calls, cb_where)
    # You can select where the callback is run
    if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
        return
    end
    if cb_where == GRB_CB_MIPNODE
        resultP = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
        if resultP[] != GRB_OPTIMAL
            return  # Solution is something other than optimal.
        end
        resultP2 = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP2)
        if resultP2[] != 0   #if the current (maybe fractional) nodal solution is not from the root node (integer feasible solutions will be checked in MIPSOL callback)  
            return 
        end
    end
    currenttime= time()-start_time
    #GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME,currenttime)
#=     if currenttime > RMP_time_limit
        print("trying to terminate RMP")
        GRBterminate(backend(RMP))
    end
 =#

    cuts_added = 0
    NG_cuts_added = 0

    # Get relaxed solution at current node
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    X_val = callback_value.(cb_data, X)
    η_val = callback_value.(cb_data, η)
    mp_sol = vcat(X_val..., η_val...)
    #check if the submitted mp solution has been found to be feasible  
    if any(all(j .- 1e-8 .<= mp_sol .<= j .+ 1e-8) for j in feasible_mp_sol_list) #  
        return
    end

    ##Memory related code
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    if any(all(j .- 1e-8 .<= X_val .<= j .+ 1e-8) for j in X_list_dsp_sol) #if the dsp has been solved for this X before
        index_of_old_dsp_solve=findfirst(all(j .- 1e-8 .<= X_val .<= j .+ 1e-8) for j in X_list_dsp_sol)
        if all(collect(η_val) .<= objval_list_dsp_sol[index_of_old_dsp_solve] + 1e-3*abs.(objval_list_dsp_sol[index_of_old_dsp_solve])) #if the solution is not better than the one found before
            need_to_do_dsp=false
            print("skipped doing dsp ")
            flush(stdout)
        else
            need_to_do_dsp=true
            print("did dsp for same X again ")
            flush(stdout)
        end
    else
        need_to_do_dsp=true
    end

        if need_to_do_dsp 
            #solve dual separation problems
            DSP_time = @timed begin
                push!(X_list_dsp_sol, X_val)
                results=let (F_vec,G_vec,d_vec,h_vec,X_val)=(F_vec,G_vec,d_vec,h_vec,X_val)
                    wp = CachingPool(workers())
                    #get primal_status, v, objective value, and ustar in pmap'd_vec function
                    pmap(s->do_DSP_optimization(s,(F_vec[s],G_vec[s],d_vec[s],h_vec[s],X_val),savedVBases[s],savedCBases[s],RMP_time_limit,currenttime),wp,S_sp; retry_delays= ExponentialBackOff(n = 3)) #continue here, make this a distributed function that solves the model. SEE SHARED ARRAYS LOOKS PROMISING. share only the model variables, solve them in workers, then get results undistributed
                end
                results=Dict(s=>results[s] for s in S_sp)
                for s in S_sp
                    if results[s][1] == INFEASIBILITY_CERTIFICATE #the solution is an unbounded ray that proves infeasiblity of the dual
                        print("INFEASIBLE DSP \n")
                        #= println(X_val)
                        flush(stdout)
                        ν=  results[s][2]
                        con=@build_constraint(ν' * (d_vec[s] - F_vec[s] * X) >= 0)
                        MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                        cuts_added = 1 + cuts_added =#
                    elseif results[s][1] == FEASIBLE_POINT
                        DSP_objval[s] = results[s][3]
                        η_val[s]
                        if η_val[s] > results[s][3] + 1e-1
                            ustar =  results[s][2]
                            con= @build_constraint(ustar' * (d_vec[s] - F_vec[s] * X) >= η[s])
                            MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                            cuts_added = 1 + cuts_added
                        end
                    else
                        print(results[s][1])
                        flush(stdout)
                    end
                end
            end
            push!(objval_list_dsp_sol,collect(values(DSP_objval)))
            DSP_times = push!(DSP_times, DSP_time)
            Normal_cuts_added = push!(Normal_cuts_added, cuts_added)
            if cuts_added >= 1
                print("ABC")#added benders cut
                flush(stdout)
                return
            end
        end
        #if no cut was added and the solution is integer
    if all([(- 1e-8 <= j <= 1e-8)||(1-1e-8 <= j <= 1+1e-8) for j in X_val])
        current_ISP_index=findfirst(all(j .- 1e-8 .<= X_val .<= j .+ 1e-8) for j in X_list_int_sol)
            if isnothing(current_ISP_index)
                push!(X_list_int_sol, X_val)
                current_ISP_index=length(X_list_int_sol)
                push!(num_of_solutions_for_ISP,0)
                push!(time_limits,ones(length(S))*time_lim_initial)
            end
            #solve integer version of subproblems to generate no-good cuts
            while num_of_solutions_for_ISP[current_ISP_index] < length(ISP_optgaps)
                integer_SP_time = @timed begin
                    results=let problem_data= (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)
                        wp = CachingPool(workers())
                        pmap(s->do_ISP_optimization(s,problem_data,ISP_optgaps,num_of_solutions_for_ISP,current_ISP_index,RMP_time_limit,currenttime,time_limits,ISP_solutions[s],X_val),wp,S_sp; retry_delays= ExponentialBackOff(n = 3)) 
                    end
                    results=Dict(s=>results[s] for s in S_sp)
                    adjusted_num_of_solutions_for_ISP=num_of_solutions_for_ISP[current_ISP_index]
                    original_num_of_solutions_for_ISP=num_of_solutions_for_ISP[current_ISP_index]
                    for s in S_sp
                        objval[s] = results[s][1]
                        bestbound[s]=results[s][2]
                        if ISP_optgaps[num_of_solutions_for_ISP[current_ISP_index]+1]==0 && results[s][3]==TIME_LIMIT
                            print("entering time limit callback part")
                            adjusted_num_of_solutions_for_ISP=length(ISP_optgaps)-2 #this will cause the current set of subproblems to be run again with a gap of 0
                            time_limits[current_ISP_index][s]=time_limits[current_ISP_index][s]*5
                            print("need to run with gap of 0 again")
                            flush(stdout)
                        end
                        #print(termination_status(ISP_models[s]))
                        #print(relative_gap(ISP_models[s]))
                        if results[s][3] == INFEASIBLE
                            #build and submit no-good feasibility cut
                            print("\n infeasible ISP \n")
                            flush(stdout)
                            con = @build_constraint(
                                sum(X[j] for j in findall(x -> x <= 0.5, X_val))
                                +
                                sum(1 - X[j] for j in findall(x -> x > 0.5, X_val)) >= 1)
                            MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                            NG_cuts_added = NG_cuts_added + 1
                        elseif results[s][3] == OPTIMAL || results[s][3] == TIME_LIMIT
                            if bestbound[s] + 1e-9 < η_val[s]
                                con= @build_constraint(
                                    η[s] <= bestbound[s] + (U[s] - bestbound[s]) *
                                    (sum(X[j] for j in findall(x -> x<=0.5,X_val))
                                    +sum(1-X[j] for j in findall(x -> x>0.5,X_val))))
                                    #NGC_RHS=bestbound[s] + (U[s] - bestbound[s]) *(sum(X_val[j] for j in findall(x -> x<=0.5,X_val))+sum(1-X_val[j] for j in findall(x -> x>0.5,X_val))) 
                                    #@show NGC_RHS
                                    #@show con
                                    #@show X_val
                                MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                                NG_cuts_added = NG_cuts_added + 1
                                #print("added ng cut")
                            end
                        else
                            print(results[s][3])
                            flush(stdout)
                            error()
                        end
                    end
                    num_of_solutions_for_ISP[current_ISP_index]=adjusted_num_of_solutions_for_ISP
                    original_num_of_solutions_for_ISP#this lets the @timed macro above record the solution number for the isp solve
                end
                push!(ISP_times, integer_SP_time)
                push!(ISP_s_times,[results[s][4] for s in S_sp])
                num_of_solutions_for_ISP[current_ISP_index] +=1
                if NG_cuts_added >=1 #if a no good cut was added to remove the current mp solution
                    print("ANGC ")#added no good cutprint("ANGC ")
                    flush(stdout)
                    #check to see if incumbent can be improved using feasible solution found in callback
                    #= incumbent_objval= Ref{Cdouble}()
                    GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST,incumbent_objval) =# #CURRENTLY BUGGED DUE TO GUROBI BUG, SEE https://support.gurobi.com/hc/en-us/community/posts/13366202203281-Wrong-values-of-best-incumbent-objective-value
                    #print("I believe a NGC was added")
                    V_val = callback_value.(cb_data, V)
                    x__val = callback_value.(cb_data, x_)
                    C_master_val = callback_value.(cb_data, C_master)
                    ISP_objval= ( - (sum(V_val[i] for i = 1:I) + sum(alpha[j] * x__val[j] + beta[j] * C_master_val[j] for j = 1:J))/1e3+sum(objval[i] for i in S_sp))
                    #@show ISP_objval
                    #@show incumbent_objval
                    #@show objval[S_sp]
                    #@show bestbound[S_sp]
                    #@show U[S_sp]-bestbound[S_sp]
                    if ISP_objval > incumbent_objval
                        MOI.submit(RMP, MOI.HeuristicSolution(cb_data), vcat(X,η...), convert.(AbstractFloat, vcat(X_val,objval[S_sp])))
                        print("new incumbent found \n")
                        flush(stdout)
                        just_improved_incumbent=true
                        incumbent_objval=ISP_objval
                        push!(feasible_mp_sol_list,vcat(X_val,objval[S_sp]))
                    end
                    No_good_cuts_added = push!(No_good_cuts_added, NG_cuts_added)
                    return
                end
                
            end
            
            print("got to end of isp part (uh oh?)")
            flush(stdout)
            
            V_val = callback_value.(cb_data, V)
            x__val = callback_value.(cb_data, x_)
            C_master_val = callback_value.(cb_data, C_master)
            ISP_objval= ( - (sum(V_val[i] for i = 1:I) + sum(alpha[j] * x__val[j] + beta[j] * C_master_val[j] for j = 1:J))/1e3+sum(objval[i] for i in S_sp))
            if ISP_objval > incumbent_objval
                print("new incumbent found \n")
                flush(stdout)
                just_improved_incumbent=true
                incumbent_objval=ISP_objval
                push!(feasible_mp_sol_list,vcat(X_val,objval[S_sp]))
            end
            No_good_cuts_added = push!(No_good_cuts_added, NG_cuts_added)
        end

    return
end

@everywhere function do_DSP_optimization(s,problem_data,savedVBasis,savedCBasis,RMP_time_limit,currenttime)
    (F_s,G_s,d_s,h_s,X_val)=problem_data
    rounded_X_val=round.(X_val;digits=6)
    #print(UB_model_s)
    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end
    
    function setup_DSP_optimization(s)
        DSP_model_s = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0, "InfUnbdInfo" => 1))
        #DSP_model_s = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env2),  "InfUnbdInfo" => 1))
        @variable(DSP_model_s, u[1:length(d_s)] >= 0)
        #The G_vec, F_vec, h_vec, and d_vec were obtained from form of constraints with free decision variables
        #i.e., the subproblems have free decision variables 
        @constraint(DSP_model_s, G_s' * u .== h_s)
        @objective(DSP_model_s, Min, (d_s - F_s * rounded_X_val)' * u)
        return DSP_model_s
    end
    DSP_model_s=setup_DSP_optimization(s)
    
    
    #set_optimizer_attribute(DSP_model_s, "MemLimit", 1)
    set_optimizer_attribute(DSP_model_s, "Threads", 8)
    if currenttime > RMP_time_limit + 100
        print("trying to terminate DSP")
        flush(stdout)
        return ["DSP_ERROR", zeros(length(d_s)), 1E9, zeros(length(d_s))]
    end
    set_optimizer_attribute(DSP_model_s, "TimeLimit", RMP_time_limit+101-currenttime)
    #set_optimizer_attribute(DSP_model_s, "NumericFocus", 3)
    #set_optimizer_attribute(DSP_model_s, "FeasibilityTol", 1e-9)
    #set_optimizer_attribute(DSP_model_s, "IntFeasTol", 1e-9)
    #set_optimizer_attribute(DSP_model_s, "Method", 1)
    #set_optimizer_attribute(DSP_model_s, "ScaleFlag", 2)
    #= if !isnothing(savedVBasis) && !isnothing(savedCBasis)
        set_optimizer_attribute(DSP_model_s, "LPWarmStart", 2)
        allvars=all_variables(DSP_model_s)
        allcons=all_constraints(DSP_model_s;include_variable_in_set_constraints=false)
        MOI.set.(DSP_model_s,Gurobi.VariableAttribute("VBasis"),allvars,savedVBasis)
        MOI.set.(DSP_model_s,Gurobi.ConstraintAttribute("CBasis"),allcons,savedCBasis)
    end =#
    try
        optimize!(DSP_model_s)
        print(termination_status(DSP_model_s))
        flush(stdout)
        allvars=all_variables(DSP_model_s)
        savedVBasis=MOI.get.(DSP_model_s,Gurobi.VariableAttribute("VBasis"),allvars)
        allcons=all_constraints(DSP_model_s;include_variable_in_set_constraints=false)
        savedCBasis=MOI.get.(DSP_model_s,Gurobi.ConstraintAttribute("CBasis"),allcons)
        @spawnat 1 savedVBases[s]=savedVBasis
        @spawnat 1 savedCBases[s]=savedCBasis
            #print(InteractiveUtils.varinfo())
        #Memory related code
        GC.gc(true)
        GC.safepoint()
        try
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        catch
        end
    
        return [primal_status(DSP_model_s),value.(DSP_model_s[:u]),objective_value(DSP_model_s),value.(DSP_model_s[:u])]
    catch
        print("error in DSP")
        flush(stdout)
        return ["DSP_ERROR", zeros(length(d_s)), 1E9, zeros(length(d_s))]
        #turn LogToConsole to 1 to see the error
#=         MOI.set(DSP_model_s, MOI.RawOptimizerAttribute("OutputFlag"), 1)
        MOI.set(DSP_model_s, MOI.RawOptimizerAttribute("LogToConsole"), 1)
        compute_conflict!(DSP_model_s)
        error("error in DSP") =#
    end

end

@everywhere function do_ISP_optimization(s,problem_data,ISP_optgaps,num_of_solutions_for_ISP,current_ISP_index,RMP_time_limit,currenttime,current_time_limits,previous_solution_for_this_problem,X_val)
    (theta_,D,I,L_i,J,K,H,θ_max,T_h,M,Cmax,C_max,Cnameplate,C_nameplate,J_,J0,epsilon,rho,eta,Bmax,M_i,CCmin,CCmax,TR_i,SQ_i,Delta,bigM,theta,n,delta,gamma,phi,psi,p_s)=problem_data
    rounded_X_val=round.(X_val;digits=8)
    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    function setup_ISP_optimization(s)
        modelMP = direct_model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
        ## first stage
        ## first stage state vars
        @variables modelMP begin
            x[i=1:I], Bin  # 1 if process built
            x_[j=1:J], Bin  # 1 if resource inventory built
            #new variables 
            u[i=1:I, k=1:K], Bin
            u_[j=1:J, k=1:K], Bin
        end
        #= @variables modelMP begin
            x[i=1:I] # 1 if process built
            x_[j=1:J]  # 1 if resource inventory built
            #new variables 
            u[i=1:I, k=1:K]
            u_[j=1:J, k=1:K]
        end =#
        ispX = vcat(x...,x_...,u...,u_...) #vector of all first stage state variables

        #= @variables modelMP begin
            w[i=1:I, l=1:1+L_i[i]], Bin  # 1 if capacity lies in inteval
            V[i=1:I] >= 0  # capital cost
            lambda[i=1:I, l=1:1+L_i[i]] >= 0  # coefficient for piecewise-linear approximation
            Cmaster[i=1:I] >= 0
            C_master[j=1:J] >= 0
        end
        Z = vcat(w...,V...,lambda...,Cmaster...,C_master...) #vector of all first stage stage variables =#
        
        
        @constraint(modelMP, Xcon, ispX .== rounded_X_val)# This fixes the first stage variables to the RMP's solution
          
          
        ## second stage
        # Variables:
        @variables modelMP begin
            B[j=1:J, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # resource input
            C[i=1:I, s=[s]] >= 0  # process capacity
            C_[j=1:J, s=[s]] >= 0  # inventory capacity
            P[i=1:I, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption of reference resource
            P_[i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h], s=[s]] >= 0  # production or consumption in operating mode
            Q[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # inventory
            Q_[s=[s], j=1:J, h=1:H]  # inventory carried over to the next season
            Sale[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # sales or waste
            y[s=[s], i=1:I, m=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if mode selected
            z[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=1:θ_max+T_h[h]], Bin  # 1 if switching from mode m1 to mode m2 at t
            Dslack[s=[s], j=1:J, h=1:H, t=1:θ_max+T_h[h]] >= 0  # slack for unmet demand (penalty applied)
        end
          
        # Network design constraints:
        @constraints modelMP begin
            ndconstr1[s=[s], i=1:I], C[i, s]/1e3 <= Cmax[i] * x[i]/1e3
            ndconstr2[s=[s], j=1:J], C_[j, s]/1e3 <= C_max[j] * x_[j]/1e3
            ndconstr5[s=[s], i=1:I], C[i, s]/1e3 == sum(u[i, k] * Cnameplate[i, k] for k = 1:K)/1e3
            ndconstr6[s=[s], j=1:J], C_[j, s]/1e3 == sum(u_[j, k] * C_nameplate[j, k] for k = 1:K)/1e3
        end
          
          # Mass balance constraitns:
        @constraints modelMP begin
            mbconstr1[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] == (1 - epsilon[j, h]) * Q[s, j, h, t-1] + sum(rho[i, j, h] * P[i, h, t, s] for i = 1:I) + B[j, h, t, s] - Sale[s, j, h, t]
            mbconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] <= eta[i, h, t, s] * C[i, s]
            mbconstr3[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], Q[s, j, h, t] <= C_[j, s]
            mbconstr4[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]], B[j, h, t, s] <= Bmax[j, h, t]
            mbconstr5[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J_[j] == 1], Sale[s, j, h, t] + Dslack[s, j, h, t] >= D[j, h, t]
            #mbconstr5[j=1:J,h=1:H,t=θ_max+1:θ_max+T_h[h]; J_[j]==1], Sale[s, j,h,t] >= D[j,h,t]
            mbconstr6[s=[s], j=1:J, h=1:H, t=θ_max+1:θ_max+T_h[h]; J0[j] == 1], Sale[s, j, h, t] == 0
        end
          
          # Mode-based operation:
        @constraints modelMP begin
            moconstr1[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(y[s, i, m, h, t] for m = 1:M if M_i[i, m] == 1) == x[i]
            moconstr2[s=[s], i=1:I, h=1:H, t=θ_max+1:θ_max+T_h[h]], P[i, h, t, s] == sum(P_[i, m, h, t, s] for m = 1:M if M_i[i, m] == 1)
            moconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(CCmin[i, m] * y[s, i, m, h, t]) <= 1e-3*(P_[i, m, h, t, s])
            moconstr4[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s]) <= 1e-3*(CCmax[i, m] * y[s, i, m, h, t])
        end
          
          # Transition constraints:
        @constraints modelMP begin
            trconstr1[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(-Delta[i, m] - bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1])) <= 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s])
            trconstr2[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], 1e-3*(P_[i, m, h, t, s] - P_[i, m, h, t-1, s]) <= 1e-3*(Delta[i, m] + bigM * (2 - y[s, i, m, h, t] - y[s, i, m, h, t-1]))
            trconstr3[s=[s], i=1:I, m=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; M_i[i, m] == 1], sum(z[s, i, m1, m, h, t-1] for m1 = 1:M if TR_i[i, m1, m] == 1) - sum(z[s, i, m, m2, h, t-1] for m2 = 1:M if TR_i[i, m, m2] == 1) == y[s, i, m, h, t] - y[s, i, m, h, t-1]
            trconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; TR_i[i, m1, m2] == 1], y[s, i, m2, h, t] >= sum(z[s, i, m1, m2, h, t-t_prime] for t_prime = 1:theta[i, m1, m2])
            trconstr5[s=[s], i=1:I, m1=1:M, m2=1:M, m3=1:M, h=1:H, t=θ_max+1:θ_max+T_h[h]; SQ_i[i, m1, m2, m3] == 1], z[s, i, m1, m2, h, t-theta_[i, m1, m2, m3]] == z[s, i, m2, m3, h, t]
        end
        #= @constraints modelMP begin
            [s=[s], h=1:H, t=θ_max+1:θ_max+T_h[h]], sum(z[s,i,m1,m2,h,t] for i=1:I, m1=1:M, m2=1:M) <= 1
          end =#
        
          # Continuity equations:
        @constraints modelMP begin
            ctconstr1[s=[s], i=1:I, m=1:M, h=1:H; M_i[i, m] == 1], y[s, i, m, h, θ_max] == y[s, i, m, h, θ_max+T_h[h]]
            ctconstr2[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t] == z[s, i, m1, m2, h, t+T_h[h]]
            ctconstr3[s=[s], i=1:I, m=1:M, h=1:H-1; M_i[i, m] == 1], y[s, i, m, h, θ_max+T_h[h]] == y[s, i, m, h+1, θ_max]
            ctconstr4[s=[s], i=1:I, m1=1:M, m2=1:M, h=1:H-1, t=θ_max-maximum(reshape(theta[i, :, :], (1, M * M)))+1:θ_max-1; TR_i[i, m1, m2] == 1], z[s, i, m1, m2, h, t+T_h[h]] == z[s, i, m1, m2, h+1, t]
            ctconstr5[s=[s], j=1:J, h=1:H], Q_[s, j, h] == Q[s, j, h, θ_max+T_h[h]] - Q[s, j, h, θ_max]
            ctconstr6[s=[s], j=1:J, h=1:H-1], Q[s, j, h, θ_max] + n[h] * Q_[s, j, h] == Q[s, j, h+1, θ_max]
            ctconstr7[s=[s], j=1:J], Q[s, j, H, θ_max] + n[H] * Q_[s, j, H] >= Q[s, j, 1, θ_max]
        end
          
          
          
        # Objective function:
        @objective(modelMP, Max, -p_s*(sum(n[h] * (sum(delta[i, m] * y[s, i, m, h, t] + gamma[i, m] * P_[i, m, h, t, s] for i = 1:I, m = 1:M if M_i[i, m] == 1)
                                                + sum(phi[j] * B[j, h, t, s] for j = 1:J) + sum(psi[j] * Sale[s, j, h, t] for j = 1:J)) for h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s])
                                + sum(1E4 * Dslack[s, j, h, t] for j = 1:J, h = 1:H, t = θ_max+1:θ_max+T_h[h], s = [s]))/1e3)
        
        return modelMP
    end
    ISP_model_s=setup_ISP_optimization(s)
    if !isnothing(previous_solution_for_this_problem)
        allvars = all_variables(ISP_model_s)
        set_start_value.(allvars, previous_solution_for_this_problem)
    end
    set_optimizer_attribute(ISP_model_s, "MIPGap", ISP_optgaps[num_of_solutions_for_ISP[current_ISP_index]+1])
    set_optimizer_attribute(ISP_model_s, "TimeLimit", current_time_limits[current_ISP_index][s])
    #set_optimizer_attribute(ISP_model_s, "MemLimit", 1)
    set_optimizer_attribute(ISP_model_s, "Threads",8)
    set_optimizer_attribute(ISP_model_s, "NumericFocus", 3)
    set_optimizer_attribute(ISP_model_s, "Method", 1)
    set_optimizer_attribute(ISP_model_s, "ScaleFlag", 2)
    #set_optimizer_attribute(ISP_model_s, "Presolve", 0)
    ISP_s_time=@timed begin
        optimize!(ISP_model_s)
    end
    flush(stdout)
    currently_used_time_limit=current_time_limits[current_ISP_index][s]
    while termination_status(ISP_model_s)==TIME_LIMIT && !has_values(ISP_model_s)
        flush(stdout)
        currently_used_time_limit=currently_used_time_limit*2
        
        set_optimizer_attribute(ISP_model_s, "TimeLimit", currently_used_time_limit)
        #@show "limit doubled"
        optimize!(ISP_model_s)
    end
    flush(stdout)
    @spawnat 1 time_limits[current_ISP_index][s]=currently_used_time_limit
    allvars = all_variables(ISP_model_s)
    allvars_solution = value.(allvars)
    @spawnat 1 ISP_solutions[s]=allvars_solution
    flush(stdout)
    #return objval,bestbound,termination_status
    #print(InteractiveUtils.varinfo())

    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    return [objective_value(ISP_model_s),objective_bound(ISP_model_s),termination_status(ISP_model_s),ISP_s_time]
end




MOI.set(RMP, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
#MOI.set(RMP, MOI.RawOptimizerAttribute("Presolve"), 0)
MOI.set(RMP, MOI.RawOptimizerAttribute("Heuristics"), 0)
MOI.set(RMP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1e-9)
MOI.set(RMP, MOI.RawOptimizerAttribute("IntFeasTol"), 1e-9)
MOI.set(RMP, MOI.RawOptimizerAttribute("MIPGap"), RMP_MIPgap)
MOI.set(RMP, MOI.RawOptimizerAttribute("DisplayInterval"), 1)
#MOI.set(RMP, MOI.RawOptimizerAttribute("MIPFocus"), 3)
MOI.set(RMP, MOI.RawOptimizerAttribute("TimeLimit"), RMP_time_limit)
set_optimizer_attribute(RMP, "NumericFocus", 3)
set_optimizer_attribute(RMP, "Method", 1)
set_optimizer_attribute(RMP, "ScaleFlag", 2)
#MOI.set(RMP, MOI.RawOptimizerAttribute("NumericFocus"), 3)
MOI.set(RMP, Gurobi.CallbackFunction(), my_callback_function)
start_time = time()
optimize!(RMP)




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
        sheet["E2"] = objective_value(RMP)
        catch
        sheet["E2"] = "no_sol"
        end
    sheet["F1"] = "objective_bound"
    sheet["F2"] = objective_bound(RMP)
    sheet["G1"] = "MIPGap"
    try
        sheet["G2"] = abs(objective_bound(RMP)-objective_value(RMP))/abs(objective_value(RMP))
        catch
        sheet["G2"] = "no_sol"
        end
    sheet["H1"] = "termination_status"
    sheet["H2"] = string(termination_status(RMP))
    sheet["J1"] = "time"
    sheet["J2"] = solve_time(RMP)
    sheet["K1"] = "algorithm"
    sheet["K2"] = "modified_L_shaped"
    
 
end




#@show dual_setup_t
@show No_good_cuts_added
@show Normal_cuts_added
#@show DSP_times
#@show ISP_times

#total_time_in_DSP = sum([DSP_times[i].time for i in 1:length(DSP_times)])
#total_time_in_int_SP = sum([ISP_times[i].time for i in 1:length(ISP_times)])
#sort([ISP_s_times[i][j].time for i in 1:length(ISP_times) for j in 4])

#plotting results
#= currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
mkdir(currenttime)
current_dir=pwd()
cd(current_dir*"/"*currenttime)
include(current_dir * "/draw_plots_design_no_shedding_fancy_colors.jl")
cd(current_dir) =#
