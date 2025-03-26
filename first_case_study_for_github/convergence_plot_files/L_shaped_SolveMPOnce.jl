

using Distributed
rmprocs(workers())
addprocs(32)
@show workers()
@everywhere using JuMP, Gurobi
using XLSX,CSV,DataFrames
using Random
using Dates
using JLD








##########JuMP specific code ###############




# Declare variables


# see https://github.com/jump-dev/JuMP.jl/blob/master/src/lp_sensitivity2.jl#L302-L412
# Converts an LP model into a standard form. This is used to find F_vec, G_vec, and d_vec automatically for the given problem
# While it might be possible to do this for any LP, it currently only works if all constraints are <= constraints

#calculate hi, Fi, Gi, and di for all i (all scenarios/subproblems)
G_vec = Vector{Any}(undef, length(S))
F_vec = Vector{Any}(undef, length(S))
h_vec = Vector{Any}(undef, length(S))
d_vec = Vector{Any}(undef, length(S))
#set up model for subproblem of scenario i (same Fi,Gi,di as subproblem i)
@everywhere const env = Gurobi.Env(;output_flag=0)
@everywhere Gurobi.GRBsetintparam(env, "OutputFlag", 0)#supress output of dual separation problems
function setup_DSP()
    global dual_setup_t,env
    dual_setup_t = @timed begin
        for s in S
            # see https://github.com/jump-dev/JuMP.jl/blob/master/src/lp_sensitivity2.jl#L302-L412
            # Converts an LP model into a standard form. This is used to find F_vec, G_vec, and d_vec automatically for the given problem
            # While it might be possible to do this for any LP, it currently only works if all constraints are <= constraints

            #set up model that includes constraints for scenario i (same Fi,Gi,di as subproblem i)
            model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
            #variables

            #first stage 
            @variable(model, y[[0],m=M,k=K[1]]) #binary for small gasifiers and turbines
            #= @constraint(model, y .<= 1)
            @constraint(model, -y.<= 0) =#
            #symmetry breaking constraints
            #@constraint(model, [m=M,k=2:length(K[1])], y[[0],m,k] .<= y[[0],m,k-1])
            X = vcat(y...) #vector of all first stage variables
            #second stage
            #z is the auxiliary investment variable
            @variable(model, z[J,M,[s]])
            @constraint(model,[j=[0],m=M,s=[s]], z[j,m,s]<=sum(y[j,m,k] for k in K[1]))
            @constraint(model,[j=[0],m=M,s=[s]], -z[j,m,s]<=-(sum(y[j,m,k] for k in K[1])))
            @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], z[j,m,s]<=0)
            @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], -z[j,m,s]<=0)
            @variable(model, x[I,J,T,[s]])
            @constraint(model, x .<= 1)
            @constraint(model, -x .<= 0)
            @variable(model, q[I,T,[s]]) #fraction of unserved demand for site i
            @constraint(model, q .<= 1)
            @constraint(model, -q .<= 0)
            @constraint(model, [i=I,t=T,s=[s]], sum(x[i,j,t,s] for j=J) +q[i,t,s] <= 1)
            @constraint(model, [i=I,t=T,s=[s]], -sum(x[i,j,t,s] for j=J) -q[i,t,s] <= -1)
            @variable(model, v[J,M,T,[s]])
            @constraint(model, -v[J,M,T,[s]] .<= 0)
            @variable(model, w[J,J,M,T,[s]])
            @constraint(model, -w[J,J,M,T,[s]] .<= 0)
            @constraint(model, [j=J,m=M,t=T,s=[s]], v[j,m,t,s] <= z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t))
            @constraint(model, [j=J,m=M,t=T,s=[s]], -v[j,m,t,s] <= -(z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t)))
            @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= sum(v[j,m,t,s]*u[m] for m in M)) #u[m] is the production rate of module type m
            @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= umax)


            @objective(model, Max, -1*(sum(p[s]*(sum(sum(c[i,j,t]*d[i,t,s]*x[i,j,t,s] for i in I) + sum(h[j,jprime,m,t]*w[j,jprime,m,t,s] for jprime in J, m in M) for j in J, t in T) + sum(penalty_cost[i,t]*q[i,t,s]*d[i,t,s] for i in I, t in T)) for s in [s])))
                        
            Yi = vcat(z.data[:],x.data[:],q.data[:],v.data[:],w.data[:])
            #second stage
           
            #use built in JuMP function to standardize form of LP defined above
            # see https://github.com/jump-dev/JuMP.jl/blob/master/src/lp_sensitivity2.jl#L302-L412
            col, low, up, A_minus_I, bo, co = JuMP._standard_form_matrix(model)
            #Extract F_vec, G_vec, and d_vec from standardized form of LP
            A = A_minus_I[:, 1:end-size(A_minus_I, 1)]
            rhs = up[1+num_variables(model):end]
            F_vec[s] = A[:, 1:length(X)]
            G_vec[s] = A[:, 1+length(X):end]
            d_vec[s] = rhs

            #now calculate hi from objective function of LP
            func = objective_function(model)
            h_vec[s] = zeros(length(Yi))
            for j in 1:length(Yi)
                try
                    h_vec[s][j] = func.terms[Yi[j]]
                catch
                    h_vec[s][j] = 0
                end
            end
        end
    end
    return F_vec,G_vec,d_vec,h_vec
end
F_vec,G_vec,d_vec,h_vec=setup_DSP()
@show dual_setup_t




UB_model=Vector{Any}(undef,num_scenarios)

#Calculate U (upper bound on subproblem cost)
function getUB()
    U = Vector{Any}(undef, length(S))

    U=let problem_data= (num_demand_sites, num_processing_sites, num_scenarios, num_time_periods, I, J, S, M,  T, K,    u, distances, points, production_centers, distances_to_production_centers, distances_between_production_centers, c, d,  h, penalty_cost, g, p, vmax, umax)
        wp = CachingPool(workers())
    pmap(s->do_UB_optimization(s,problem_data),wp,S; retry_delays= ExponentialBackOff(n = 3))
    #continue here, make this a distributed function that solves the model. SEE SHARED ARRAYS LOOKS PROMISING. share only the model variables, solve them in workers, then get results undistributed
    end
    return U
end

@everywhere function do_UB_optimization(s,problem_data)
    (num_demand_sites, num_processing_sites, num_scenarios, num_time_periods, I, J, S, M,  T, K,    u, distances, points, production_centers, distances_to_production_centers, distances_between_production_centers, c, d,  h, penalty_cost, g, p, vmax, umax)=problem_data
    #print(UB_model_s)

    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    function setup_UB_optimization(s)
        model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
       #first stage 
        @variable(model, y[[0],m=M,k=K[1]], Bin) #binary for small gasifiers and turbines
        @constraint(model, y .<= 1)
        @constraint(model, -y.<= 0)
        #symmetry breaking constraints
        @constraint(model, [m=M,k=2:length(K[1])], y[[0],m,k] .<= y[[0],m,k-1])
        #second stage
        #z is the auxiliary investment variable
        @variable(model, z[J,M,[s]], Int)
        @constraint(model,[j=[0],m=M,s=[s]], z[j,m,s]<=sum(y[j,m,k] for k in K[1]))
        @constraint(model,[j=[0],m=M,s=[s]], -z[j,m,s]<=-(sum(y[j,m,k] for k in K[1])))
        @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], z[j,m,s]<=0)
        @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], -z[j,m,s]<=0)
        @variable(model, x[I,J,T,[s]])
        @constraint(model, x .<= 1)
        @constraint(model, -x .<= 0)
        @variable(model, q[I,T,[s]]) #fraction of unserved demand for site i
        @constraint(model, q .<= 1)
        @constraint(model, -q .<= 0)
        @constraint(model, [i=I,t=T,s=[s]], sum(x[i,j,t,s] for j=J) +q[i,t,s] <= 1)
        @constraint(model, [i=I,t=T,s=[s]], -sum(x[i,j,t,s] for j=J) -q[i,t,s] <= -1)
        @variable(model, v[J,M,T,[s]], Int)
        @constraint(model, -v[J,M,T,[s]] .<= 0)
        @variable(model, w[J,J,M,T,[s]], Int)
        @constraint(model, -w[J,J,M,T,[s]] .<= 0)
        @constraint(model, [j=J,m=M,t=T,s=[s]], v[j,m,t,s] <= z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t))
        @constraint(model, [j=J,m=M,t=T,s=[s]], -v[j,m,t,s] <= -(z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t)))
        @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= sum(v[j,m,t,s]*u[m] for m in M)) #u[m] is the production rate of module type m
        @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= umax)


        @objective(model, Max, -1*(sum(p[s]*(sum(sum(c[i,j,t]*d[i,t,s]*x[i,j,t,s] for i in I) + sum(h[j,jprime,m,t]*w[j,jprime,m,t,s] for jprime in J, m in M) for j in J, t in T) + sum(penalty_cost[i,t]*q[i,t,s]*d[i,t,s] for i in I, t in T)) for s in [s])))

        return model
    end


    global env
    UB_model_s=setup_UB_optimization(s)
    set_optimizer_attribute(UB_model_s, "MIPGap", 0.01)
    set_optimizer_attribute(UB_model_s, "TimeLimit", 360)
    set_optimizer_attribute(UB_model_s, "Threads", 1)
    optimize!(UB_model_s)

    #Memory related code
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    return objective_bound(UB_model_s)
end
@time U=getUB()












RMP = direct_model( Gurobi.Optimizer())
#first stage 
@variable(RMP, y[[0],m=M,k=K[1]], Bin) #binary for small gasifiers and turbines
@constraint(RMP, y .<= 1)
@constraint(RMP, -y.<= 0)
#symmetry breaking constraints
@constraint(RMP, [m=M,k=2:length(K[1])], y[0,m,k] <= y[0,m,k-1])
X=vcat(y...)


@variable(RMP, η[s=S] .<= U[s]) #upper bound placed on ηi so problem is always bounded
@objective(RMP, Max, sum(η[i] for i in S) - sum(g[m]*y[0,m,k] for m in M for k in K[1]))
RMP_MIPgap=0.005
RMP_time_limit=3600*3










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
DSP_objval = Vector{Any}(undef, length(S))
bestbound = Vector{Any}(undef, length(S))
incumbent_objval=-Inf
just_improved_incumbent=false

ISP_optgaps=[0.00]
time_lim_initial=3600*3
global JUSTADDEDCUTS=false
global boundsvec=Vector{Any}(undef, 0)
global boundstimes=Vector{Any}(undef, 0)
global cuttypes=Vector{Any}(undef, 0)
global incumbenttimes=Vector{Any}(undef, 0)
global incumbentobjvals=Vector{Any}(undef, 0)
function my_callback_function(cb_data, cb_where::Cint)
    # You can reference variables outside the function as normal
    global DSP_times, ISP_times, Normal_cuts_added, No_good_cuts_added, X, η, incumbent_objval,num_of_solutions_for_ISP, ISP_solutions, time_limits, ISP_s_times, objval, DSP_objval, just_improved_incumbent,time_lim_initial,JUSTADDEDCUTS,boundsvec,boundstimes,cuttypes,incumbenttimes,incumbentobjvals
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
    currenttime= Ref{Cdouble}()
    GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME,currenttime)
    if currenttime[] > RMP_time_limit
        GRBterminate(backend(RMP))
    end


    cuts_added = 0
    NG_cuts_added = 0

    # Get relaxed solution at current node
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    X_val = callback_value.(cb_data, X)
    η_val = callback_value.(cb_data, η)
    mp_sol = vcat(X_val, η_val)
    #check if the submitted mp solution has been found to be feasible  
    if any(all(j .- 1e-9 .<= mp_sol .<= j .+ 1e-9) for j in feasible_mp_sol_list) #  
        return
    end

    ##Memory related code
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    #record bound if needed
    if JUSTADDEDCUTS==true
        currentbound=Ref{Cdouble}()
        if cb_where == GRB_CB_MIPSOL
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBND,currentbound)
        elseif cb_where == GRB_CB_MIPNODE
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND,currentbound)
        end
        push!(boundsvec,currentbound[])
        JUSTADDEDCUTS=false
    end
    if any(all(j .- 1e-9 .<= X_val .<= j .+ 1e-9) for j in X_list_dsp_sol) #if the dsp has been solved for this X before
        index_of_old_dsp_solve=findfirst(all(j .- 1e-9 .<= X_val .<= j .+ 1e-9) for j in X_list_dsp_sol)
        if all(η_val .<= objval_list_dsp_sol[index_of_old_dsp_solve] .+ 1e-6)
            need_to_do_dsp=false
            print("skipped doing dsp ")
        else
            need_to_do_dsp=true
            print("did dsp for same X again ")
        end
    else
        need_to_do_dsp=true
    end

    if need_to_do_dsp 
        #solve dual separation problems
        DSP_time = @timed begin
            results=let (F_vec,G_vec,d_vec,h_vec,X_val)=(F_vec,G_vec,d_vec,h_vec,X_val)
                wp = CachingPool(workers())
                #get primal_status, v, objective value, and ustar in pmap'd_vec function
                pmap(s->do_DSP_optimization(s,(F_vec[s],G_vec[s],d_vec[s],h_vec[s],X_val),savedVBases[s],savedCBases[s]),wp,S; retry_delays= ExponentialBackOff(n = 3)) #continue here, make this a distributed function that solves the model. SEE SHARED ARRAYS LOOKS PROMISING. share only the model variables, solve them in workers, then get results undistributed
            end

            for s in S
                if results[s][1] == INFEASIBILITY_CERTIFICATE #the primal solution is an unbounded ray
                    print("INFEASIBLE DSP \n")
                    ν=  results[s][2]
                    con=@build_constraint(ν' * (d_vec[s] - F_vec[s] * X) >= 0)
                    MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                    cuts_added = 1 + cuts_added
                elseif results[s][1] == FEASIBLE_POINT
                    DSP_objval[s] = results[s][3]
                    η_val[s]
                    if η_val[s] > results[s][3] + 1e-6
                        ustar =  results[s][2]
                        con= @build_constraint(ustar' * (d_vec[s] - F_vec[s] * X) >= η[s])
                        MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                        cuts_added = 1 + cuts_added
                    end
                else
                    print(results[s][1])
                end
            end
        end
        push!(X_list_dsp_sol, X_val)
        push!(objval_list_dsp_sol,DSP_objval)
        DSP_times = push!(DSP_times, DSP_time)
        Normal_cuts_added = push!(Normal_cuts_added, cuts_added)
        if cuts_added >= 1
            JUSTADDEDCUTS=true
            push!(boundstimes,time()-start_time)
            push!(cuttypes,"Benders")
            return
        end
    end
    cuttype_code="none"
    #if no cut was added and the solution is integer
    if all([(- 1e-9 <= j <= 1e-9)||(1-1e-9 <= j <= 1+1e-9) for j in X_val])
        current_ISP_index=findfirst(all(j .- 1e-9 .<= X_val .<= j .+ 1e-9) for j in X_list_int_sol)
        if isnothing(current_ISP_index)
            push!(X_list_int_sol, X_val)
            current_ISP_index=length(X_list_int_sol)
            push!(num_of_solutions_for_ISP,0)
            push!(time_limits,ones(length(S))*time_lim_initial)
        end
        #solve integer version of subproblems to generate no-good cuts
        while num_of_solutions_for_ISP[current_ISP_index] < length(ISP_optgaps)
            integer_SP_time = @timed begin
                results=let problem_data= (num_demand_sites, num_processing_sites, num_scenarios, num_time_periods, I, J, S, M,  T, K,    u, distances, points, production_centers, distances_to_production_centers, distances_between_production_centers, c, d,  h, penalty_cost, g, p, vmax, umax)
                    wp = CachingPool(workers())
                    pmap(s->do_ISP_optimization(s,problem_data,ISP_optgaps,num_of_solutions_for_ISP,current_ISP_index,time_limits,ISP_solutions[s],X_val),wp,S; retry_delays= ExponentialBackOff(n = 3)) 
                end
                adjusted_num_of_solutions_for_ISP=num_of_solutions_for_ISP[current_ISP_index]
                original_num_of_solutions_for_ISP=num_of_solutions_for_ISP[current_ISP_index]
                if ISP_optgaps[num_of_solutions_for_ISP[current_ISP_index]+1]==0
                    cuttype_code="full_opt"
                else
                    cuttype_code="subopt"
                end
                for s in S
                    objval[s] = results[s][1]
                    bestbound[s]=results[s][2]
                    if ISP_optgaps[num_of_solutions_for_ISP[current_ISP_index]+1]==0 && results[s][3]==TIME_LIMIT
                        adjusted_num_of_solutions_for_ISP=length(ISP_optgaps)-2 #this will cause the current set of subproblems to be run again with a gap of 0
                        time_limits[current_ISP_index][s]=time_limits[current_ISP_index][s]*5
                        print("need to run with gap of 0 again")
                        cuttype_code="some_full_opt_some_time_limit"
                    end
                    #print(termination_status(ISP_models[s]))
                    #print(relative_gap(ISP_models[s]))
                    if results[s][3] == INFEASIBLE
                        #build and submit no-good feasibility cut
                        print("\n infeasible ISP \n")
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
                            MOI.submit(RMP, MOI.LazyConstraint(cb_data), con)
                            NG_cuts_added = NG_cuts_added + 1
                            print("added ng cut")
                        end
                    else
                        print(results[s][3])
                        error()
                    end
                end
                num_of_solutions_for_ISP[current_ISP_index]=adjusted_num_of_solutions_for_ISP
                original_num_of_solutions_for_ISP#this lets the @timed macro above record the solution number for the isp solve
            end
            push!(ISP_times, integer_SP_time)
            push!(ISP_s_times,[results[s][4] for s in S])
            num_of_solutions_for_ISP[current_ISP_index] +=1
            if NG_cuts_added >=1 #if a no good cut was added to remove the current mp solution
                #check to see if incumbent can be improved using feasible solution found in callback
                #= incumbent_objval= Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST,incumbent_objval) =# #CURRENTLY BUGGED DUE TO GUROBI BUG, SEE https://support.gurobi.com/hc/en-us/community/posts/13366202203281-Wrong-values-of-best-incumbent-objective-value
               
                y_val = callback_value.(cb_data, y)
                ISP_objval= (- (sum(g[m]*y_val[0,m,k] for m in M for k in K[1])))+sum(objval[i] for i in S)
                #@show ISP_objval
                #@show incumbent_objval
                if ISP_objval > incumbent_objval
                    MOI.submit(RMP, MOI.HeuristicSolution(cb_data), vcat(X,η), convert.(AbstractFloat, vcat(X_val,objval)))
                    print("new incumbent found \n")
                    just_improved_incumbent=true
                    incumbent_objval=ISP_objval
                    push!(feasible_mp_sol_list,vcat(X_val,objval))
                    push!(incumbenttimes,time()-start_time)
                    push!(incumbentobjvals,ISP_objval)
                end
                No_good_cuts_added = push!(No_good_cuts_added, NG_cuts_added)
                JUSTADDEDCUTS=true
                push!(boundstimes,time()-start_time)
                push!(cuttypes,cuttype_code)
                return
            end
            
        end
        
        print("got to end of isp part (uh oh?)")
        
        y_val = callback_value.(cb_data, y)
        ISP_objval= (- (sum(g[m]*y_val[0,m,k] for m in M for k in K[1])))+sum(objval[i] for i in S)
        if ISP_objval > incumbent_objval
            print("new incumbent found \n")
            just_improved_incumbent=true
            incumbent_objval=ISP_objval
            push!(feasible_mp_sol_list,vcat(X_val,objval))
        end
        No_good_cuts_added = push!(No_good_cuts_added, NG_cuts_added)
    end

    return
end

@everywhere function do_DSP_optimization(s,problem_data,savedVBasis,savedCBasis)
    (F_s,G_s,d_s,h_s,X_val)=problem_data
    #print(UB_model_s)
    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end
    
    function setup_DSP_optimization(s)
        DSP_model_s = direct_model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0, "InfUnbdInfo" => 1))
        @variable(DSP_model_s, u[1:length(d_s)] >= 0)
        #The G_vec, F_vec, h_vec, and d_vec were obtained from form of constraints with free decision variables
        #i.e., the subproblems have free decision variables 
        @constraint(DSP_model_s, G_s' * u .== h_s)
        @objective(DSP_model_s, Min, (d_s - F_s * X_val)' * u)
        return DSP_model_s
    end
    DSP_model_s=setup_DSP_optimization(s)
    
    
    #set_optimizer_attribute(DSP_model_s, "MemLimit", 1)
    set_optimizer_attribute(DSP_model_s, "Threads", 1)
    if !isnothing(savedVBasis) && !isnothing(savedCBasis)
        set_optimizer_attribute(DSP_model_s, "LPWarmStart", 2)
        allvars=all_variables(DSP_model_s)
        allcons=all_constraints(DSP_model_s;include_variable_in_set_constraints=false)
        MOI.set.(DSP_model_s,Gurobi.VariableAttribute("VBasis"),allvars,savedVBasis)
        MOI.set.(DSP_model_s,Gurobi.ConstraintAttribute("CBasis"),allcons,savedCBasis)
    end
    optimize!(DSP_model_s)
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
end

@everywhere function do_ISP_optimization(s,problem_data,ISP_optgaps,num_of_solutions_for_ISP,current_ISP_index,time_limits,previous_solution_for_this_problem,X_val)
    (num_demand_sites, num_processing_sites, num_scenarios, num_time_periods, I, J, S, M,  T, K,    u, distances, points, production_centers, distances_to_production_centers, distances_between_production_centers, c, d,  h, penalty_cost, g, p, vmax, umax)=problem_data

    #memory related
    GC.gc(true)
    GC.safepoint()
    try
    ccall(:malloc_trim, Cvoid, (Cint,), 0)
    catch
    end

    function setup_ISP_optimization(s)
        model = direct_model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "LogToConsole" => 0))
        #first stage 
        @variable(model, y[[0],m=M,k=K[1]], Bin) #binary for small gasifiers and turbines
        @constraint(model, y .<= 1)
        @constraint(model, -y.<= 0)
        #symmetry breaking constraints
        @constraint(model, [m=M,k=2:length(K[1])], y[[0],m,k] .<= y[[0],m,k-1])
        ispX = vcat(y...) #vector of all first stage variables
        @constraint(model, Xcon, ispX .== X_val)# This fixes the first stage variables to the RMP's solution
        #second stage
        #z is the auxiliary investment variable
        @variable(model, z[J,M,[s]], Int)
        @constraint(model,[j=[0],m=M,s=[s]], z[j,m,s]<=sum(y[j,m,k] for k in K[1]))
        @constraint(model,[j=[0],m=M,s=[s]], -z[j,m,s]<=-(sum(y[j,m,k] for k in K[1])))
        @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], z[j,m,s]<=0)
        @constraint(model,[j=1:num_processing_sites,m=M,s=[s]], -z[j,m,s]<=0)
        @variable(model, x[I,J,T,[s]])
        @constraint(model, x .<= 1)
        @constraint(model, -x .<= 0)
        @variable(model, q[I,T,[s]]) #fraction of unserved demand for site i
        @constraint(model, q .<= 1)
        @constraint(model, -q .<= 0)
        @constraint(model, [i=I,t=T,s=[s]], sum(x[i,j,t,s] for j=J) +q[i,t,s] <= 1)
        @constraint(model, [i=I,t=T,s=[s]], -sum(x[i,j,t,s] for j=J) -q[i,t,s] <= -1)
        @variable(model, v[J,M,T,[s]], Int)
        @constraint(model, -v[J,M,T,[s]] .<= 0)
        @variable(model, w[J,J,M,T,[s]], Int)
        @constraint(model, -w[J,J,M,T,[s]] .<= 0)
        @constraint(model, [j=J,m=M,t=T,s=[s]], v[j,m,t,s] <= z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t))
        @constraint(model, [j=J,m=M,t=T,s=[s]], -v[j,m,t,s] <= -(z[j,m,s]+sum(-sum(w[j,jprime,m,tau,s] for jprime in J)+sum(w[jprime,j,m,tau,s] for jprime in J) for tau=1:t)))
        @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= sum(v[j,m,t,s]*u[m] for m in M)) #u[m] is the production rate of module type m
        @constraint(model, [j=J,t=T,s=[s]], sum(d[i,t,s]*x[i,j,t,s] for i=I) <= umax)


        @objective(model, Max, -1*(sum(p[s]*(sum(sum(c[i,j,t]*d[i,t,s]*x[i,j,t,s] for i in I) + sum(h[j,jprime,m,t]*w[j,jprime,m,t,s] for jprime in J, m in M) for j in J, t in T) + sum(penalty_cost[i,t]*q[i,t,s]*d[i,t,s] for i in I, t in T)) for s in [s])))

        return model
    end
    ISP_model_s=setup_ISP_optimization(s)
    if !isnothing(previous_solution_for_this_problem)
        allvars = all_variables(ISP_model_s)
        set_start_value.(allvars, previous_solution_for_this_problem)
    end
    set_optimizer_attribute(ISP_model_s, "MIPGap", ISP_optgaps[num_of_solutions_for_ISP[current_ISP_index]+1])
    set_optimizer_attribute(ISP_model_s, "TimeLimit", time_limits[current_ISP_index][s])
    #set_optimizer_attribute(ISP_model_s, "MemLimit", 1)
    set_optimizer_attribute(ISP_model_s, "Threads",1)
    ISP_s_time=@timed begin
        optimize!(ISP_model_s)
    end
    while termination_status(ISP_model_s)==TIME_LIMIT && !has_values(ISP_model_s)
        time_limits[current_ISP_index][s]=time_limits[current_ISP_index][s]*2
        set_optimizer_attribute(ISP_model_s, "TimeLimit", time_limits[current_ISP_index][s])
        @show "limit doubled"
        optimize!(ISP_model_s)
    end
    newtimelimit=time_limits[current_ISP_index][s]
    @spawnat 1 time_limits[current_ISP_index][s]=newtimelimit
    allvars = all_variables(ISP_model_s)
    allvars_solution = value.(allvars)
    @spawnat 1 ISP_solutions[s]=allvars_solution
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
MOI.set(RMP, MOI.RawOptimizerAttribute("MIPFocus"), 0)
MOI.set(RMP, MOI.RawOptimizerAttribute("TimeLimit"), RMP_time_limit)
MOI.set(RMP, Gurobi.CallbackFunction(), my_callback_function)
global start_time = time()
optimize!(RMP)
push!(boundsvec,objective_bound(RMP))



currenttime = Dates.format(Dates.now(), "dd_u_yyyy_HH-MM-SS")
dirstring="$(num_demand_sites)_$(num_processing_sites)_$(num_scenarios)_$(num_time_periods)_"*currenttime
using Base.Filesystem
# Get the full path of the current file
full_path = @__FILE__
# Get just the file name
#= file_name = basename(full_path)
open(file_name[1:end-3]*"_solution.txt", "w") do io
    for x in all_variables(RMP)
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
    sheet["E2"] = objective_value(RMP)
    sheet["F1"] = "objective_bound"
    sheet["F2"] = objective_bound(RMP)
    sheet["G1"] = "MIPGap"
    sheet["G2"] = abs(objective_bound(RMP)-objective_value(RMP))/abs(objective_value(RMP))
    sheet["H1"] = "termination_status"
    sheet["H2"] = string(termination_status(RMP))
    sheet["J1"] = "time"
    sheet["J2"] = solve_time(RMP)
    sheet["K1"] = "algorithm"
    sheet["K2"] = "L_shaped"
 
end
filestring=dirstring*"results.jld"
@save filestring boundstimes boundsvec cuttypes incumbenttimes incumbentobjvals 


@show value.(y)
@show dual_setup_t
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
