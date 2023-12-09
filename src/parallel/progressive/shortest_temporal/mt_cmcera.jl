function upper_bound_average_diameter(delta::Float64,diam::Int64,tdd::Array{Int64},sample_size::Int64,verbose::Bool=false,norm::Float64 = 1.0)::Float64
    avg_dist::Float64 = 0.0
    # Diam is the vertex diameter it's already diam +1 
    for i in 1:diam
        avg_dist += (i-1) * tdd[i]
    end
    avg_dist = avg_dist/sample_size * norm
    # Upper bound using Bernstein bound
    log_term_avg_dist::Float64 = log(1. /delta)
    c_term_avg_dist::Float64 = (diam - 2)*log_term_avg_dist/sample_size
    average_diam_ub_b = avg_dist + c_term_avg_dist + sqrt(2*c_term_avg_dist*diam + c_term_avg_dist^2)
    var_estimate_diam::Float64 = 0.0
    # Upper bound using Empirical Bernstein bound
    for i in 1:diam
        var_estimate_diam+= (tdd[i] - avg_dist)^2
    end
    var_estimate_diam = var_estimate_diam/(sample_size-1) * norm
    #=
    for i in 1:diam
        for j in (i+1):diam
            var_estimate_diam += ((i-1)-(j-1))^2 * tdd[i]/sample_size * tdd[j]/(sample_size-1) * norm
        end
    end
    =#
    log_term_avg_dist = log(2/delta)
    average_diam_ub_eb::Float64 = avg_dist + 7/3 * (diam -2) * log_term_avg_dist/sample_size + sqrt(2*var_estimate_diam*log_term_avg_dist / sample_size)
    avg_diam_upperbound = min(average_diam_ub_b,average_diam_ub_eb)
    if verbose
        println("Average diameter "*string(avg_dist))
        println("Average diameter UB (Bernstein) "*string(average_diam_ub_b))
        println("Average diameter UB (Emp-Bernstein) "*string(average_diam_ub_eb))
        println("Variance estimate average diameter "*string(var_estimate_diam))
        flush(stdout)
    end
    if diam -2 >0
        return min(avg_diam_upperbound,diam-2)
    end
    return avg_diam_upperbound 
end
function upper_bound_top_1_tbc(top1_est_bc::Float64,delta::Float64,sample_size::Int64)::Float64
    log_term_top_1_tbc::Float64 = log(1/delta)
    const_term_top_1_bc::Float64 =  log_term_top_1_tbc/sample_size
    top1bc_upperbound::Float64  = top1_est_bc +  const_term_top_1_bc + sqrt(2*const_term_top_1_bc*top1_est_bc+const_term_top_1_bc^2)
    return min(1.0,top1bc_upperbound)
end
function number_samples_bound_function(x::Float64,rho::Float64,delta::Float64,eps::Float64)::Float64
    v_x::Float64 = x*(1-x)
    arg_h::Float64 = eps/v_x
    denom::Float64 = v_x * ((1+arg_h) * log(1+arg_h)-arg_h)
    return log(2*rho / (x*delta))/denom

end
function upper_bound_samples(max_tbc::Float64,max_var::Float64, avg_dist::Float64,eps::Float64, delta_bound::Float64,debug::Bool = false)::Float64
    x_hat::Float64 = 0.0
    x_hat_l::Float64 = 0.5-sqrt(eps/3)
    x_hat_l = max(x_hat_l,0.0)
    x_hat_h::Float64 = 0.5 
    v_x::Float64 = 0.0
    arg_h::Float64 = 0.0
    f_val::Float64 = 0.0
    while x_hat_h - x_hat_l > 0.0001
        x_hat = (x_hat_h + x_hat_l) /2.0
        if debug
            println(" X_hat "*string(x_hat))
        end
        v_x = x_hat * (1-x_hat)
        arg_h = eps/v_x
        f_val = v_x * ((1+arg_h) * log(1+arg_h)-arg_h)
        if f_val <= 2*eps^2
            x_hat_h = x_hat
        else
            x_hat_l = x_hat
        end
    end
    x_hat = x_hat_h 
    x_h::Float64 = min(x_hat,max_tbc)
    x_h_var::Float64 = 1.0
    if max_var < 0.25
        x_h_var = 0.5 - sqrt(0.25-max_var)
        x_h = min(x_h,x_h_var)
    end
    x_l::Float64 = x_h
    step::Float64 = x_h /1000.
    num_samples_bound_high::Float64 = number_samples_bound_function(x_h,avg_dist,delta_bound,eps) 
    num_samples_bound::Float64 =num_samples_bound_high+1
    while num_samples_bound > num_samples_bound_high
        x_l = x_h - step
        if x_l > 0
            num_samples_bound = number_samples_bound_function(x_l,avg_dist,delta_bound,eps)
            if num_samples_bound > num_samples_bound_high
                x_h = x_l
                num_samples_bound_high = num_samples_bound
            end
        else
            num_samples_bound = num_samples_bound_high -1
        end
    end
    return num_samples_bound_high
end

function _check_stopping_condition_topk!(approx_top_k::Array{Tuple{Int64,Float64}},betweenness::Array{Float64},wv::Array{Float64},last_stopping_samples::Float64,num_samples::Int64,eps::Float64,delta::Float64,iteration::Int64,second_phase::Bool,diam::Int64,tdd::Array{Int64},sample_size::Int64,mc_trials::Int64,partition_index::Array{Int64},partitions_ids_map::Dict{Int64,Int64},emp_wimpy_vars::Array{Float64},mcrade::Array{Float64},number_of_non_empty_partitions::Int64,omega::Vector{Float64},norm::Float64,has_to_stop::Vector{Bool},conv_numb::Vector{Int64})
    n::Int64 = lastindex(betweenness)
    num_samples_d::Float64 = num_samples
    delta_for_progressive_bound::Float64 = delta/2^iteration
   
    sup_bcest_partition::Array{Float64} = zeros(n)
    sup_empwvar_partition::Array{Float64} = zeros(n)
    epsilon_partition::Array{Float64} = ones(n)
    max_mcera_partition::Array{Float64} = [-num_samples for i in 1:(mc_trials*n) ]
    # Update MCERA
    for i in 1:n
        v_rade_idx = i*mc_trials
        node_partition_idx = partition_index[i]
        mapped_partition_index = partitions_ids_map[node_partition_idx]
        sup_bcest_partition[mapped_partition_index] = max(sup_bcest_partition[mapped_partition_index],betweenness[i] )
        sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index],emp_wimpy_vars[i])
        mcera_partition_index = mc_trials*mapped_partition_index
        for j in 1:mc_trials
            max_mcera_partition[j+mcera_partition_index] = max(max_mcera_partition[j+mcera_partition_index] ,mcrade[v_rade_idx+j])
        end
    end
    mcera_partition_avg::Array{Float64} = zeros(number_of_non_empty_partitions)
    mcera_avg::Float64 = 0.0
    mcera_partition_index::Int64 = 0.0 
    delta_each_partition::Float64 = delta_for_progressive_bound/number_of_non_empty_partitions
    for i in 1:number_of_non_empty_partitions
        mcera_avg = 0.0
        mcera_partition_index = mc_trials *i
        for j in 1:mc_trials
            mcera_avg+=max_mcera_partition[mcera_partition_index+j]/mc_trials
        end
        mcera_avg = mcera_avg/num_samples_d
        mcera_partition_avg[i] = mcera_avg
        sup_emp_wimpy_var = sup_empwvar_partition[i]/num_samples_d
        current_eps = epsilon_mcrade(sup_emp_wimpy_var,mcera_avg,delta_each_partition,num_samples_d,mc_trials)
        epsilon_partition[i] = current_eps
    end
    sup_eps::Float64 = 0.0
    converged::Bool = true
    number_of_converged::Int64 = 0
    for i in 1:lastindex(approx_top_k)
        u = approx_top_k[i][1]
        node_partition_idx = partition_index[u]
        map_node_partition_index = partitions_ids_map[node_partition_idx]
        eps_current_node = epsilon_partition[map_node_partition_index]
        converged = converged & (eps_current_node <= eps)
        #println("ξ_cur "*string(u)*" = "*string(eps_current_node))
        if eps_current_node <= eps 
            number_of_converged+=1
        end
    end
   
    has_to_stop[1]= converged
    conv_numb[1] = number_of_converged
    return nothing
end
function _check_stopping_condition!(betweenness::Array{Float64},wv::Array{Float64},last_stopping_samples::Float64,num_samples::Int64,eps::Float64,delta::Float64,iteration::Int64,second_phase::Bool,diam::Int64,tdd::Array{Int64},sample_size::Int64,mc_trials::Int64,partition_index::Array{Int64},partitions_ids_map::Dict{Int64,Int64},emp_wimpy_vars::Array{Float64},mcrade::Array{Float64},number_of_non_empty_partitions::Int64,omega::Vector{Float64},norm::Float64,has_to_stop::Vector{Bool})
    n::Int64 = lastindex(betweenness)
    num_samples_d::Float64 = num_samples
    delta_for_progressive_bound::Float64 = delta/2^iteration
    #println("Checking stopping condition at iteration "*string(iteration)*" sample size "*string(num_samples)*" δ = "*string(delta_for_progressive_bound))
    #=
    if second_phase
        avg_diam_upperbound::Float64 = upper_bound_average_diameter(delta_for_progressive_bound,diam,tdd,sample_size,false,norm)
        top1_est_bc::Float64 = maximum(betweenness)/num_samples_d
        top1bc_upperbound::Float64 = upper_bound_top_1_tbc(top1_est_bc,delta_for_progressive_bound,sample_size)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(maximum(wv)/num_samples_d,delta_for_progressive_bound,sample_size)
        max_num_samples::Int64 = trunc(Int,upper_bound_samples(top1bc_upperbound,wimpy_var_upper_bound,avg_diam_upperbound,eps,delta_for_progressive_bound))
        if last_stopping_samples > max_num_samples
            last_stopping_samples = max_num_samples
            
            omega[1] = last_stopping_samples
            println("New stopping condition update, last stopping samples "*string(last_stopping_samples))
        end
        if max_num_samples <= num_samples
            println("New stopping condition TRUE")
        end
    end
    =#
    sup_bcest_partition::Array{Float64} = zeros(n)
    sup_empwvar_partition::Array{Float64} = zeros(n)
    epsilon_partition::Array{Float64} = ones(n)
    max_mcera_partition::Array{Float64} = [-num_samples for i in 1:(mc_trials*n) ]
    # Update MCERA
    for i in 1:n
        v_rade_idx = i*mc_trials
        node_partition_idx = partition_index[i]
        mapped_partition_index = partitions_ids_map[node_partition_idx]
        sup_bcest_partition[mapped_partition_index] = max(sup_bcest_partition[mapped_partition_index],betweenness[i] )
        sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index],emp_wimpy_vars[i])
        mcera_partition_index = mc_trials*mapped_partition_index
        for j in 1:mc_trials
            max_mcera_partition[j+mcera_partition_index] = max(max_mcera_partition[j+mcera_partition_index] ,mcrade[v_rade_idx+j])
        end
    end
    mcera_partition_avg::Array{Float64} = zeros(number_of_non_empty_partitions)
    mcera_avg::Float64 = 0.0
    mcera_partition_index::Int64 = 0.0 
    delta_each_partition::Float64 = delta_for_progressive_bound/number_of_non_empty_partitions
    for i in 1:number_of_non_empty_partitions
        mcera_avg = 0.0
        mcera_partition_index = mc_trials *i
        for j in 1:mc_trials
            mcera_avg+=max_mcera_partition[mcera_partition_index+j]/mc_trials
        end
        mcera_avg = mcera_avg/num_samples_d
        mcera_partition_avg[i] = mcera_avg
        sup_emp_wimpy_var = sup_empwvar_partition[i]/num_samples_d
        current_eps = epsilon_mcrade(sup_emp_wimpy_var,mcera_avg,delta_each_partition,num_samples_d,mc_trials)
        epsilon_partition[i] = current_eps
    end
    sup_eps::Float64 = 0.0
    for i in 1:number_of_non_empty_partitions
        sup_eps = max(sup_eps,epsilon_partition[i])
    end
    if sup_eps <= eps
        println("MCRADE STOPS with ξ : "*string(sup_eps))
        println("MCRADE STOPS at iteration : "*string(iteration))
        println("MCRADE STOPS at sample size : "*string(num_samples))
        flush(stdout)
    else
     println("MC-RADE ξ "*string(sup_eps)*" target ε  "*string(eps) )
     flush(stdout)
    end
    has_to_stop[1]= (sup_eps <= eps)
    return nothing
end
function epsilon_mcrade(sup_emp_wimpy_var,mcera,delta,num_samples,mc_trials)
    mcera = max(mcera,0.0)
    log_term_mcrade::Float64 = log(5/delta)/num_samples
    var_ub::Float64 = sup_emp_wimpy_var +log_term_mcrade+sqrt(log_term_mcrade^2 +2 * log_term_mcrade*sup_emp_wimpy_var) 
    era_ub::Float64 = mcera + sqrt(4*sup_emp_wimpy_var *log_term_mcrade / mc_trials)
    ra_ub::Float64 = era_ub + log_term_mcrade + sqrt(log_term_mcrade^2 + 2*log_term_mcrade * era_ub)
    eps_ub::Float64 = 2* ra_ub+ sqrt(2*log_term_mcrade*(var_ub+4*ra_ub))
    return eps_ub
end 

function get_next_stopping_sample(ss::Float64,iteration_index::Int64)
    ss = ss * 1.2
    iteration_index +=1
    return ss,iteration_index
end

function threaded_progressive_cmcera(tg::temporal_graph,eps::Float64,delta::Float64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,force_gc::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    ntasks = nthreads()
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    wv::Array{Float64} = Array{Float64}([])
    emp_w_node::Float64 = 0.0
    min_inv_w_node::Float64 = 0.0
    # Partitions
    number_of_non_empty_partitions::Int64 = 0.0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(tg.num_nodes)
    part_idx::Int64 = 1
    # TBC
    tmp_has_to_stop::Array{Bool} = Array{Bool}([false])
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0
    if (diam == -1)
        println("Approximating (sh)-Temporal Diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,0,0.9,false)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        flush(stdout)
    end
    start_time_bootstrap = time()

    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap using Variance")
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    if vc_upper_bund == true
        println("Bootstrap using VC dimension")
    else
        println("Bootstrap using Variance")
    end
    flush(stdout)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "rtb"
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
            if (Sys.free_memory() / Sys.total_memory() < 0.1)
                clean_gc()
                sleep(10)
            end
        end
    end

    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    println("Empirical peeling phase:")
    flush(stdout)
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0

    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i])
        max_wv = max(max_wv,wv[i])
        emp_w_node = wv[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    if !vc_upper_bund
        # Upper bound on the average distance
        avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
        # Upper bound on the top-1 temporal betweenness
        top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
        # define delta_for_progressive_bound
        println("Average (sh)-temporal path (upper bound) "*string(avg_diam_ub))
        # Upper limit on number of samples
        max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
        omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
        println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
        println("Sup tbc estimation "*string(max_tbc))
        println("Sup empirical wimpy variance "*string(max_wv/tau))
        flush(stdout)
    end
    iteration_index::Int64 =1 
    
    local_wv = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for _ in 1:ntasks]
    local_sp_lengths = [zeros(tg.num_nodes) for _ in 1:ntasks]
    
    omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    if max_num_samples > 0
        omega = max_num_samples
    end
    
    max_num_samples = omega
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    finish_bootstrap = string(round(time() - start_time_bootstrap; digits=4))
    println("Bootstrap completed in "*finish_bootstrap)
    println("Inferring initial sample size for the geometric sampler")
    flush(stdout)
    while first_sample_upper - first_sample_lower> 10
        num_samples = (first_sample_upper+first_sample_lower)÷2
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples) + log(2/delta)/num_samples/3
        if eps_guess > eps
            first_sample_lower = num_samples
        else
            first_sample_upper = num_samples
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("Maximum number of iterations "*string(last_stopping_samples))
    println("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stdout)
    next_stopping_samples::Float64 = first_stopping_samples
    prev_stopping_samples::Float64 = 0.0
    has_to_stop::Bool = false
    num_samples = 0
    sample_i::Int64 = 0
    number_of_gc_calls::Int64 = 0
    overall_cleaning_time::Float64 = 0.0
    while !has_to_stop
        sample_i = trunc(Int,next_stopping_samples-prev_stopping_samples)
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
    
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "rtb"
                    _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
                if (Sys.free_memory() / Sys.total_memory() < 0.1)
                    clean_gc()
                    sleep(10)
                end
            end 
        end
        

        num_samples += sample_i

        if num_samples >= omega
            has_to_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(num_samples)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end
    
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            
            betweenness = reduce(+, local_temporal_betweenness)
     
            wv = reduce(+,local_wv)
            #sp_lengths = reduce(+,local_sp_lengths) 
            r_mcrade = reduce(+,mcrade)
                
            tmp_omega = Vector{Float64}([omega])
            tmp_has_to_stop = Vector{Bool}([false])
            _check_stopping_condition!(betweenness,wv,last_stopping_samples,num_samples,eps,delta,iteration_index,true,diam,sp_lengths,num_samples,mc_trials,partition_index,partitions_ids_map,wv,r_mcrade,number_of_non_empty_partitions,tmp_omega,norm,tmp_has_to_stop)
            
            omega = tmp_omega[1]
            has_to_stop = tmp_has_to_stop[1]
            if has_to_stop
                println("Progressive sampler converged!")
                flush(stdout)
            else
                prev_stopping_samples = next_stopping_samples
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
                flush(stdout)
                if (force_gc)
                    clean_gc()
                end
                
            end
                    
        end


    end

    println("(SH)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness.*[1/num_samples],num_samples,max_num_samples,time()-start_time


end


function threaded_progressive_cmcera_relative_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,geo::Float64 = 1.2,empirical_peeling_a::Float64 = 2.0)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    ntasks = nthreads()
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    union_sample::Int64 = min(tg.num_nodes,max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/ntasks),k+20))
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    wv::Array{Float64} = Array{Float64}([])
    emp_w_node::Float64 = 0.0
    min_inv_w_node::Float64 = 0.0
    # Partitions
    number_of_non_empty_partitions::Int64 = 0.0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(tg.num_nodes)
    part_idx::Int64 = 1
    # TBC
    tmp_has_to_stop::Array{Bool} = Array{Bool}([false])
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    approx_top_k =  Array{Tuple{Int64,Float64}}([])
    top_k_result::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])

    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0
    local_distinct_nodes_top_k::Array{Array{Int64}} = [[0] for _ in 1:ntasks]
    if (diam == -1)
        println("Approximating (sh)-Temporal Diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,0,0.9,false)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        flush(stdout)
    end
    start_time_bootstrap = time()
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap using Variance")
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    bootstrap_sampling_completed::Bool = false
    n_pairs_bootstrap::Int64 = 0
    while !bootstrap_sampling_completed
        @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    #_sh_accumulate_trk_rel_topk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sh_accumulate_onbra_rel_topk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t],local_distinct_nodes_top_k[t])
                elseif algo == "rtb"
                    #_sh_accumulate_rtb_rel_topk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
            end
        end
        n_pairs_bootstrap += tau
        bootstrap_sampling_completed = (reduce(+,local_distinct_nodes_top_k)[1] >= 1.5*k)
    end
    tau = n_pairs_bootstrap
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    println("Empirical peeling phase:")
    flush(stdout)
    est_kth_bc::Float64 = 10.0
    for i in 1:tg.num_nodes
        push!(approx_top_k,(i,betweenness[i]))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    est_kth_bc = approx_top_k[k][2] +1/tau
    approx_top_k =  Array{Tuple{Int64,Float64}}([])
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0
    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i])
        max_wv = max(max_wv,wv[i])
        emp_w_node = wv[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    # Upper bound on the average distance
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
    # Upper bound on the top-1 temporal betweenness
    top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
    # define delta_for_progressive_bound
    println("Average (sh)-temporal path (upper bound) "*string(avg_diam_ub))
    # Upper limit on number of samples
    #max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    #omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    
    iteration_index::Int64 =1 
    
    local_wv = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for _ in 1:ntasks]
    local_sp_lengths = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega = 10^15
    if max_num_samples > 0
        omega = max_num_samples
    end
    max_num_samples = omega
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc estimation "*string(max_tbc))
    println("Sup empirical wimpy variance "*string(max_wv/tau))
    flush(stdout)
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    violated_condition_lower_1::Bool = false
    violated_condition_lower_2::Bool = false
    violated_condition_upper_1::Bool = false
    while first_sample_upper - first_sample_lower > 10
        num_samples = (first_sample_upper+first_sample_lower)/2
        eps_guess = sqrt(2*est_kth_bc*log(2/delta/est_kth_bc) /num_samples) + log(2/delta/est_kth_bc)/num_samples/3
        #println("Eps guess "*string(eps_guess))
        violated_condition_lower_1 = (est_kth_bc + eps_guess > est_kth_bc/(1-eps))
        violated_condition_upper_1 = (est_kth_bc - eps_guess < est_kth_bc/(1+eps))
        violated_condition_lower_2 = (eps_guess >= (est_kth_bc *eps*((1-eps)/(1+eps))^2))
        #println(violated_condition_lower_1, " ",violated_condition_upper_1," ",violated_condition_lower_2)
        if violated_condition_lower_1 || violated_condition_upper_1 || violated_condition_lower_2
            first_sample_lower = num_samples
            #println(" (top-k) lower increased at: "*string(first_sample_lower))
        else
            first_sample_upper = num_samples
            #println(" (top-k) upped increased at: "*string(first_sample_upper))
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("Maximum number of iterations "*string(last_stopping_samples))
    println("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stdout)
    next_stopping_samples::Int64 = trunc(Int,first_stopping_samples)
    prev_stopping_samples::Float64 = 0.0
    has_to_stop::Bool = false
    num_samples = 0
    sample_i::Int64 = 0
    while !has_to_stop
        sample_i = trunc(Int,next_stopping_samples-prev_stopping_samples)
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "rtb"
                    _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
            end
        end


        num_samples += sample_i
        println("num s ",num_samples)
        if num_samples >= omega
            has_to_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(num_samples)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end

        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            betweenness = reduce(+, local_temporal_betweenness)
            approx_top_k =  Array{Tuple{Int64,Float64}}([])

            wv = reduce(+,local_wv)
            r_mcrade = reduce(+,mcrade)
            for i in 1:tg.num_nodes
                push!(approx_top_k,(i,betweenness[i]/num_samples))
            end
            sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
            tmp_has_to_stop = Vector{Bool}([false])
            has_to_stop,top_k_result  = check_stopping_condition_topk(betweenness,k,union_sample,num_samples,eps,delta,iteration_index,mc_trials,r_mcrade,norm,partition_index, number_of_non_empty_partitions, non_empty_partitions,partitions_ids_map,wv,diam,avg_diam_ub)
            if has_to_stop
                println("Progressive sampler converged!")
                flush(stdout)
            else
                prev_stopping_samples = next_stopping_samples
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
                flush(stdout)
            end
        end 

    end
    println("(SH)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))

    return top_k_result,num_samples,time()-start_time
end

function threaded_progressive_cmcera_topk(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,k::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    ntasks = nthreads()
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    union_sample::Int64 = min(tg.num_nodes,max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/ntasks),k+20))
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    wv::Array{Float64} = Array{Float64}([])
    emp_w_node::Float64 = 0.0
    min_inv_w_node::Float64 = 0.0
    # Partitions
    number_of_non_empty_partitions::Int64 = 0.0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(tg.num_nodes)
    part_idx::Int64 = 1
    # TBC
    tmp_has_to_stop::Array{Bool} = Array{Bool}([false])
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating (sh)-Temporal Diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,0,0.9,false)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        flush(stdout)
    end
    start_time_bootstrap = time()

    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap using Variance")
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "rtb"
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
        end
    end

    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    println("Empirical peeling phase:")
    flush(stdout)
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0

    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i])
        max_wv = max(max_wv,wv[i])
        emp_w_node = wv[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    #=
    println("Number of non empty partitions "*string(number_of_non_empty_partitions))
    println("Bootstrap completed in "*string(round(time() - start_time; digits=4)))
    for key in keys(non_empty_partitions)
        println(" Part w. index "*string(key)*" has "*string(non_empty_partitions[key])*" elements, map to "*string(partitions_ids_map[key]))
    end
    flush(stdout)
    =#
    # Upper bound on the average distance
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
    # Upper bound on the top-1 temporal betweenness
    top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
    # define delta_for_progressive_bound
    println("Average (sh)-temporal path (upper bound) "*string(avg_diam_ub))
    # Upper limit on number of samples
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    prov = omega
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc estimation "*string(max_tbc))
    println("Sup empirical wimpy variance "*string(max_wv/tau))
    flush(stdout)
    iteration_index::Int64 =1 
    
    local_wv = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for _ in 1:ntasks]
    local_sp_lengths = [zeros(tg.num_nodes) for _ in 1:ntasks]
    
    omega = 10^15
    if max_num_samples > 0
        omega = max_num_samples
    end
    
    max_num_samples = omega
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    finish_bootstrap = string(round(time() - start_time_bootstrap; digits=4))
    println("Bootstrap completed in "*finish_bootstrap)
    println("Inferring initial sample size for the geometric sampler")
    flush(stdout)
    while first_sample_upper - first_sample_lower> 10
        num_samples = (first_sample_upper+first_sample_lower)÷2
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples) + log(2/delta)/num_samples/3
        if eps_guess > eps
            first_sample_lower = num_samples
        else
            first_sample_upper = num_samples
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stdout)
    next_stopping_samples::Float64 = first_stopping_samples
    prev_stopping_samples::Float64 = 0.0
    has_to_stop::Bool = false
    num_samples = 0
    sample_i::Int64 = 0
    while !has_to_stop
        approx_top_k =  Array{Tuple{Int64,Float64}}([])
        sample_i = trunc(Int,next_stopping_samples-prev_stopping_samples)
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "rtb"
                    _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
            end
        end


        num_samples += sample_i

        if num_samples >= omega
            has_to_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(num_samples)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end
    
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            
            betweenness = reduce(+, local_temporal_betweenness)
     
            wv = reduce(+,local_wv)
            #sp_lengths = reduce(+,local_sp_lengths) 
            r_mcrade = reduce(+,mcrade)
                
            tmp_omega = Vector{Float64}([omega])
            tmp_has_to_stop = Vector{Bool}([false])
            tmp_conv_number = Vector{Int64}([0])
            for u in 1:tg.num_nodes
                push!(approx_top_k,(u,betweenness[u]))
            end
            sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
            _check_stopping_condition_topk!(approx_top_k[begin:union_sample],betweenness,wv,last_stopping_samples,num_samples,eps,delta,iteration_index,true,diam,sp_lengths,num_samples,mc_trials,partition_index,partitions_ids_map,wv,r_mcrade,number_of_non_empty_partitions,tmp_omega,norm,tmp_has_to_stop,tmp_conv_number)
            println("Converged "*string(tmp_conv_number[1])*"/"*string(union_sample))
            omega = tmp_omega[1]
            has_to_stop = tmp_has_to_stop[1]
            if has_to_stop
                println("Progressive sampler converged!")
                flush(stdout)
            else
                prev_stopping_samples = next_stopping_samples
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
                flush(stdout)
            end
                    
        end


    end
    
    println("(SH)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    betweenness = reduce(+, local_temporal_betweenness)
    approx_top_k =  Array{Tuple{Int64,Float64}}([])
    for u in 1:tg.num_nodes
        push!(approx_top_k,(u,betweenness[u]/num_samples))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    return approx_top_k,num_samples,max_num_samples,time()-start_time


end
function threaded_progressive_cmcera_dep(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    wv::Array{Float64} = Array{Float64}([])
    emp_w_node::Float64 = 0.0
    min_inv_w_node::Float64 = 0.0
    # Partitions
    number_of_non_empty_partitions::Int64 = 0.0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(tg.num_nodes)
    part_idx::Int64 = 1
    # TBC
    tmp_has_to_stop::Array{Bool} = Array{Bool}([false])
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for i in 1:nthreads()]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0

    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step,0.9,false)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        diam+=1
    end
    
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0

    println("Bootstrap phase "*string(tau)*" iterations")
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
        elseif algo == "ob"
            _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
        elseif algo == "rtb"
            _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    println("Empirical peeling phase:")
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0

    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i])
        max_wv = max(max_wv,wv[i])
        emp_w_node = wv[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    println("Number of non empty partitions "*string(number_of_non_empty_partitions))
    println("Bootstrap completed in "*string(round(time() - start_time; digits=4)))
    for key in keys(non_empty_partitions)
        println(" Part w. index "*string(key)*" has "*string(non_empty_partitions[key])*" elements, map to "*string(partitions_ids_map[key]))
    end
   
    # Upper bound on the average distance
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
    # Upper bound on the top-1 temporal betweenness
    top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
    # define delta_for_progressive_bound
    println("AVERAGE DIAM UB "*string(avg_diam_ub))
    # Upper limit on number of samples
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    omega = 0.5/eps/eps * (log2(diam-2)+1+log(2/delta))
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc est "*string(max_tbc))
    println("Sup emp wimpy variance "*string(max_wv/tau))
    iteration_index::Int64 =1 
    
    local_wv = [zeros(tg.num_nodes) for i in 1:nthreads()]
    local_temporal_betweenness = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for i in 1:nthreads()]
    local_sp_lengths = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega = 10^15
    if max_num_samples > 0
        omega = max_num_samples
    end
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau

    while first_sample_upper - first_sample_lower> 10
        num_samples = (first_sample_upper+first_sample_lower)÷2
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples) + log(2/delta)/num_samples/3
        if eps_guess > eps
            first_sample_lower = num_samples
        else
            first_sample_upper = num_samples
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("First stopping samples "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("First stopping samples dropped to "*string(first_stopping_samples))
    end
    next_stopping_samples::Float64 = first_stopping_samples
    
    has_to_stop::Bool = false
    num_samples = 0
    while !has_to_stop
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            elseif algo == "rtb"
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            end
        end
        num_samples += sample_step
       
        if num_samples >= omega
            println("Num samples/ω : "*string(num_samples)*"/"*string(omega))
            has_to_stop = true
        end
        #println("Checking stopping condition")
        #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
        #println(" num_samples ",num_samples,"  ",next_stopping_samples)
        # & (num_samples >= next_stopping_samples)
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            betweenness = reduce(+, local_temporal_betweenness)
     
            wv = reduce(+,local_wv)
            sp_lengths = reduce(+,local_sp_lengths) 
            r_mcrade = reduce(+,mcrade)
            println("Checking stopping condition")
            #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
            #println(" num_samples ",num_samples,"  ",next_stopping_samples)
            tmp_omega = Vector{Float64}([omega])
            tmp_has_to_stop = Vector{Bool}([false])
            _check_stopping_condition!(betweenness,wv,last_stopping_samples,num_samples,eps,delta,iteration_index,true,diam,sp_lengths,num_samples,mc_trials,partition_index,partitions_ids_map,wv,r_mcrade,number_of_non_empty_partitions,tmp_omega,norm,tmp_has_to_stop)
            omega = tmp_omega[1]
            has_to_stop = tmp_has_to_stop[1]
            if has_to_stop
                println("Progressive sampler converged!")
            else
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
            end
                    
        end


    end

    println("(SH)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))

    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness.*[1/num_samples],num_samples,max_num_samples,time()-start_time


end


function _sh_accumulate_trk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SRTP_DS(tg.num_nodes, indexes)
    end
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    t_z::Int64 = tg.temporal_edges[lastindex(tg.temporal_edges)][3]+1
    index_z::Int64 = tn_index[(z,t_z)]
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64,Int64} = (-1, -1,-1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
    cur_w::Tuple{Int64,Int64} = (-1,-1)
    for u in 1:tg.num_nodes
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0,tni))
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = temporal_node[3]
        if bfs_ds.dist_t[tni] < d_z_min
            for neig in next_temporal_neighbors(tal, u, t)
                w = neig[1]
                t_w = neig[2]
                tni_w = tn_index[(w, t_w)]
                if bfs_ds.dist_t[tni_w] == -1
                    bfs_ds.dist_t[tni_w] = bfs_ds.dist_t[tni] + 1
                    if bfs_ds.dist[w] == -1
                        bfs_ds.dist[w] = bfs_ds.dist_t[tni] + 1

                        if w == z
                            d_z_min = bfs_ds.dist[w]
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, (w,t_w,tni_w))
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], (temporal_node[1],temporal_node[2],tni))
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
                            return [], 0.0
                        end
                        bfs_ds.sigma[w] += bfs_ds.sigma_t[tni]
                    end
                    if w == z
                        push!(bfs_ds.predecessors[index_z], (neig[1],neig[2],tni_w))
                    end
                end
            end
        end
    end
    if bfs_ds.dist[z] > 0
        sp_lengths[bfs_ds.dist[z]+1] +=1
        totalWeight = 0
        randomEdge = 0
        curWeight = 0
        totalWeight = bfs_ds.sigma[z]
        cur_w = (z,t_z)
        tni = index_z
        while cur_w[1] != s
            if cur_w == (z,t_z)
                totalWeight = bfs_ds.sigma[z]
            else
                totalWeight = bfs_ds.sigma_t[tni]
            end
            randomEdge = rand(0:totalWeight-1)
            curEdge = 0
            for pred in bfs_ds.predecessors[tni]
                pred_i = pred[3]
                curEdge += bfs_ds.sigma_t[pred_i]
                cur_w = (pred[1],pred[2])
                tni = pred_i
                if curEdge > randomEdge
                    if pred[1]!= s && pred[1] != z
                        temporal_betweenness_centrality[pred[1]] += 1
                        wimpy_variance[pred[1]] += 1
                        if !boostrap_phase
                            #mcrade deve essere usato da ogni thread in modo indipendente
                            v_idx = pred[1]*mc_trials
                            for j in 1:mc_trials
                                mcrade[v_idx + j] += lambdas[j] * 1
                            end
                        end
                    end
                    break
                end
            end
            
        end
        
    end
    if (bigint)
        bfs_ds = BI_BFS_SRTP_DS(0,0)
    else
        bfs_ds = BFS_SRTP_DS(0,0)
    end
    return nothing

end





#=
ONBRA
=#

function _sh_accumulate_onbra!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    end
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min
            for neig in next_temporal_neighbors(tal, u, t)
                w = neig[1]
                t_w = neig[2]
                tni_w = tn_index[(w, t_w)]
                if bfs_ds.dist_t[tni_w] == -1
                    bfs_ds.dist_t[tni_w] = bfs_ds.dist_t[tni] + 1
                    if bfs_ds.dist[w] == -1
                        bfs_ds.dist[w] = bfs_ds.dist_t[tni] + 1

                        if w == z
                            d_z_min = bfs_ds.dist[w]
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, neig)
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], temporal_node)
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
                            return [], 0.0
                        end
                        bfs_ds.sigma[w] += bfs_ds.sigma_t[tni]
                    end
                end
            end
        end
    end
    if bfs_ds.dist[z] > 0
        sp_lengths[bfs_ds.dist[z]+1] +=1
        for tn in 1:lastindex(bfs_ds.dist_t)
            bfs_ds.sigma_z[tn] = 0
            bfs_ds.boolean_matrix[tn] = false
        end
        tni = tn_index[(s, 0)]
        bfs_ds.sigma_z[tni] = 1
    
        for t in 1:lastindex(tg.file_time)
            tni = get(tn_index, (z, t), 0)
            if tni > 0 && bfs_ds.sigma_t[tni] > 0
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] == typemax(UInt128))
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += 1
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
        while length(bfs_ds.backward_queue) > 0
            temporal_node = dequeue!(bfs_ds.backward_queue)
            tni = tn_index[(temporal_node[1], temporal_node[2])]
            if temporal_node[1] != s
                temporal_betweenness_centrality[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                wimpy_variance[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))^2
                if !boostrap_phase
                    #mcrade deve essere usato da ogni thread in modo indipendente
                    v_idx = temporal_node[1]*mc_trials
                    for j in 1:mc_trials
                        mcrade[v_idx + j] += lambdas[j] * (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                    end
                end
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] > typemax(UInt128) - bfs_ds.sigma_z[tni])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += bfs_ds.sigma_z[tni]
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(0, 0)
    else
        bfs_ds = BFS_ONBRA_DS(0, 0)
    end
    return nothing

end


# RTB


function _sh_accumulate_rtb!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})

    if (bigint)
        bfs_ds = BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_DS(tg.num_nodes, length(keys(tn_index)))
    end
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    v::Int64 = -1
    t_v::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:(length(tn_index))
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.delta_sh[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    sp_lengths[bfs_ds.dist[s]+1]+=1
    enqueue!(bfs_ds.queue, (s, 0))
    while length(bfs_ds.queue) != 0
        temporal_node = dequeue!(bfs_ds.queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        for neig in next_temporal_neighbors(tal, u, t)
            w = neig[1]
            t_w = neig[2]
            tni_w = tn_index[(w, t_w)]
            if bfs_ds.dist_t[tni_w] == -1
                bfs_ds.dist_t[tni_w] = bfs_ds.dist_t[tni] + 1
                if bfs_ds.dist[w] == -1
                    bfs_ds.dist[w] = bfs_ds.dist_t[tni_w]
                    sp_lengths[bfs_ds.dist[w]+1] += 1
                end
                enqueue!(bfs_ds.queue, neig)
                push!(bfs_ds.stack, neig)
            end
            if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                    println("Overflow occurred with source ", s)
                    return [], 0.0
                end
                bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                push!(bfs_ds.predecessors[tni_w], temporal_node)
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                        println("Overflow occurred with source ", s)
                        return [], 0.0
                    end
                    bfs_ds.sigma[w] += bfs_ds.sigma_t[tni]
                end
            end
        end
    end
    temporal_betweenness_centrality[s] -= (1/(tg.num_nodes-1))* (count(x -> x >= 0, bfs_ds.dist) - 1)
    while length(bfs_ds.stack) != 0
        temporal_node = pop!(bfs_ds.stack)
        w = temporal_node[1]
        t_w = temporal_node[2]
        tni_w = tn_index[(w, t_w)]
        if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
            bfs_ds.delta_sh[tni_w] += bfs_ds.sigma_t[tni_w] / bfs_ds.sigma[w]
        end
        for pred in bfs_ds.predecessors[tni_w]
            v = pred[1]
            t_v = pred[2]
            tni_v = tn_index[(v, t_v)]
            bfs_ds.delta_sh[tni_v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w]
            temporal_betweenness_centrality[v] += (1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w]
            wimpy_variance[v] += ((1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w])^2
            if !boostrap_phase
                #mcrade deve essere usato da ogni thread in modo indipendente
                v_idx = v*mc_trials
                for j in 1:mc_trials
                    mcrade[v_idx + j] += lambdas[j] * ((1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w])
                end
            end
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_DS(0,0)
    else
        bfs_ds = BFS_DS(0,0)
    end
    return nothing
end



function _sh_accumulate_onbra_rel_topk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64},distinct_nodes_top_k::Array{Int64})
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    end
    topk_bootstrap::Array{Int64} = zeros(tg.num_nodes)
    apx_to_add_bootstrap::Array{Int8} = zeros(tg.num_nodes)

    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min
            for neig in next_temporal_neighbors(tal, u, t)
                w = neig[1]
                t_w = neig[2]
                tni_w = tn_index[(w, t_w)]
                if bfs_ds.dist_t[tni_w] == -1
                    bfs_ds.dist_t[tni_w] = bfs_ds.dist_t[tni] + 1
                    if bfs_ds.dist[w] == -1
                        bfs_ds.dist[w] = bfs_ds.dist_t[tni] + 1

                        if w == z
                            d_z_min = bfs_ds.dist[w]
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, neig)
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], temporal_node)
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
                            return [], 0.0
                        end
                        bfs_ds.sigma[w] += bfs_ds.sigma_t[tni]
                    end
                end
            end
        end
    end
    if bfs_ds.dist[z] > 0
        sp_lengths[bfs_ds.dist[z]+1] +=1
        for tn in 1:lastindex(bfs_ds.dist_t)
            bfs_ds.sigma_z[tn] = 0
            bfs_ds.boolean_matrix[tn] = false
        end
        tni = tn_index[(s, 0)]
        bfs_ds.sigma_z[tni] = 1
    
        for t in 1:lastindex(tg.file_time)
            tni = get(tn_index, (z, t), 0)
            if tni > 0 && bfs_ds.sigma_t[tni] > 0
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] == typemax(UInt128))
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += 1
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
        while length(bfs_ds.backward_queue) > 0
            temporal_node = dequeue!(bfs_ds.backward_queue)
            tni = tn_index[(temporal_node[1], temporal_node[2])]
            if temporal_node[1] != s
                temporal_betweenness_centrality[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                wimpy_variance[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))^2
                topk_bootstrap[temporal_node[1]] += 1
                if (topk_bootstrap[temporal_node[1]] >= 3 && apx_to_add_bootstrap[temporal_node[1]] == 0)
                    apx_to_add_bootstrap[temporal_node[1]] = 1
                    distinct_nodes_top_k[1]+=1
                end
                if !boostrap_phase
                    #mcrade deve essere usato da ogni thread in modo indipendente
                    v_idx = temporal_node[1]*mc_trials
                    for j in 1:mc_trials
                        mcrade[v_idx + j] += lambdas[j] * (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                    end
                end
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] > typemax(UInt128) - bfs_ds.sigma_z[tni])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += bfs_ds.sigma_z[tni]
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(0, 0)
    else
        bfs_ds = BFS_ONBRA_DS(0, 0)
    end
    return nothing

end