function progressive_cmcera(tg::temporal_graph,eps::Float64,delta::Float64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,force_gc::Bool = false)
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
    betweenness::Array{Float64} = zeros(Float64,tg.num_nodes)
    r_mcrade::Array{Float64} = zeros(Float64,(tg.num_nodes+1)*mc_trials)
    sp_lengths::Array{Int64} = zeros(Int64,tg.num_nodes)
    wv::Array{Float64} = zeros(Float64,tg.num_nodes)

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
    if vc_upper_bund == true
        println("Bootstrap using VC dimension")
    else
        println("Bootstrap using Variance")
    end
    flush(stdout)
    for _ in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,betweenness,wv,r_mcrade,sp_lengths)
        elseif algo == "ob"
            _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,betweenness,wv,r_mcrade,sp_lengths)
        elseif algo == "rtb"
            _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,betweenness,wv,r_mcrade,sp_lengths)
        end
    end
    
    println("Empirical peeling phase:")
    flush(stdout)
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0
    
    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i]/tau)
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
    betweenness = zeros(Float64,tg.num_nodes)
    r_mcrade = zeros(Float64,(tg.num_nodes+1)*mc_trials)
    sp_lengths = zeros(Int64,tg.num_nodes)
    wv = zeros(Float64,tg.num_nodes)
    next_stopping_samples::Float64 = first_stopping_samples
    prev_stopping_samples::Float64 = 0.0
    has_to_stop::Bool = false
    num_samples = 0
    sample_i::Int64 = 0
    number_of_gc_calls::Int64 = 0
    overall_cleaning_time::Float64 = 0.0
    while !has_to_stop
        sample_i = trunc(Int,next_stopping_samples-prev_stopping_samples)
        for _ in 1:sample_i
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,betweenness,wv,r_mcrade,sp_lengths)
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,betweenness,wv,r_mcrade,sp_lengths)
            elseif algo == "rtb"
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,betweenness,wv,r_mcrade,sp_lengths)
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
            
            


            #betweenness = reduce(+, local_temporal_betweenness)
     
            #wv = reduce(+,local_wv)
            #sp_lengths = reduce(+,local_sp_lengths) 
            #r_mcrade = reduce(+,mcrade)
                
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
    return betweenness.*[1/num_samples],num_samples,max_num_samples,time()-start_time


end







