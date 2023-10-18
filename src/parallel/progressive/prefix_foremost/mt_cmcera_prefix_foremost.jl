
function threaded_progressive_cmcera_prefix_foremost(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,algo::String = "trk",diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10)
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
        println("Approximating (pfm)-Temporal Diameter ")
        _,_,_,_,_,diam_v,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step,0.9)
        diam = trunc(Int,diam_v)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
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
                _pfm_accumulate_trk!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            else
                _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
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
    println("Average (pfm)-temporal path (upper bound) "*string(avg_diam_ub))
    # Upper limit on number of samples
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc est "*string(max_tbc))
    println("Sup emp wimpy variance "*string(max_wv/tau))
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
        sample_i = trunc(Int,next_stopping_samples)-trunc(Int,prev_stopping_samples)
        #println("Iteration "*string(iteration_index)*" Sample size "*string(sample_i))
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _pfm_accumulate_trk!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                else
                    _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
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
        
        if !has_to_stop & (num_samples < trunc(Int,last_stopping_samples))&(num_samples >= trunc(Int,next_stopping_samples))
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
            end
        end


    end

    println("(PFM)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness.*[1/num_samples],num_samples,max_num_samples,time()-start_time


end


function _pfm_accumulate_trk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    bigint::Bool  = false
    bfs_ds = BFS_PFM_SRTP_DS(tg.num_nodes)
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    cur_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
    dist::Array{Int64} = zeros(tg.num_nodes)
    for u in 1:tg.num_nodes
        bfs_ds.t_hit[u] = -1
        bfs_ds.sigma[u] = 0
        bfs_ds.predecessors[u] = Set{Int64}()
    end
    
    bfs_ds.sigma[s] = 1
    bfs_ds.t_hit[s] = 0
    dist[s] = 0
    for tn in tal[s]
        enqueue!(bfs_ds.priority_queue,(s,tn[1],tn[2]),tn[2])
    end
    t_z_min::Float64 = Inf
    while length(bfs_ds.priority_queue) != 0
        temporal_edge = dequeue!(bfs_ds.priority_queue)
        v = temporal_edge[1]
        w = temporal_edge[2]
        t_w = temporal_edge[3]
        if t_w < t_z_min
            if bfs_ds.t_hit[w] == -1
                bfs_ds.t_hit[w] = t_w
                dist[w] = dist[v] +1
                for neig in next_temporal_neighbors(tal,w,t_w)
                    enqueue!(bfs_ds.priority_queue,(w,neig[1],neig[2]),neig[2])
                end
                if z == w &&  t_z_min > t_w
                    t_z_min = t_w
                end
            end
            if bfs_ds.t_hit[w] == t_w
                if (!bigint && bfs_ds.sigma[v] > typemax(UInt128) - bfs_ds.sigma[w])
                    println("Overflow occurred with source ", s)
                    return [], 0.0
                end
                bfs_ds.sigma[w] += bfs_ds.sigma[v]
                push!(bfs_ds.predecessors[w], v)
            end
        end
    end
    if bfs_ds.sigma[z] > 0
        sp_lengths[dist[z]+1] +=1
        cur_w = z
        while cur_w != s
            totalWeight = 0
            if cur_w == z
                totalWeight = bfs_ds.sigma[z]
            else
                for pred in bfs_ds.predecessors[cur_w]
                    totalWeight+=bfs_ds.sigma[pred]
                end
            end
            randomEdge = rand(0:totalWeight-1)
            curEdge = 0
            for pred in bfs_ds.predecessors[cur_w]
                curEdge+= bfs_ds.sigma[pred]
                if curEdge > randomEdge
                    cur_w = pred
                    if pred[1] != s && pred[1] != z
                        temporal_betweenness_centrality[pred] += 1
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
    bfs_ds = BFS_PFM_SRTP_DS(0)

    return nothing

end





#=
ONBRA
=#

function _pfm_accumulate_onbra!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    bigint::Bool  = false
    bfs_ds::BFS_ONBRA_PFM_DS = BFS_ONBRA_PFM_DS(tg.num_nodes)
    
    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    dist::Array{Int64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    for u in 1:tg.num_nodes
        bfs_ds.t_hit[u] = -1
        bfs_ds.sigma[u] = 0
        bfs_ds.predecessors[u] = Set{Int64}()
    end
    
    bfs_ds.sigma[s] = 1
    bfs_ds.t_hit[s] = 0
    dist[s] = 0
    for tn in tal[s]
        enqueue!(bfs_ds.priority_queue,(s,tn[1],tn[2]),tn[2])
    end
    t_z_min = Inf
    while length(bfs_ds.priority_queue) != 0
        temporal_edge = dequeue!(bfs_ds.priority_queue)
        v = temporal_edge[1]
        w = temporal_edge[2]
        t_w = temporal_edge[3]
        if t_w < t_z_min
            if bfs_ds.t_hit[w] == -1
                bfs_ds.t_hit[w] = t_w
                dist[w] = dist[v] +1
                for neig in next_temporal_neighbors(tal,w,t_w)
                    enqueue!(bfs_ds.priority_queue,(w,neig[1],neig[2]),neig[2])
                end
                if z == w &&  t_z_min > t_w
                    t_z_min = t_w
                end
            end
            if bfs_ds.t_hit[w] == t_w
                if (!bigint && bfs_ds.sigma[v] > typemax(UInt128) - bfs_ds.sigma[w])
                    println("Overflow occurred with source ", s)
                    return [], 0.0
                end
                bfs_ds.sigma[w] += bfs_ds.sigma[v]
                push!(bfs_ds.predecessors[w], v)
            end
        end
    end
    if bfs_ds.sigma[z] > 0
        sp_lengths[dist[z]+1] +=1

        for v in 1:lastindex(bfs_ds.sigma)
            bfs_ds.sigma_z[v] = 0
            bfs_ds.boolean_matrix[v] = false
        end
        #bfs_ds.sigma_z[s] = 1
    
        for pred in bfs_ds.predecessors[z]
            bfs_ds.sigma_z[pred] =1
            if !bfs_ds.boolean_matrix[pred]
                enqueue!(bfs_ds.backward_queue,pred)
                bfs_ds.boolean_matrix[pred] = true
            end
        end
        while length(bfs_ds.backward_queue) > 0
            v = dequeue!(bfs_ds.backward_queue)
            if v != s
                temporal_betweenness_centrality[v] += (bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z]))
                wimpy_variance[v] += ((bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z])))^2
                if !boostrap_phase
                    #mcrade deve essere usato da ogni thread in modo indipendente
                    v_idx = v*mc_trials
                    for j in 1:mc_trials
                        mcrade[v_idx + j] += lambdas[j] * (bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z]))
                    end
                end
                for pred in bfs_ds.predecessors[v]
                    if (!bigint && bfs_ds.sigma_z[pred] > typemax(UInt128) - bfs_ds.sigma_z[pred])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[pred] += bfs_ds.sigma_z[v]
                    if !bfs_ds.boolean_matrix[pred]
                        enqueue!(bfs_ds.backward_queue, pred) 
                        bfs_ds.boolean_matrix[pred] = true
                    end                    
                end
            end
        end
    end
    bfs_ds = BFS_ONBRA_PFM_DS(0)

    return nothing

end


# RTB


function _pfm_accumulate_rtb!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    bigint::Bool = false
    bfs_ds::BFS_prefix_foremost_betweenness = BFS_prefix_foremost_betweenness(tg.num_nodes)

    lambdas::Array{Int64} = zeros(mc_trials)
    maxval_lambdas::Int64 = 100000000
    maxval_half::Int64 = maxval_lambdas/2
    if !boostrap_phase
        for j in 1:mc_trials
            lambdas[j] = (sample(0:maxval_lambdas-1) >= maxval_half)*2-1
        end
    end
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    dist::Array{Int64} = zeros(tg.num_nodes)
    for u in 1:tg.num_nodes
        bfs_ds.delta[u] = 1
        bfs_ds.sigma[u] = 0
        bfs_ds.predecessors[u] = Set{Int64}()
        bfs_ds.t_min[u] = -1
    end
    bfs_ds.t_min[s] = 0
    bfs_ds.sigma[s] = 1
    dist[s] = 0
    for tn in tal[s]
        enqueue!(bfs_ds.priority_queue,(s,tn[1],tn[2]),tn[2])
    end
    while length(bfs_ds.priority_queue) != 0
        temporal_edge = dequeue!(bfs_ds.priority_queue)
        v = temporal_edge[1]
        w = temporal_edge[2]
        t_w = temporal_edge[3]
        if bfs_ds.t_min[w] == -1
            bfs_ds.t_min[w] = t_w
            dist[w] = dist[v] + 1
            sp_lengths[dist[w]+1]+=1
            push!(bfs_ds.stack,w)
            for neig in next_temporal_neighbors(tal,w,t_w)
                enqueue!(bfs_ds.priority_queue,(w,neig[1],neig[2]),neig[2])
            end
        end
        if bfs_ds.t_min[w] == t_w
            if (bfs_ds.sigma[v] > typemax(UInt128) - bfs_ds.sigma[w])
                println("Overflow occurred with source ", s)
                return [], 0.0
            end
            bfs_ds.sigma[w] += bfs_ds.sigma[v]
            push!(bfs_ds.predecessors[w], v)
        end
    end


    temporal_betweenness_centrality[s] -=(1/(tg.num_nodes-1))* (count(x -> x >= 0, bfs_ds.t_min) - 1)
    while length(bfs_ds.stack) != 0
        w = pop!(bfs_ds.stack)
        for v in bfs_ds.predecessors[w]
            summand = (Float64(bfs_ds.sigma[v] / bfs_ds.sigma[w])) * bfs_ds.delta[w]
            bfs_ds.delta[v] += summand
            temporal_betweenness_centrality[v] += (1/(tg.num_nodes-1))*summand
            wimpy_variance[v] += ((1/(tg.num_nodes-1))*summand)^2
            if !boostrap_phase
                #mcrade deve essere usato da ogni thread in modo indipendente
                v_idx = v*mc_trials
                for j in 1:mc_trials
                    mcrade[v_idx + j] += lambdas[j] * ((1/(tg.num_nodes-1))*summand)
                end
            end
        end
    end
    bfs_ds = BFS_prefix_foremost_betweenness(0)

    return nothing
end