function threaded_progressive_silvan_prefix_foremost_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,verbose_step::Int64,diam::Int64 = -1,algo::String = "trk",empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    union_sample::Int64 = min(tg.num_nodes,max(sqrt(lastindex(tg.temporal_edges))/nthreads(),k+20))
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
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
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for i in 1:nthreads()]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0
    top_k_result::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
    #diam::Float64 = 0.0
    if (diam == -1) 
        println("Approximating diameter ")
        _,_,_,_,_,diam_v,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step,0.9)
        diam = trunc(Int,diam_v)
    end
    #max_num_samples = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0
    bootstrap_finished::Bool = false
    approx_to_add::Array{Int64} = zeros(tg.num_nodes)
    distinct_nodes_top_k::Array{Int64} = [0]
    tbc::Array{Float64} = zeros(tg.num_nodes)
    println("Bootstrap phase "*string(tau)*" iterations")
    n_pairs = 0
    while !bootstrap_finished
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _pfm_accumulate_trk_topk!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],distinct_nodes_top_k,approx_to_add)
            elseif algo == "ob"
                _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            else
                _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            end
        end
        #betweenness = reduce(+, local_temporal_betweenness)
        #betweenness = betweenness .* [1/tau]
        #
        n_pairs += sample_step
        bootstrap_finished = (distinct_nodes_top_k[1] >= 1.5*k )&&(n_pairs >= tau)
        println(" distinct_nodes_top_k ",distinct_nodes_top_k[1])
    end
    tau = n_pairs
    betweenness = tbc .*[1/tau]
    println("Empirical peeling phase:")
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0
    wv = reduce(+,local_wv)

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
    #=
    # Upper bound on the average distance
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
    # Upper bound on the top-1 temporal betweenness
    top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
    # define delta_for_progressive_bound
    println("AVERAGE DIAM UB "*string(avg_diam_ub))
    # Upper limit on number of samples
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc est "*string(max_tbc))
    println("Sup emp wimpy variance "*string(max_wv/tau))
    =#
    iteration_index::Int64 =1 
    est_kth_bc::Float64 = 10.0
    for i in 1:tg.num_nodes
        push!(approx_top_k,(i,betweenness[i]))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    est_kth_bc = approx_top_k[k][2] +1/tau
    approx_top_k =  Array{Tuple{Int64,Float64}}([])
    println("Est k-th tbc "*string(est_kth_bc))
    local_wv = [zeros(tg.num_nodes) for i in 1:nthreads()]
    local_temporal_betweenness = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for i in 1:nthreads()]
    local_sp_lengths = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega = compute_vapnik_chervonenkis_bound(trunc(Int,diam_v),eps,delta,0.5)
    if max_num_samples > 0
        omega = max_num_samples
    end
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    println("First sample lower "*string(first_sample_lower)*"First sample upper "*string(first_sample_upper))
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    violated_condition_lower_1::Bool = false
    violated_condition_lower_2::Bool = false
    violated_condition_upper_1::Bool = false
    while first_sample_upper - first_sample_lower > 10
        num_samples = (first_sample_upper+first_sample_lower)/2
        eps_guess = sqrt(2*est_kth_bc*log(2/delta/est_kth_bc) /num_samples) + log(2/delta/est_kth_bc)/num_samples/3
        println("Eps guess "*string(eps_guess))
        violated_condition_lower_1 = (est_kth_bc + eps_guess > est_kth_bc/(1-eps))
        violated_condition_upper_1 = (est_kth_bc - eps_guess < est_kth_bc/(1+eps))
        violated_condition_lower_2 = (eps_guess >= (est_kth_bc *eps*((1-eps)/(1+eps))^2))
        println(violated_condition_lower_1, " ",violated_condition_upper_1," ",violated_condition_lower_2)
        if violated_condition_lower_1 || violated_condition_upper_1 || violated_condition_lower_2
            first_sample_lower = num_samples
            println(" (top-k) lower increased at: "*string(first_sample_lower))
        else
            first_sample_upper = num_samples
            println(" (top-k) upped increased at: "*string(first_sample_upper))
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("First stopping samples "*string(first_stopping_samples)* " Last stopping samples "*string(last_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("First stopping samples dropped to "*string(first_stopping_samples))
    end
    next_stopping_samples::Float64 = first_stopping_samples
    
    has_to_stop::Bool = false
    num_samples::Float64 = 0.0
    while !has_to_stop
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _pfm_accumulate_trk!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            elseif algo == "ob"
                _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            else
                _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            end
        end
        num_samples += sample_step
        #println("Num s "*string(num_samples)*" "*string(has_to_stop)* " "*string((num_samples < last_stopping_samples))*"  "*string((num_samples >= next_stopping_samples)))
        #println(" "*string(num_samples)*" "*string(last_stopping_samples)*" FRST = "*string(first_stopping_samples))
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
            approx_top_k =  Array{Tuple{Int64,Float64}}([])

            wv = reduce(+,local_wv)
            sp_lengths = Array{Float64}(reduce(+,local_sp_lengths))
            r_mcrade = reduce(+,mcrade)
            println("Checking stopping condition")
            for i in 1:tg.num_nodes
                push!(approx_top_k,(i,betweenness[i]/num_samples))
            end
            sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
            #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
            #println(" num_samples ",num_samples,"  ",next_stopping_samples)
            #(approx_top_k,union_sample,num_samples,eps,delta,iteration_index,partition_index,partitions_ids_map,number_of_non_empty_partitions,mc_trials,sup_bcest_partition,sup_empwvar_partition,norm)
            has_to_stop,top_k_result  = check_stopping_condition_topk(betweenness,k,union_sample,num_samples,eps,delta,iteration_index,mc_trials,r_mcrade,norm,partition_index, number_of_non_empty_partitions, non_empty_partitions,partitions_ids_map,wv,diam,sp_lengths)

            if has_to_stop
                println("Progressive sampler converged!")
            else
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
            end
                    
        end


    end

    println("(PFM)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))

    return top_k_result,num_samples,time()-start_time


end



function _pfm_accumulate_trk_topk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},distinct_nodes_top_k::Vector{Int64},approx_to_add::Array{Int64})
    bigint::Bool  = false
    bfs_ds = BFS_PFM_SRTP_DS(tg.num_nodes)
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
                        lk = ReentrantLock()
                        lock(lk)
                        try 
                            temporal_betweenness_centrality[pred[1]] += 1
                       
                            if (approx_to_add[pred[1]] == 0) && (temporal_betweenness_centrality[pred[1]] >= 3)
                                approx_to_add[pred[1]] = 1
                                distinct_nodes_top_k[1] +=1
                            end
                        finally 
                            unlock(lk)
                        end
                    end
                    break
                end
            end

        end    
    end
   
    return nothing

end

