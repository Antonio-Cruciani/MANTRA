
function threaded_progressive_cmcera_shortest_foremost(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10,hb::Bool = false)
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
    end    # Wimpy variance
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
        println("Approximating diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_foremost_diameter(tg,64,verbose_step,0.9)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        diam+=1
        flush(stdout)
    end
    
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0

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
                _sfm_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _sfm_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "rtb"
                _sfm_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
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
    println("Number of non empty partitions "*string(number_of_non_empty_partitions))
    println("Bootstrap completed in "*string(round(time() - start_time; digits=4)))
    for key in keys(non_empty_partitions)
        println(" Part w. index "*string(key)*" has "*string(non_empty_partitions[key])*" elements, map to "*string(partitions_ids_map[key]))
    end
    flush(stdout)
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
    flush(stdout)
    next_stopping_samples::Float64 = first_stopping_samples
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
                    _sfm_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sfm_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "rtb"
                    _sfm_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
            end
        end

        #=
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
        =#
        num_samples += sample_i
       
        if num_samples >= omega
            println("Num samples/ω : "*string(num_samples)*"/"*string(omega))
            has_to_stop = true
            flush(stdout)
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
            flush(stdout)
            #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
            #println(" num_samples ",num_samples,"  ",next_stopping_samples)
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

    println("(SFM)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))
    flush(stdout)
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness.*[1/num_samples],num_samples,max_num_samples,time()-start_time


end


function _sfm_accumulate_trk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
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
        bfs_ds.t_hit[u] = -1
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
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0,tni))
    d_z_min = Inf
    t_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = temporal_node[3]
        if bfs_ds.dist_t[tni] < d_z_min && t < t_z_min
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
                            if t_z_min > t_w
                                t_z_min = t_w
                            end
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
                if bfs_ds.t_hit[w] == -1 || bfs_ds.t_hit[w] > t_w
                    bfs_ds.t_hit[w] = t_w
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
        bfs_ds = BI_BFS_SFM_SRTP_DS(0,0)
    else
        bfs_ds = BFS_SFM_SRTP_DS(0,0)
    end
    return nothing

end





#=
ONBRA
=#

function _sfm_accumulate_onbra!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_SFM_DS(tg.num_nodes, length(keys(tn_index)))
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
        bfs_ds.t_hit[u] = -1
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
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    t_z_min = Inf
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min && t < t_z_min
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
                            if t_z_min > t_w
                                t_z_min = t_w
                            end
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

                if bfs_ds.t_hit[w] == -1 || bfs_ds.t_hit[w] > t_w
                    bfs_ds.t_hit[w] = t_w
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
    end
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
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(0,0)
    else
        bfs_ds = BFS_ONBRA_SFM_DS(0,0)
    end
    return nothing

end


# RTB


function _sfm_accumulate_rtb!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},mcrade::Vector{Float64},sp_lengths::Vector{Int64})

    if (bigint)
        bfs_ds = BI_BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BI_BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
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
        bfs_ds.t_hit[u] = -1
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:(length(tn_index))
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.delta_fm[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    bfs_ds.t_hit[s] = 0
    sp_lengths[bfs_ds.dist[s]+1] +=1
    enqueue!(bfs_ds.queue, (s, 0))
    iter = 0
    while length(bfs_ds.queue) != 0
        iter += 1
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
                    sp_lengths[bfs_ds.dist[w]+1]+=1
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
            if bfs_ds.t_hit[w] == -1 || bfs_ds.t_hit[w] > t_w
                bfs_ds.t_hit[w] = t_w
            end
        end
    end
    temporal_betweenness_centrality[s] -= (1/(tg.num_nodes-1))* (count(x -> x >= 0, bfs_ds.dist) - 1)
    while length(bfs_ds.stack) != 0
        temporal_node = pop!(bfs_ds.stack)
        w = temporal_node[1]
        t_w = temporal_node[2]
        tni_w = tn_index[(w, t_w)]
        if bfs_ds.t_hit[w] == t_w
            bfs_ds.delta_fm[tni_w] += 1
        end
        for pred in bfs_ds.predecessors[tni_w]
            v = pred[1]
            t_v = pred[2]
            tni_v = tn_index[(v, t_v)]
            bfs_ds.delta_fm[tni_v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
            temporal_betweenness_centrality[v] += (1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
            wimpy_variance[v] += ((1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w])^2
            if !boostrap_phase
                #mcrade deve essere usato da ogni thread in modo indipendente
                v_idx = v*mc_trials
                for j in 1:mc_trials
                    mcrade[v_idx + j] += lambdas[j] * ((1/(tg.num_nodes-1))*(bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w])
                end
            end
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_SFM_DS(0,0)
    else
        bfs_ds = BI_BFS_SFM_DS(0,0)
    end
    return nothing
end