function threaded_trk_shortest_foremost(tg::temporal_graph,sample_size::Int64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}

    start_time = time()
    sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    z::Int64 = 0
    println("Using ",nthreads()," Trheads")

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i][1]
        z = sample[i][2]
        _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TRK-SFM. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size]
    return betweenness,time()-start_time
end


function _trk_sfm_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,temporal_betweenness_centrality::Vector{Float64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
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
        
        totalWeight = 0
        randomEdge = 0
       
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
                    end
                    break
                end
            end
        end
       
    end
    return nothing
end







function threaded_progressive_trk_shortest_foremost(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)

    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
    else
        omega = trunc(Int,(1.0/(2*eps^2))*log2(2*tg.num_nodes/delta))
    end
    println("Maximum sample size "*string(omega))
    tau::Int64 = trunc(Int64,omega/start_factor)
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap phase "*string(tau)*" iterations")
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    eps_lb::Array{Float64} = zeros(tg.num_nodes)
    eps_ub::Array{Float64} = zeros(tg.num_nodes)
    delta_lb_min_guess::Array{Float64} = [0.0]
    delta_ub_min_guess::Array{Float64} = [0.0]
    delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
    delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
    _compute_δ_guess!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess) 
    println("Bootstrap completed ")

    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    while sampled_so_far <= omega && !stop[1]
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        end
        sampled_so_far += sample_step
        betweenness = reduce(+, local_temporal_betweenness)
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-TRK-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
    end
    return betweenness,eps_lb,eps_ub,sampled_so_far,time()-start_time


end