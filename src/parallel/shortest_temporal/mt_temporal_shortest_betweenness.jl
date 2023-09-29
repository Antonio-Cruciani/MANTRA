

function distributed_temporal_shortest_betweenness(tg::temporal_graph,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    #local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]


    betweenness = @distributed (+) for s in 1:tg.num_nodes
        temp_betweenness = zeros(tg.num_nodes)
        _sstp_accumulate!(tg,tal,tn_index,s,bigint,temp_betweenness)
        #=
        if (verbose_step > 0 && s % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            time_to_finish = string(round((tg.num_nodes*(time() - start_time) / s )-(time() - start_time) ; digits=4))
            println("TSB. Processed " * string(s) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
            flush(stdout)
        end
        =#
        temp_betweenness
    end
 

   
    return betweenness,time()-start_time
end


function threaded_temporal_shortest_betweenness_prova(tg::temporal_graph,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    processed_so_far::Int64 = 0
    println("Using ",nthreads()," Trheads")
    
    
    vs_active = [i for i in 1:tg.num_nodes]
    d, r = divrem(tg.num_nodes, nthreads())
    ntasks = d == 0 ? r : nthreads()
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    if (bigint)
        bfs_ds = [BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))  for _ in 1:ntasks ]
    else
        bfs_ds = [BFS_DS(tg.num_nodes, length(keys(tn_index))) for _ in 1:ntasks ]
    end
    task_size = cld(tg.num_nodes, ntasks)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tg.num_nodes, task_size))
        Threads.@spawn for s in @view(vs_active[task_range])
            _sstp_accumulate_prova!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t],bfs_ds[t])
            processed_so_far = processed_so_far + 1
            if (verbose_step > 0 && processed_so_far % verbose_step == 0)
                finish_partial::String = string(round(time() - start_time; digits=4))
                time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
                println("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
                flush(stdout)
            end
        end
    end
    #=
    Base.Threads.@threads for s in 1:tg.num_nodes
        _sstp_accumulate_prova!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()],bfs_ds[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
            flush(stdout)
        end
    end
    =#
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness,time()-start_time
end

function threaded_temporal_shortest_betweenness(tg::temporal_graph,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    processed_so_far::Int64 = 0
    println("Using ",nthreads()," Trheads")
    vs_active = [i for i in 1:tg.num_nodes]
    d, r = divrem(tg.num_nodes, nthreads())
    ntasks = d == 0 ? r : nthreads()
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    task_size = cld(tg.num_nodes, ntasks)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tg.num_nodes, task_size))
        Threads.@spawn for s in @view(vs_active[task_range])
            _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t])
            processed_so_far = processed_so_far + 1
            if (verbose_step > 0 && processed_so_far % verbose_step == 0)
                finish_partial::String = string(round(time() - start_time; digits=4))
                time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
                println("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
                flush(stdout)
            end
        end
    end
 
    #=
    Base.Threads.@threads for s in 1:tg.num_nodes
        _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
            flush(stdout)
        end
    end
    =#
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness,time()-start_time
end


function _sstp_accumulate_prova!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},s::Int64,bigint::Bool,temporal_betweenness_centrality::Vector{Float64},bfs_ds)

    
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
    temporal_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.dist) - 1)
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
            temporal_betweenness_centrality[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w]
        end
    end

    return nothing
end

function _sstp_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},s::Int64,bigint::Bool,temporal_betweenness_centrality::Vector{Float64})

    if (bigint)
        bfs_ds = BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_DS(tg.num_nodes, length(keys(tn_index)))
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
    temporal_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.dist) - 1)
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
            temporal_betweenness_centrality[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w]
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_DS(0, 0)
    else
        bfs_ds = BFS_DS(0, 0)
    end
    return nothing
end