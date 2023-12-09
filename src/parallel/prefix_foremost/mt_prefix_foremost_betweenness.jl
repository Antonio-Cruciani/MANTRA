function distributed_temporal_prefix_foremost_betweenness(tg::temporal_graph,verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    #local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]

    betweenness = @distributed (+) for s in 1:tg.num_nodes
        temp_betweenness = zeros(tg.num_nodes)
        _ssptp_accumulate!(tg,tal,s,temp_betweenness)
        #=
        if (verbose_step > 0 && s % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            time_to_finish = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TPFM. Processed " * string(s) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
            flush(stdout)
        end
        =#
        temp_betweenness
    end

    return betweenness,time()-start_time
end


function threaded_temporal_prefix_foremost_betweenness(tg::temporal_graph,verbose_step::Int64,force_gc::Bool = false)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    processed_so_far::Int64 = 0
    println("Using ",nthreads()," Trheads")
    vs_active = [i for i in 1:tg.num_nodes]
    d, r = divrem(tg.num_nodes, nthreads())
    ntasks = d == 0 ? r : nthreads()
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    task_size = cld(tg.num_nodes, ntasks)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tg.num_nodes, task_size))
        Threads.@spawn for s in @view(vs_active[task_range])
            _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
            if (Sys.free_memory() / Sys.total_memory() < 0.1)
                clean_gc()
                sleep(0.01)
            end
            processed_so_far = processed_so_far + 1
            if (verbose_step > 0 && processed_so_far % verbose_step == 0)
                finish_partial::String = string(round(time() - start_time; digits=4))
                time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
                println("TPFM. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
                flush(stdout)
                if force_gc
                    clean_gc()
                end
            end
        end
    end
    #=
    Base.Threads.@threads for s in 1:tg.num_nodes
        _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TPFM. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
            flush(stdout)
        end
    end
    =#
    betweenness = reduce(+, local_temporal_betweenness)
    return betweenness,time()-start_time
end

function _ssptp_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,temporal_betweenness_centrality::Vector{Float64})
    bfs_ds::BFS_prefix_foremost_betweenness = BFS_prefix_foremost_betweenness(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    for u in 1:tg.num_nodes
        bfs_ds.delta[u] = 1
        bfs_ds.sigma[u] = 0
        bfs_ds.predecessors[u] = Set{Int64}()
        bfs_ds.t_min[u] = -1
    end
    bfs_ds.t_min[s] = 0
    bfs_ds.sigma[s] = 1
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
    temporal_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.t_min) - 1)
    while length(bfs_ds.stack) != 0
        w = pop!(bfs_ds.stack)
        for v in bfs_ds.predecessors[w]
            summand = (Float64(bfs_ds.sigma[v] / bfs_ds.sigma[w])) * bfs_ds.delta[w]
            bfs_ds.delta[v] += summand
            temporal_betweenness_centrality[v] += summand
        end
    end
    bfs_ds = BFS_prefix_foremost_betweenness(0)
    return nothing

end
function temporal_prefix_foremost_betweenness(tg::temporal_graph, verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    bfs_ds::BFS_prefix_foremost_betweenness = BFS_prefix_foremost_betweenness(tg.num_nodes)
    temporal_prefix_foremost_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    processed_so_far::Int64 = 0
    for s in 1:tg.num_nodes
        for u in 1:tg.num_nodes
            bfs_ds.delta[u] = 1
            bfs_ds.sigma[u] = 0
            bfs_ds.predecessors[u] = Set{Int64}()
            bfs_ds.t_min[u] = -1
        end
        bfs_ds.t_min[s] = 0
        bfs_ds.sigma[s] = 1
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
        temporal_prefix_foremost_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.t_min) - 1)
        while length(bfs_ds.stack) != 0
            w = pop!(bfs_ds.stack)
            for v in bfs_ds.predecessors[w]
                summand = (Float64(bfs_ds.sigma[v] / bfs_ds.sigma[w])) * bfs_ds.delta[w]
                bfs_ds.delta[v] += summand
                temporal_prefix_foremost_betweenness_centrality[v] += summand
            end
        end
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("PREFIX. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
            flush(stdout)
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_prefix_foremost_betweenness_centrality, finish_total
end

