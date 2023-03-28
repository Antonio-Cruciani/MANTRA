using DataStructures

struct BFS_SFM_DS
    t_hit::Array{Int64}
    sigma::Array{UInt128}
    dist::Array{Int64}
    sigma_t::Array{UInt128}
    dist_t::Array{Int64}
    delta_fm::Array{Float64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    queue::Queue{Tuple{Int64,Int64}}
    stack::Stack{Tuple{Int64,Int64}}
    function BFS_SFM_DS(nn::Int64, ntn::Int64)
        return new(Array{Int64}(undef,nn),Array{UInt128}(undef, nn), Array{Int64}(undef, nn), Array{UInt128}(undef, ntn), Array{Int64}(undef, ntn), Array{Float64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), Queue{Tuple{Int64,Int64}}(), Stack{Tuple{Int64,Int64}}())
    end
end

struct BI_BFS_SFM_DS
    t_hit::Array{Int64}
    sigma::Array{BigInt}
    dist::Array{Int64}
    sigma_t::Array{BigInt}
    dist_t::Array{Int64}
    delta_fm::Array{Float64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    queue::Queue{Tuple{Int64,Int64}}
    stack::Stack{Tuple{Int64,Int64}}
    function BI_BFS_SFM_DS(nn::Int64, ntn::Int64)
        return new(Array{Int64}(undef,nn),Array{BigInt}(undef, nn), Array{Int64}(undef, nn), Array{BigInt}(undef, ntn), Array{Int64}(undef, ntn), Array{Float64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), Queue{Tuple{Int64,Int64}}(), Stack{Tuple{Int64,Int64}}())
    end
end




function temporal_shortest_foremost_betweenness(tg::temporal_graph, verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    if (bigint)
        bfs_ds = BI_BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    end
    temporal_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    v::Int64 = -1
    t_v::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    processed_so_far::Int64 = 0
    for s in 1:tg.num_nodes
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
        temporal_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.dist) - 1)
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
                temporal_betweenness_centrality[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
            end
        end
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("TSFB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
            flush(stdout)
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_betweenness_centrality, finish_total
end


function temporal_shortest_foremost_betweenness_multithread(tg::temporal_graph, verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    
    temporal_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
   
    processed_so_far::Int64 = 0
    lk = ReentrantLock()
    Threads.@threads for s in 1:tg.num_nodes
        if (bigint)
            bfs_ds = BI_BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
        else
            bfs_ds = BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
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
        try lock(lk) 
            temporal_betweenness_centrality[s] -= (count(x -> x >= 0, bfs_ds.dist) - 1)
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
                    temporal_betweenness_centrality[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
                end
            end
        finally
            unlock(lk)
        end
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TSFB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_betweenness_centrality, finish_total
end
