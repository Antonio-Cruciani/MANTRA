using DataStructures

function rtb_sample(tg::temporal_graph, sample_size::Int64)::Array{Int64}
    sample_nodes::Array{Int64} = []
    s::Int64 = 0
    while length(sample_nodes) < sample_size
        s = sample(1:tg.num_nodes)
        push!(sample_nodes, s)
    end
    return sample_nodes
end

function rtb(tg::temporal_graph, sample_size::Int64,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Tuple{Float64,Float64,Float64}}
    start_time = time()
    sample = test_sample
    if (length(sample) == 0 || length(sample) != sample_size)
        sample::Array{Int64} = rtb_sample(tg, sample_size)
    end
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    if (bigint)
        bfs_ds = BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_DS(tg.num_nodes, length(keys(tn_index)))
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
    exec_time::Array{Float64} = zeros(sample_size)
    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i]
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
        processed_so_far = processed_so_far + 1
        exec_time[i] = time() - exec_time[i]

        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((tg.num_nodes*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("RTB. Processed " * string(processed_so_far) * "/" * string(sample_size) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
       
    end
    for k in 1:lastindex(temporal_betweenness_centrality)
        temporal_betweenness_centrality[k] = temporal_betweenness_centrality[k]/(tg.num_nodes-1)
    end
    return temporal_betweenness_centrality, (mean(exec_time), std(exec_time), time() - start_time)
end





function sum_all_results(tg::temporal_graph,temporal_betweenness_centrality::Array{Float64},sample_size::Int64)
    apx::Array{Float64} = [0.0 for i in 1:tg.num_nodes]

    for l in 1:sample_size
        for u in 1:tg.num_nodes
            apx[u] += temporal_betweenness_centrality[(l-1)*tg.num_nodes+u]
        end
    end
    return apx
end




#------ Progressive 

function progressive_rtb(tg::temporal_graph, c::Int64, verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Int64,Tuple{Float64,Float64,Float64}}
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    if (bigint)
        bfs_ds = BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_DS(tg.num_nodes, length(keys(tn_index)))
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
    repeat::Bool = true
    k::Int64 = 0
    max_betweenness::Float64 = 0.0
    exec_time::Array{Float64} = Array{Float64}([])
    value::Float64 = 0
    while repeat
        push!(exec_time,time())
        s =  sample(1:tg.num_nodes)
        k+=1
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
                value = (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w]
                bfs_ds.delta_sh[tni_v] += value
                temporal_betweenness_centrality[v] += value
                if temporal_betweenness_centrality[v] >= c * tg.num_nodes
                    repeat = false
                end
                if temporal_betweenness_centrality[v] > max_betweenness
                    max_betweenness = temporal_betweenness_centrality[v]
                end
            end
        end
        processed_so_far = processed_so_far + 1
        exec_time[k] = time() - exec_time[k]

        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            if repeat
                println("P-RTB. Processed " * string(processed_so_far) *" nodes in " * finish_partial * " seconds | max betwenness ",round(max_betweenness; digits=4), " threshold ",c*tg.num_nodes, " | sampling new nodes...")
                flush(stdout)
            else
                println("P-RTB. Processed " * string(processed_so_far) *" nodes in " * finish_partial * " seconds | max betwenness ",round(max_betweenness; digits=4), " threshold " ,c*tg.num_nodes," | converged.")
                flush(stdout)
            end
        end
    end
    for y in 1:lastindex(temporal_betweenness_centrality)
        temporal_betweenness_centrality[y] = temporal_betweenness_centrality[y] *(1/(k*(tg.num_nodes-1)))
    end
    return temporal_betweenness_centrality,k, (mean(exec_time), std(exec_time), time() - start_time)

end

