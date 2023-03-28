using DataStructures


function rtb_prefix_foremost(tg::temporal_graph, sample_size::Int64,verbose_step::Int64; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Tuple{Float64,Float64,Float64}}
    start_time = time()
    sample = test_sample
    if (length(sample) == 0 || length(sample) != sample_size)
        sample::Array{Int64} = rtb_sample(tg, sample_size)
    end
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    bfs_ds = BFS_prefix_foremost_betweenness(tg.num_nodes)
    
    temporal_prefix_foremost_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    processed_so_far::Int64 = 0
    exec_time::Array{Float64} = zeros(sample_size)
    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i]
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
        exec_time[i] = time() - exec_time[i]
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("RTB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
        end
    end
    for k in 1:lastindex(temporal_prefix_foremost_betweenness_centrality)
        temporal_prefix_foremost_betweenness_centrality[k] = temporal_prefix_foremost_betweenness_centrality[k]/(tg.num_nodes-1)
    end
    return temporal_prefix_foremost_betweenness_centrality, (mean(exec_time), std(exec_time), time() - start_time)

end



#------- Progressive ------------


function progressive_rtb_prefix_foremost(tg::temporal_graph, c::Int64, verbose_step::Int64)::Tuple{Array{Float64},Int64,Tuple{Float64,Float64,Float64}}
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    bfs_ds = BFS_prefix_foremost_betweenness(tg.num_nodes)

    temporal_prefix_foremost_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    processed_so_far::Int64 = 0
    repeat::Bool = true
    k::Int64 = 0
    to_sub::Int64 = 0
    max_betweenness::Float64 = 0.0
    exec_time::Array{Float64} = Array{Float64}([])
    while repeat
        push!(exec_time,time())
        s = sample(1:tg.num_nodes)
        k+=1
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

                if temporal_prefix_foremost_betweenness_centrality[v] >= c * tg.num_nodes
                    repeat = false
                end
                if temporal_prefix_foremost_betweenness_centrality[v] > max_betweenness
                    max_betweenness = temporal_prefix_foremost_betweenness_centrality[v]
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
    for y in 1:lastindex(temporal_prefix_foremost_betweenness_centrality)
        temporal_prefix_foremost_betweenness_centrality[y] = temporal_prefix_foremost_betweenness_centrality[y] *(1/(k * (tg.num_nodes-1)))
    end
    return temporal_prefix_foremost_betweenness_centrality,k, (mean(exec_time), std(exec_time), time() - start_time)

end




