using DataStructures
struct BFS_prefix_foremost_betweenness
    t_min::Array{Int64}
    sigma::Array{UInt128}
    delta::Array{Float64}
    predecessors::Array{Set{Int64}}
    # priority queue Q and Stack S
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    stack::Stack{Int64}
    function BFS_prefix_foremost_betweenness(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{UInt128}(undef,nn),Array{Float64}(undef,nn),Array{Set{Int64}}(undef,nn), PriorityQueue{Tuple{Int64,Int64,Int64},Int64}(),Stack{Int64}() )
    end
end


struct BFS_prefix_foremost_betweenness_BI
    t_min::Array{Int64}
    sigma::Array{BigInt}
    delta::Array{Float64}
    predecessors::Array{Set{Int64}}
    # priority queue Q and Stack S
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    stack::Stack{Int64}
    function BFS_prefix_foremost_betweenness_BI(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{BigInt}(undef,nn),Array{Float64}(undef,nn),Array{Set{Int64}}(undef,nn), PriorityQueue{Tuple{Int64,Int64,Int64},Int64}(),Stack{Int64}() )
    end
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



