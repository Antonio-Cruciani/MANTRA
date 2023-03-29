using DataStructures

struct BFS_DS
    sigma::Array{UInt128}
    dist::Array{Int64}
    sigma_t::Array{UInt128}
    dist_t::Array{Int64}
    delta_sh::Array{Float64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    queue::Queue{Tuple{Int64,Int64}}
    stack::Stack{Tuple{Int64,Int64}}
    function BFS_DS(nn::Int64, ntn::Int64)
        return new(Array{UInt128}(undef, nn), Array{Int64}(undef, nn), Array{UInt128}(undef, ntn), Array{Int64}(undef, ntn), Array{Float64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), Queue{Tuple{Int64,Int64}}(), Stack{Tuple{Int64,Int64}}())
    end
end

struct BI_BFS_DS
    sigma::Array{BigInt}
    dist::Array{Int64}
    sigma_t::Array{BigInt}
    dist_t::Array{Int64}
    delta_sh::Array{Float64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    queue::Queue{Tuple{Int64,Int64}}
    stack::Stack{Tuple{Int64,Int64}}
    function BI_BFS_DS(nn::Int64, ntn::Int64)
        return new(Array{BigInt}(undef, nn), Array{Int64}(undef, nn), Array{BigInt}(undef, ntn), Array{Int64}(undef, ntn), Array{Float64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), Queue{Tuple{Int64,Int64}}(), Stack{Tuple{Int64,Int64}}())
    end
end




function next_temporal_neighbors(tal::Array{Array{Tuple{Int64,Int64}}}, u::Int64, t::Int64)::Array{Tuple{Int64,Int64}}
    neig::Array{Tuple{Int64,Int64}} = tal[u]
    left::Int64 = 1
    right::Int64 = length(neig)
    pos::Int64 = length(neig) + 1
    mid::Int64 = -1
    while (left <= right)
        mid = (left + right) ÷ 2
        if (neig[mid][2] <= t)
            left = mid + 1
        else
            pos = mid
            right = mid - 1
        end
    end
    return neig[pos:end]
end

function temporal_node_index(tg::temporal_graph)::Dict{Tuple{Int64,Int64},Int64}
    d::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    current_index = 1
    for edge in tg.temporal_edges
        if (get(d, (edge[2], edge[3]), 0) == 0)
            d[(edge[2], edge[3])] = current_index
            current_index = current_index + 1
        end
    end
    for s in 1:tg.num_nodes
        d[(s, 0)] = current_index
        current_index = current_index + 1
    end
    return d
end

function temporal_shortest_betweenness(tg::temporal_graph, verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
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
    for s in 1:tg.num_nodes
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
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("TSB. Processed " * string(processed_so_far) * "/" * string(tg.num_nodes) * " nodes in " * finish_partial * " seconds")
            flush(stdout)
        end
    end
    finish_total::Float64 = round(time() - start_time; digits=4)
    return temporal_betweenness_centrality, finish_total
end



