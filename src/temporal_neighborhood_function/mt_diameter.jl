using DataStructures
using Graphs
struct BFS_DIAM_SH
    dist::Array{Int64}
    dist_t::Array{Int64}
    queue::Queue{Tuple{Int64,Int64}}
    function BFS_DIAM_SH(nn::Int64, ntn::Int64)
        return new(Array{Int64}(undef, nn), Array{Int64}(undef, ntn),  Queue{Tuple{Int64,Int64}}())
    end
end

function threaded_temporal_shortest_diameter(tg::temporal_graph,sample_size::Int64,verbose_step::Int64,threshold::Float64 = 0.9)::Tuple{Int64,Float64,Float64,Float64,Float64,Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    local_temporal_distance_distribution::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    temporal_hop_table::Array{Float64} = Array{Float64}([])
    processed_so_far::Int64 = 0
    diameter::Int64 = 0
    accum::Int64 = 0
    h::Int64 = 1
    s::Int64 = 0
    sample_space::Array{Int64} = Array{Int64}([])
    println("Using ",nthreads()," Trheads")
    if sample_size == 0
        println("Computing exact distance measures")
        sample_size = tg.num_nodes
        sample_space = [i for i in 1:tg.num_nodes]
    else
        println("Computing approximation using ",sample_size," seeds")
        sample_space = rtb_sample(tg, sample_size)
    end

    Base.Threads.@threads for i in 1:sample_size
         s = sample_space[i]
        _sstp_diameter!(tg,tal,tn_index,s,local_temporal_distance_distribution[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TS-Diam. Processed " * string(processed_so_far) * "/" * string(sample_size) * " nodes in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    dd = reduce(+, local_temporal_distance_distribution)
    diameter = get_diameter(dd)
    temporal_hop_table = zeros(diameter+1)
    for h in 1:(diameter+1)
        accum += dd[h]
        temporal_hop_table[h] = tg.num_nodes * accum/sample_size
    end
    avg_dist::Float64 = average_distance(temporal_hop_table)
    #println("Average distance ",avg_dist)
    eff_diam::Float64 = effective_diameter(temporal_hop_table,threshold)
    total_couples::Float64 = total_reachable_couples(temporal_hop_table)
    alpha::Float64 = total_couples / (tg.num_nodes*(tg.num_nodes-1))
    #println("Diameter ",diameter, " Effective Diameter ",eff_diam," Average Distance ",avg_dist, " #Couples ",total_couples, " α ",alpha)
    return diameter,avg_dist,eff_diam,total_couples,alpha,time()-start_time
end
function get_diameter(dd::Array{Float64})::Int64
    diam::Int64 = 0
    for i in 1:lastindex(dd)
        if dd[i] != 0
            diam = i
        else
            return diam -1
        end
    end
    return diam -1
end
function average_distance(ht::Array{Float64})::Float64
    if length(ht) == 0
        return 0.0
    end
    distance::Array{Float64} = distance_function(ht)
    m::Float64 = 0.0
    for i in 1:lastindex(distance)
        m += (distance[i] * i)
    end
    return m/ht[length(ht)]
end

function distance_function(ht::Array{Float64})::Array{Float64}
    table::Array{Float64} = copy(ht)
    for i in lastindex(table):-1:2
        table[i] -= table[i-1]
    end
    return table
end

function effective_diameter(h::Array{Float64},threshold::Float64)::Float64
    couples::Float64 = total_reachable_couples(h)
    d::Int64 = 1
    while (h[d]/couples) < threshold
        d+=1
    end
    return (d-1)+ interpolate(h[d-1],h[d],threshold*couples)
end

function total_reachable_couples(h::Array{Float64})::Float64
    return h[lastindex(h)]
end

function interpolate(y0::Float64,y1::Float64,y::Float64)::Float64
    return (y-y0)/(y1-y0)
end


function reduce_max(local_diameters::Vector{Int64})::Int64
    return maximum(local_diameters)
end

function _sstp_diameter!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},s::Int64,temporal_distance_distribution::Vector{Float64})

    bfs_ds = BFS_DIAM_SH(tg.num_nodes, length(keys(tn_index)))
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    u::Int64 = 0
    t::Int64 = 0
    tni::Int64 = 0
    w::Int64 = 0
    t_w::Int64 = 0
    tni_w::Int64 = 0
    for u in 1:tg.num_nodes
        bfs_ds.dist[u] = -1
    end
    for tn in 1:(length(tn_index))
        bfs_ds.dist_t[tn] = -1
    end
    tni = tn_index[(s, 0)]
    bfs_ds.dist_t[tni] = 0
    bfs_ds.dist[s] = 0 
    temporal_distance_distribution[bfs_ds.dist[s]+1] += 1
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
                    temporal_distance_distribution[bfs_ds.dist[w]+1] += 1
                end
                enqueue!(bfs_ds.queue, neig)
            end
        end
    end
    return nothing
end