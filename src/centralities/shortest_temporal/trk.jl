using DataStructures
using StatsBase

struct BFS_SRTP_DS
    sigma::Array{UInt128}
    dist::Array{Int64}
    sigma_t::Array{UInt128}
    sigma_z::Array{UInt128}
    dist_t::Array{Int64}
    predecessors::Array{Set{Tuple{Int64,Int64,Int64}}}
    boolean_array::Array{Bool}
    forward_queue::Queue{Tuple{Int64,Int64,Int64}}
    function BFS_SRTP_DS(nn::Int64, ntn::Int64)
        return new(Array{UInt128}(undef, nn), Array{Int64}(undef, nn), Array{UInt128}(undef, ntn), zeros(Int64, ntn), Array{Int64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn+1), falses(ntn), Queue{Tuple{Int64,Int64,Int64}}())
    end
end

struct BI_BFS_SRTP_DS
    sigma::Array{BigInt}
    dist::Array{Int64}
    sigma_t::Array{BigInt}
    sigma_z::Array{BigInt}
    dist_t::Array{Int64}
    predecessors::Array{Set{Tuple{Int64,Int64,Int64}}}
    boolean_array::Array{Bool}
    forward_queue::Queue{Tuple{Int64,Int64,Int64}}
    function BI_BFS_SRTP_DS(nn::Int64, ntn::Int64)
        return new(Array{BigInt}(undef, nn), Array{Int64}(undef, nn), Array{BigInt}(undef, ntn), zeros(Int64, ntn), Array{Int64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn+1), falses(ntn), Queue{Tuple{Int64,Int64,Int64}}())
    end
end


function trk_sample(tg::temporal_graph, sample_size::Int64)::Array{Tuple{Int64,Int64}}
    sample_pairs::Array{Tuple{Int64,Int64}} = []
    s::Int64 = 0
    z::Int64 = 0
    while length(sample_pairs) < sample_size
        s, z = sample(1:tg.num_nodes, 2, replace=false)
        push!(sample_pairs, (s, z))
    end
    return sample_pairs
end


function empirical_variance(tilde_b::Array{Float64}, sample_size::Int64, v::Int64)::Float64
    n::Int64 = div(length(tilde_b), sample_size)
    variance::Float64 = 0
    for i in 1:sample_size
        for j in (i+1):sample_size
            variance += (((tilde_b[(i-1)*n+v] - tilde_b[(j-1)*n+v]) / sample_size)^2)
        end
    end
    return variance / (sample_size * (sample_size - 1))
end

function initialize_weights(bigint::Bool)
    if bigint
        totalWeight_bi::BigInt = 0
        randomEdge_bi::BigInt = 0
        curWeight_bi::BigInt = 0
        curEdge_bi::BigInt = 0
        return(totalWeight_bi,randomEdge_bi,curWeight_bi,curEdge_bi)
    else
        totalWeight::UInt128 = 0
        randomEdge::UInt128 = 0
        curWeight::UInt128 = 0
        curEdge::UInt128 = 0
        return(totalWeight,randomEdge,curWeight,curEdge)
    end
end
function trk(tg::temporal_graph, sample_size::Int64,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Int64,Float64}
    start_time = time()
    sample = test_sample
    if (length(sample) == 0 || length(sample) != sample_size)
        sample::Array{Tuple{Int64,Int64}} = trk_sample(tg, sample_size)
    end

    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SRTP_DS(tg.num_nodes, indexes)
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64,Int64} = (-1, -1,-1)
    processed_so_far::Int64 = 0
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
    cur_w::Tuple{Int64,Int64} = (-1,-1)
    temporal_betweenness_centrality::Array{Float64} = zeros(tg.num_nodes)
    exec_time::Array{Float64} = zeros(sample_size)
    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i][1]
        z = sample[i][2]
        t_z::Int64 = tg.temporal_edges[lastindex(tg.temporal_edges)][3]+1
        tn_index[(z,t_z)] = indexes+1
        for u in 1:tg.num_nodes
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
        
        bfs_ds.predecessors[indexes+1] = Set{Tuple{Int64,Int64,Int64}}()
        enqueue!(bfs_ds.forward_queue, (s, 0,tni))
        d_z_min = Inf
        while length(bfs_ds.forward_queue) != 0
            temporal_node = dequeue!(bfs_ds.forward_queue)
            u = temporal_node[1]
            t = temporal_node[2]
            #tni = tn_index[(u, t)]
            tni = temporal_node[3]
            if bfs_ds.dist_t[tni] < d_z_min
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
                            push!(bfs_ds.predecessors[indexes+1], (neig[1],neig[2],tni_w))
                        end
                    end
                end
            end
        end
        if bfs_ds.dist[z] > 0
           
            totalWeight = 0
            randomEdge = 0
            curWeight = 0
            totalWeight = bfs_ds.sigma[z]
            cur_w = (z,t_z)
            tni = indexes+1
            #println("PATH  s = ",s, " z = ",z)
            while cur_w[1] != s
                #tni = get(tn_index, cur_w, 0)
                if cur_w == (z,t_z)
                    totalWeight = bfs_ds.sigma[z]
                else
                    totalWeight = bfs_ds.sigma_t[tni]
                end
                #println(" TOTAL WE ",totalWeight)
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
        exec_time[i] = time() - exec_time[i]
        delete!(tn_index,(z,t_z))
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("RTK. Processed " * string(processed_so_far) * "/" * string(sample_size) * " pairs in " * finish_partial * " seconds")
        end
       

    end
   
    return temporal_betweenness_centrality,  time() - start_time

end

function average_results(temporal_betweenness_centrality::Array{Float64},tg::temporal_graph,num_samples::Int64)
    apx::Array{Float64} = [0.0 for i in 1:tg.num_nodes]
    for l in 1:num_samples
        for u in 1:tg.num_nodes
            apx[u] += temporal_betweenness_centrality[(l-1)*tg.num_nodes+u]
        end
    end
    for u in 1:tg.num_nodes
        apx[u] = apx[u]/num_samples
    end 
    return apx
end


