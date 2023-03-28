struct BFS_PFM_SRTP_DS
    t_hit::Array{Int64}
    sigma::Array{UInt128}
    predecessors::Array{Set{Int64}}
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    function BFS_PFM_SRTP_DS(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{UInt128}(undef,nn),Array{Set{Int64}}(undef,nn), PriorityQueue{Tuple{Int64,Int64,Int64},Int64}() )
    end
end
struct BFS_PFM_SRTP_DS_BI
    t_hit::Array{Int64}
    sigma::Array{BigInt}
    predecessors::Array{Set{Int64}}
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    function BFS_PFM_SRTP_DS_BI(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{BigInt}(undef,nn),Array{Set{Int64}}(undef,nn), PriorityQueue{Tuple{Int64,Int64,Int64},Int64}() )
    end
end

function trk_prefix_foremost(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Int64,Tuple{Float64,Float64,Float64}}

    start_time = time()
   
    
    sample::Array{Tuple{Int64,Int64}} = trk_sample(tg, sample_size)
    
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    if (bigint)
        bfs_ds = BFS_PFM_SRTP_DS_BI(tg.num_nodes)
    else
        bfs_ds = BFS_PFM_SRTP_DS(tg.num_nodes)
    end
    tilde_b::Array{Float64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    cur_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    processed_so_far::Int64 = 0
    exec_time::Array{Float64} = zeros(sample_size)
    avg_path_length::Float64 = 0
    path_sampled::Int64 = 0
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)

    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i][1]
        z = sample[i][2]
        for u in 1:tg.num_nodes
            bfs_ds.t_hit[u] = -1
            bfs_ds.sigma[u] = 0
            bfs_ds.predecessors[u] = Set{Int64}()
        end
        
        bfs_ds.sigma[s] = 1
        bfs_ds.t_hit[s] = 0
        for tn in tal[s]
            enqueue!(bfs_ds.priority_queue,(s,tn[1],tn[2]),tn[2])
        end
        t_z_min = Inf
        while length(bfs_ds.priority_queue) != 0
            temporal_edge = dequeue!(bfs_ds.priority_queue)
            v = temporal_edge[1]
            w = temporal_edge[2]
            t_w = temporal_edge[3]
            if t_w < t_z_min
                if bfs_ds.t_hit[w] == -1
                    bfs_ds.t_hit[w] = t_w
                    for neig in next_temporal_neighbors(tal,w,t_w)
                        enqueue!(bfs_ds.priority_queue,(w,neig[1],neig[2]),neig[2])
                    end
                    if z == w &&  t_z_min > t_w
                        t_z_min = t_w
                    end
                end
                if bfs_ds.t_hit[w] == t_w
                    if (!bigint && bfs_ds.sigma[v] > typemax(UInt128) - bfs_ds.sigma[w])
                        println("Overflow occurred with source ", s)
                        return [], 0.0
                    end
                    bfs_ds.sigma[w] += bfs_ds.sigma[v]
                    push!(bfs_ds.predecessors[w], v)
                end
            end
        end
        if bfs_ds.sigma[z] > 0
            cur_w = z
            while cur_w != s
                totalWeight = 0
                if cur_w == z
                    totalWeight = bfs_ds.sigma[z]
                else
                    for pred in bfs_ds.predecessors[cur_w]
                        totalWeight+=bfs_ds.sigma[pred]
                    end
                end
                randomEdge = rand(0:totalWeight-1)
                curEdge = 0
                for pred in bfs_ds.predecessors[cur_w]
                    curEdge+= bfs_ds.sigma[pred]
                    if curEdge > randomEdge
                        cur_w = pred
                        if pred[1] != s && pred[1] != z
                            tilde_b[pred] += 1
                        end
                        break
                    end
                end

            end

            path_sampled += 1
        
        end
        exec_time[i] = time() - exec_time[i]
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("TRK. Processed " * string(processed_so_far) * "/" * string(sample_size) * " pairs in " * finish_partial * " seconds")
            flush(stdout)
        end
    end
    avg_path_length = 0

    return tilde_b,path_sampled ,(mean(exec_time), std(exec_time), time() - start_time)

end

