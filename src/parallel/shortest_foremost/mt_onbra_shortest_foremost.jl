function threaded_onbra_shortest_foremost(tg::temporal_graph,sample_size::Int64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}

    start_time = time()
    sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    z::Int64 = 0
    println("Using ",nthreads()," Trheads")

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i][1]
        z = sample[i][2]
        _onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("ONBRA-SFM. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size]
    return betweenness,time()-start_time
end


function _onbra_sfm_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,temporal_betweenness_centrality::Vector{Float64})
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.t_hit[u] = -1
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    t_z_min = Inf
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min && t < t_z_min
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
                            if t_z_min > t_w
                                t_z_min = t_w
                            end
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, neig)
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], temporal_node)
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
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
    end
    if bfs_ds.dist[z] > 0
        for tn in 1:lastindex(bfs_ds.dist_t)
            bfs_ds.sigma_z[tn] = 0
            bfs_ds.boolean_matrix[tn] = false
        end
        tni = tn_index[(s, 0)]
        bfs_ds.sigma_z[tni] = 1
    end
    for t in 1:lastindex(tg.file_time)
        tni = get(tn_index, (z, t), 0)
        if tni > 0 && bfs_ds.sigma_t[tni] > 0
            for pred in bfs_ds.predecessors[tni]
                tni_w = tn_index[(pred[1], pred[2])]
                if (!bigint && bfs_ds.sigma_z[tni_w] == typemax(UInt128))
                    println("Overflow occurred with sample (", s, ",", z, ")")
                    return [], 0.0
                end
                bfs_ds.sigma_z[tni_w] += 1
                if !bfs_ds.boolean_matrix[tni_w]
                    enqueue!(bfs_ds.backward_queue, pred)
                    bfs_ds.boolean_matrix[tni_w] = true
                end
            end
        end
    end
    while length(bfs_ds.backward_queue) > 0
        temporal_node = dequeue!(bfs_ds.backward_queue)
        tni = tn_index[(temporal_node[1], temporal_node[2])]
        if temporal_node[1] != s
            temporal_betweenness_centrality[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
            for pred in bfs_ds.predecessors[tni]
                tni_w = tn_index[(pred[1], pred[2])]
                if (!bigint && bfs_ds.sigma_z[tni_w] > typemax(UInt128) - bfs_ds.sigma_z[tni])
                    println("Overflow occurred with sample (", s, ",", z, ")")
                    return [], 0.0
                end
                bfs_ds.sigma_z[tni_w] += bfs_ds.sigma_z[tni]
                if !bfs_ds.boolean_matrix[tni_w]
                    enqueue!(bfs_ds.backward_queue, pred)
                    bfs_ds.boolean_matrix[tni_w] = true
                end
            end
        end
    end
    return nothing

end

function threaded_progressive_onbra_shortest_foremost(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]

    local_B::Vector{Dict{Float64,Float64}} =  [Dict{Float64,Float64}() for i in 1:nthreads()]
    local_B_2::Vector{Vector{Float64}} =  [zeros(tg.num_nodes) for i in 1:nthreads()]
    B_vectorized::Array{Float64} =Array{Float64}([])
    s::Int64 = 0
    z::Int64 = 0
    sample_size_schedule::Array{Int64} = [0,initial_sample]
    xi::Float64 = 0
    sampled_so_far::Int64 = 0
    new_sample::Int64 = 0
    keep_sampling::Bool = true
    k::Int64 = 0
    j::Int64 = 2
    finish_partial::String = ""
    println("Using ",nthreads()," Trheads")
    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end
        Base.Threads.@threads for i in 1:(sample_size_schedule[j]-sample_size_schedule[j-1])
            sampled_so_far+=1
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            _p_onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()],local_B[Base.Threads.threadid()],local_B_2[Base.Threads.threadid()])
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial)
            end
        end
        B_vectorized = reduce_dictionary(local_B)
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized,sample_size_schedule[j])
            delta_i = delta/2^k
            xi = 2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j]))
        else
            xi = Inf
        end
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-ONBRA-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)

        
        if xi <= epsilon
            keep_sampling = false
        else
            j+=1
        end

    end
    
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size_schedule[j]]
    return betweenness,sample_size_schedule,xi,time()-start_time
end



function _p_onbra_sfm_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,B_1::Vector{Float64},B::Dict{Float64,Float64},B_2::Vector{Float64})
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    end
    summand,b ,b_1 = def_summand(bigint)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.t_hit[u] = -1
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    t_z_min = Inf
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min && t < t_z_min
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
                            if t_z_min > t_w
                                t_z_min = t_w
                            end
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, neig)
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], temporal_node)
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
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
    end
    if bfs_ds.dist[z] > 0
        for tn in 1:lastindex(bfs_ds.dist_t)
            bfs_ds.sigma_z[tn] = 0
            bfs_ds.boolean_matrix[tn] = false
        end
        tni = tn_index[(s, 0)]
        bfs_ds.sigma_z[tni] = 1
    
        for t in 1:lastindex(tg.file_time)
            tni = get(tn_index, (z, t), 0)
            if tni > 0 && bfs_ds.sigma_t[tni] > 0
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] == typemax(UInt128))
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += 1
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
        while length(bfs_ds.backward_queue) > 0
            temporal_node = dequeue!(bfs_ds.backward_queue)
            tni = tn_index[(temporal_node[1], temporal_node[2])]
            if temporal_node[1] != s
                if bigint
                    summand = Float64(bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z])  , RoundUp)
                else
                    summand = (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                end
                 # Updating phase
                 b = B_2[temporal_node[1]]
                 b_1 = b + summand^2
                 if !haskey(B,b_1) 
                     B[b_1] = 1
                 else
                     B[b_1] += 1
                 end
                 if b > 0 && B[b] >= 1
                     B[b] -= 1
                 end
                 if b > 0 && B[b] == 0
                     delete!(B, b)
                 end
                 B_1[temporal_node[1]] += summand
                 B_2[temporal_node[1]] += summand^2
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] > typemax(UInt128) - bfs_ds.sigma_z[tni])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += bfs_ds.sigma_z[tni]
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
    end
    return nothing
end



#------------------------------------------------------
# Progressive TRK using Bernstein Bound to compute ξ
#------------------------------------------------------


function threaded_progressive_onbra_shortest_foremost_bernstein(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    print_algorithm_status("ONBRA","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)

    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for i in 1:nthreads()]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    betweenness::Vector{Float64} = zeros(tg.num_nodes)

    s::Int64 = 0
    z::Int64 = 0
    sample_size_schedule::Array{Int64} = [0,initial_sample]
    xi::Float64 = 0
    sampled_so_far::Int64 = 0
    new_sample::Int64 = 0
    keep_sampling::Bool = true
    k::Int64 = 0
    j::Int64 = 2
    finish_partial::String = ""
    println("Using ",nthreads()," Trheads")
    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end

        Base.Threads.@threads for i in 1:(sample_size_schedule[j]-sample_size_schedule[j-1])
            sampled_so_far+=1
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            _p_onbra_sfm_bernstein_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            end
        end
       
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)

        xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc),sample_size_schedule[j],delta/2^k)


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-ONBRA. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        if xi <= epsilon
            keep_sampling = false
        else
            j+=1
            local_temporal_betweenness = [[zeros(tg.num_nodes)] for i in 1:nthreads()]
        end

    end
    
    #_reduce_list_of_arrays!(local_temporal_betweenness,betweenness,sample_size_schedule[j],tg.num_nodes)
    betweenness = reduce(+,t_bc) .* [1/sample_size_schedule[j]]
    return betweenness,sample_size_schedule,xi,time()-start_time
end



function _p_onbra_sfm_bernstein_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    end
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)
    summand,b ,b_1 = def_summand(bigint)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    for u in 1:tg.num_nodes
        bfs_ds.t_hit[u] = -1
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:lastindex(bfs_ds.dist_t)
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0))
    t_z_min = Inf
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = tn_index[(u, t)]
        if bfs_ds.dist_t[tni] < d_z_min && t < t_z_min
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
                            if t_z_min > t_w
                                t_z_min = t_w
                            end
                        end
                    end
                    enqueue!(bfs_ds.forward_queue, neig)
                end
                if bfs_ds.dist_t[tni_w] == bfs_ds.dist_t[tni] + 1
                    if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma_t[tni_w])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_t[tni_w] += bfs_ds.sigma_t[tni]
                    push!(bfs_ds.predecessors[tni_w], temporal_node)
                    if bfs_ds.dist_t[tni_w] == bfs_ds.dist[w]
                        if (!bigint && bfs_ds.sigma_t[tni] > typemax(UInt128) - bfs_ds.sigma[w])
                            println("Overflow occurred with sample (", s, ",", z, ")")
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
    end
    if bfs_ds.dist[z] > 0
        for tn in 1:lastindex(bfs_ds.dist_t)
            bfs_ds.sigma_z[tn] = 0
            bfs_ds.boolean_matrix[tn] = false
        end
        tni = tn_index[(s, 0)]
        bfs_ds.sigma_z[tni] = 1
    
        for t in 1:lastindex(tg.file_time)
            tni = get(tn_index, (z, t), 0)
            if tni > 0 && bfs_ds.sigma_t[tni] > 0
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] == typemax(UInt128))
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += 1
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
        while length(bfs_ds.backward_queue) > 0
            temporal_node = dequeue!(bfs_ds.backward_queue)
            tni = tn_index[(temporal_node[1], temporal_node[2])]
            if temporal_node[1] != s
                if bigint
                    summand = Float64(bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z])  , RoundUp)
                else
                    summand = (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
                end
                B_1[p][temporal_node[1]] += summand
                B_2[temporal_node[1]] += summand
                for pred in bfs_ds.predecessors[tni]
                    tni_w = tn_index[(pred[1], pred[2])]
                    if (!bigint && bfs_ds.sigma_z[tni_w] > typemax(UInt128) - bfs_ds.sigma_z[tni])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[tni_w] += bfs_ds.sigma_z[tni]
                    if !bfs_ds.boolean_matrix[tni_w]
                        enqueue!(bfs_ds.backward_queue, pred)
                        bfs_ds.boolean_matrix[tni_w] = true
                    end
                end
            end
        end
    end
    return nothing
end