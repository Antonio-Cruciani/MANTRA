function threaded_rtb(tg::temporal_graph,sample_size::Int64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    sample::Array{Int64} = rtb_sample(tg, sample_size)
    println("Using ",nthreads()," Trheads")

 

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i]
        _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / i )-(time() - start_time) ; digits=4))
            println("RTB-SH. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/(sample_size*(tg.num_nodes-1))]
    return betweenness,time()-start_time
end



#------------------------------------------------------
# Progressive RTB using Bernstein Bound to compute ξ
#------------------------------------------------------


function threaded_progressive_rtb_bernstein(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    print_algorithm_status("RTB","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)

    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for i in 1:nthreads()]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    betweenness::Vector{Float64} = zeros(tg.num_nodes)
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]

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
            s = sample(1:tg.num_nodes)
            _sstp_accumulate_bernstein!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
        end
        sampled_so_far+=sample_size_schedule[j]-sample_size_schedule[j-1]

        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        xi = theoretical_error_bound(reduced_betweenness .* [1/(tg.num_nodes-1)],reduce(+,t_bc) .* [1/(tg.num_nodes-1)] ,sample_size_schedule[j],delta/2^k)


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-RTB. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        if xi <= epsilon
            keep_sampling = false
        else
            j+=1
            local_temporal_betweenness = [[zeros(tg.num_nodes)] for i in 1:nthreads()]
        end

    end
    
    _reduce_list_of_arrays!(local_temporal_betweenness,betweenness,sample_size_schedule[j],tg.num_nodes)
    betweenness = betweenness .* [1/(sample_size_schedule[j]*(tg.num_nodes-1))]
    return betweenness,sample_size_schedule,xi,time()-start_time
end



function _sstp_accumulate_bernstein!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},s::Int64,bigint::Bool,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})

    if (bigint)
        bfs_ds = BI_BFS_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_DS(tg.num_nodes, length(keys(tn_index)))
    end
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)

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
    sub = (count(x -> x >= 0, bfs_ds.dist) - 1) 
    B_1[p][s] -= sub * (1/(tg.num_nodes-1))
    B_2[s] -= sub
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
            B_1[p][v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w] * 1/(tg.num_nodes-1)
            B_2[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_sh[tni_w] 
        end
    end

    return nothing
end