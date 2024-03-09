function threaded_rtb_prefix_foremost(tg::temporal_graph,sample_size::Int64,verbose_step::Int64)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    sample::Array{Int64} = rtb_sample(tg, sample_size)
    println("Using ",nthreads()," Trheads")

 

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i]
        _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("RTB-PFM. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/(sample_size*(tg.num_nodes-1))]
    return betweenness,time()-start_time
end


function threaded_progressive_rtb_prefix_foremost_bernstein(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    print_algorithm_status("RTB","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for i in 1:nthreads()]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    betweenness::Vector{Float64} = zeros(tg.num_nodes)

    s::Int64 = 0
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
            s = sample(1:tg.num_nodes)
                      
            _ssptp_accumulate_bernstein!(tg,tal,s,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-RTB-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            end
        end
       
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        xi = theoretical_error_bound(reduced_betweenness .* [1/(tg.num_nodes-1)],reduce(+,t_bc) .* [1/(tg.num_nodes-1)] ,sample_size_schedule[j],delta/2^k)


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-RTB-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
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


function _ssptp_accumulate_bernstein!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    bfs_ds::BFS_prefix_foremost_betweenness = BFS_prefix_foremost_betweenness(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    sub::Float64 = 0.0
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)
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
    sub  = (count(x -> x >= 0, bfs_ds.t_min) - 1)
    B_1[p][s] -= sub * (1/(tg.num_nodes-1))
    B_2[s] -= sub
    while length(bfs_ds.stack) != 0
        w = pop!(bfs_ds.stack)
        for v in bfs_ds.predecessors[w]
            summand = (Float64(bfs_ds.sigma[v] / bfs_ds.sigma[w])) * bfs_ds.delta[w]
            bfs_ds.delta[v] += summand
            B_1[p][v] += summand * (1/(tg.num_nodes-1))
            B_2[v] += summand
        end
    end
    bfs_ds = nothing

    return nothing

end