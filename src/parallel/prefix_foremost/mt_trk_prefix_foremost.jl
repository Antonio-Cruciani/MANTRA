function threaded_trk_prefix_foremost(tg::temporal_graph,sample_size::Int64,verbose_step::Int64)::Tuple{Array{Float64},Float64}

    start_time = time()
    sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    z::Int64 = 0
    println("Using ",nthreads()," Trheads")

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i][1]
        z = sample[i][2]
        _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TRK-PFM. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size]
    return betweenness,time()-start_time
end


function _trk_pfm_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,temporal_betweenness_centrality::Vector{Float64})
    bigint::Bool  = false
    bfs_ds = BFS_PFM_SRTP_DS(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    cur_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
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
    t_z_min::Float64 = Inf
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
                        temporal_betweenness_centrality[pred] += 1
                    end
                    break
                end
            end

        end    
    end
    return nothing
end










function threaded_progressive_trk_prefix_foremost(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,algo::String="trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". VD = "*string(diam))
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
        println("ω = ",omega)
    else
        omega = trunc(Int,(1.0/(2*eps^2))*log2(2*tg.num_nodes/delta))
    end
    println("Maximum sample size "*string(omega))
    tau::Int64 = trunc(Int64,omega/start_factor)
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap phase "*string(tau)*" iterations")
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        end

    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    eps_lb::Array{Float64} = zeros(tg.num_nodes)
    eps_ub::Array{Float64} = zeros(tg.num_nodes)
    delta_lb_min_guess::Array{Float64} = [0.0]
    delta_ub_min_guess::Array{Float64} = [0.0]
    delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
    delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
    _compute_δ_guess!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess) 
    println("Bootstrap completed ")
    local_temporal_betweenness = [zeros(tg.num_nodes) for i in 1:nthreads()]

    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    while sampled_so_far < omega && !stop[1]
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            end
        end
        sampled_so_far += sample_step
        betweenness = reduce(+, local_temporal_betweenness)
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-TRK-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
    end
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end








#------------------------------------------------------
# Progressive TRK using Bernstein Bound to compute ξ
#------------------------------------------------------

function threaded_progressive_trk_prefix_foremost_bernstein(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    print_algorithm_status("TRK","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

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
            _p_trk_pfm_bernstein_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-TRK. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            end
        end
       
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc),sample_size_schedule[j],delta/2^k)


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-TRK. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        if xi <= epsilon
            keep_sampling = false
        else
            j+=1
            local_temporal_betweenness = [[zeros(tg.num_nodes)] for i in 1:nthreads()]
        end

    end
    
    _reduce_list_of_arrays!(local_temporal_betweenness,betweenness,sample_size_schedule[j],tg.num_nodes)
    betweenness = betweenness .* [1/sample_size_schedule[j]]
    return betweenness,sample_size_schedule,xi,time()-start_time
end



function _p_trk_pfm_bernstein_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    bigint::Bool  = false
    bfs_ds = BFS_PFM_SRTP_DS(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    cur_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
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
    t_z_min::Float64 = Inf
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
                        B_1[p][pred] += 1
                        B_2[pred] += 1
                    end
                    break
                end
            end

        end    
    end
    return nothing
end
