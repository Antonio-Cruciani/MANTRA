function threaded_progressive_det_era_prefix_foremost(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64,algo::String="ob")::Tuple{Array{Float64},Array{Int64},Float64,Float64}
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"

    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

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
    omega_max::Int64 = 1000
    t_diam::Float64 = 0.0
    diam::Float64 = -1.0
    if (diam == -1) 
        println("Approximating diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". VD = "*string(diam))
        flush(stdout)
    end
    
    omega_max = trunc(Int,(0.5/epsilon^2) * ((floor(log2(diam-2)))+log(1/delta)))
    println("ω = ",omega_max)
    
    println("Maximum sample size "*string(omega_max))
    println("Using ",nthreads()," Trheads")
    flush(stdout)
    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end
        Base.Threads.@threads for i in 1:(sample_size_schedule[j]-sample_size_schedule[j-1])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "ob"
                _p_onbra_pfm_det_era_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()],local_B[Base.Threads.threadid()],local_B_2[Base.Threads.threadid()])
            elseif algo == "trk"
                println("TODO")
            elseif algo == "rtb"
                println("TODO")
            end
                if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial )
            end
        end
        sampled_so_far+=sample_size_schedule[j]-sample_size_schedule[j-1]

        B_vectorized = reduce_dictionary(local_B)
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized,sample_size_schedule[j])
            delta_i = delta/2^k
            xi = 2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j]))
        else
            xi = Inf
        end
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-ONBRA-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        
        if xi <= epsilon || sampled_so_far >= omega_max
            keep_sampling = false
        else
            j+=1
        end

    end
    
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size_schedule[j]]
    return betweenness,sample_size_schedule,xi,time()-start_time
end


function _p_onbra_pfm_det_era_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,B_1::Vector{Float64},B::Dict{Float64,Float64},B_2::Vector{Float64})
    bigint::Bool = false
    bfs_ds::BFS_ONBRA_PFM_DS = BFS_ONBRA_PFM_DS(tg.num_nodes)
    summand,b ,b_1 = def_summand(bigint)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)

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
        for v in 1:lastindex(bfs_ds.sigma)
            bfs_ds.sigma_z[v] = 0
            bfs_ds.boolean_matrix[v] = false
        end
        
        #bfs_ds.sigma_z[s] = 1
        
        for pred in bfs_ds.predecessors[z]
            bfs_ds.sigma_z[pred] = 1
            if !bfs_ds.boolean_matrix[pred]
                enqueue!(bfs_ds.backward_queue,pred)
                bfs_ds.boolean_matrix[pred] = true
            end
        end
        while length(bfs_ds.backward_queue) > 0
            v = dequeue!(bfs_ds.backward_queue)
            if v != s
                if bigint
                    summand = Float64(bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z])   , RoundUp)
                else
                    summand = (bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z]))
                end

                # Updating phase
                b = B_2[v]
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
                B_1[v] += summand
                B_2[v] += summand^2
                for pred in bfs_ds.predecessors[v]
                    if (!bigint && bfs_ds.sigma_z[pred] > typemax(UInt128) - bfs_ds.sigma_z[pred])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[pred] +=  bfs_ds.sigma_z[v]
                    if !bfs_ds.boolean_matrix[pred]
                        enqueue!(bfs_ds.backward_queue, pred) 
                        bfs_ds.boolean_matrix[pred] = true
                    end
                end
            end
        end               
    end
    bfs_ds = nothing

    return nothing


end