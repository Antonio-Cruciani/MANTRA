


function threaded_progressive_onbra_prefix_foremost_set(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,group::Array{Int64},theta::Float64,verbose_step::Int64)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]

    local_B::Vector{Dict{Float64,Float64}} =  [Dict{Float64,Float64}() for i in 1:nthreads()]
    local_B_2::Vector{Vector{Float64}} =  [zeros(tg.num_nodes) for i in 1:nthreads()]
    C::Dict{Int64,Int16} = Dict{Int64,Int16}([])
    for i in 1:lastindex(group)
        C[group[i]] = 0
    end
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
            _p_onbra_pfm_accumulate_set!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()],local_B[Base.Threads.threadid()],local_B_2[Base.Threads.threadid()],C)
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial )
            end
        end
        B_vectorized = reduce_dictionary(local_B)
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized,sample_size_schedule[j])
            delta_i = delta/2^k
            xi = 1/theta* (2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j])))
        else
            xi = Inf
        end
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-ONBRA-GROUP-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        
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


function _p_onbra_pfm_accumulate_set!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,B_1::Vector{Float64},B::Dict{Float64,Float64},B_2::Vector{Float64},C::Dict{Int64,Int16})
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
                if haskey(C,v)
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
                end
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
    return nothing


end



function threaded_progressive_onbra_prefix_foremost_topk_dep(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,k::Int64,verbose_step::Int64)
    start_time::Float64 = time()
    delta_1::Float64 = 1/11
    delta_2::Float64 = (10/11 -1)/(10*(1/11-1))
    println("Delta "*string((1-delta_1)*(1-delta_2)))

    betweenness::Array{Float64} = zeros(tg.num_nodes)
    group::Array{Int64} = Array{Int64}([])
    ranking::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
    CI::Array{Tuple{Int64,Float64,Float64}} = Array{Tuple{Int64,Float64,Float64}}([])

    final_ranking::Array{Tuple{Int64,Float64,Float64,Float64}} = Array{Tuple{Int64,Float64,Float64,Float64}}([])

    xi_1::Float64 = 0.0
    xi_2::Float64 = 0.0
    y_1::Float64 = 0.0
    y_2::Float64 = 0.0
    b_k::Float64 = 0.0
    b_k_2::Float64 = 0.0
    z::Float64 = 0.0

    betweenness,_,xi_1,_ = threaded_progressive_onbra_prefix_foremost(tg,initial_sample,epsilon,delta_1,geo,0)

    for u in 1:tg.num_nodes
        push!(CI,(u,betweenness[u]-xi_1,betweenness[u]+xi_1))
    end
    sort!(CI, by=CI->-CI[2])
    println(CI)
    println("V_K "*string(CI[k][2]))
    for u in 1:tg.num_nodes
        if CI[u][3] >= CI[k][2]
            push!(final_ranking,(u,betweenness[u],CI[u][2],CI[u][3]))
        end
    end
    println(final_ranking)
    #=
    sort!(ranking, by=ranking->-ranking[2])
    b_k = ranking[k][2]
    y_1 = b_k -xi_1
    println("y' = "*string(y_1))
    println("First phase completed in "*string(round(time() - start_time; digits=4)))
    for i in 1:lastindex(ranking)
        if ranking[i][2] >= b_k - 2*xi_1
            push!(group,ranking[i][1])
        end
    end
    println("Group retreived | size = "*string(lastindex(group)))
    betweenness,_,xi_2,_ = threaded_progressive_onbra_prefix_foremost_set(tg,initial_sample,epsilon,delta_2,geo,group,y_1,0)
    println("Second phase completed in "*string(round(time() - start_time; digits=4)))
    ranking = Array{Tuple{Int64,Float64}}([])
    for u in 1:tg.num_nodes
        push!(ranking,(u,betweenness[u]))
    end
    sort!(ranking, by=ranking->-ranking[2])

    b_k_2 = ranking[k][2]
    y_2 = b_k_2/(1+xi_2)
    z = max(y_1,y_2)

    for i in 1:lastindex(ranking)
        if ranking[i][2] >= z*(1-xi_2)
            push!(final_ranking,(ranking[i][1],ranking[i][2]))
        end
    end
    =#
    return final_ranking,time()-start_time
end



# Bernstein



function threaded_progressive_onbra_prefix_foremost_bernstein_set(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,group::Array{Int64},theta::Float64,verbose_step::Int64)::Tuple{Array{Float64},Array{Int64},Float64,Float64}


    start_time = time()
    print_algorithm_status("ONBRA","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for i in 1:nthreads()]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    betweenness::Vector{Float64} = zeros(tg.num_nodes)
    C::Dict{Int64,Int16} = Dict{Int64,Int16}([])
    for i in 1:lastindex(group)
        C[group[i]] = 0
    end
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
            _p_onbra_pfm_bernstein_accumulate_set!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()],C)
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            end
        end
       
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        xi = 1/theta * theoretical_error_bound(reduced_betweenness,reduce(+,t_bc),sample_size_schedule[j],delta/2^k)


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-ONBRA. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
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


function _p_onbra_pfm_bernstein_accumulate_set!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},s::Int64,z::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64},C::Dict{Int64,Int16})
    bigint::Bool = false
    bfs_ds::BFS_ONBRA_PFM_DS = BFS_ONBRA_PFM_DS(tg.num_nodes)
    summand,b ,b_1 = def_summand(bigint)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)

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

                if haskey(C,v)
                    B_1[p][v] += summand
                    B_2[v] += summand
                end
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
    return nothing


end


function threaded_progressive_onbra_prefix_foremost_bernstein_topk(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,k::Int64,verbose_step::Int64)
    start_time::Float64 = time()
    delta_1::Float64 = 1/11
    delta_2::Float64 = (10/11 -1)/(10*(1/11-1))
    println("Delta "*string((1-delta_1)*(1-delta_2)))

    betweenness::Array{Float64} = zeros(tg.num_nodes)
    group::Array{Int64} = Array{Int64}([])
    ranking::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
    final_ranking::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])

    xi_1::Float64 = 0.0
    xi_2::Float64 = 0.0
    y_1::Float64 = 0.0
    y_2::Float64 = 0.0
    b_k::Float64 = 0.0
    b_k_2::Float64 = 0.0
    z::Float64 = 0.0
    betweenness,_,xi_1,_ = threaded_progressive_onbra_prefix_foremost_bernstein(tg,initial_sample,epsilon,delta_1,geo,0)

    for u in 1:tg.num_nodes
        push!(ranking,(u,betweenness[u]))
    end
    sort!(ranking, by=ranking->-ranking[2])
    b_k = ranking[k][2]
    y_1 = b_k -xi_1
    println("b_k ",b_k)
    for i in 1:lastindex(ranking)
        if ranking[i][2] >= b_k - 2*xi_1
            push!(group,ranking[i][1])
        end
    end
    betweenness,_,xi_2,_ = threaded_progressive_onbra_prefix_foremost_bernstein_set(tg,initial_sample,epsilon,delta_2,geo,group,y_1,0)
    ranking = Array{Tuple{Int64,Float64}}([])
    for u in 1:tg.num_nodes
        push!(ranking,(u,betweenness[u]))
    end
    sort!(ranking, by=ranking->-ranking[2])

    b_k_2 = ranking[k][2]
    y_2 = b_k_2/(1+xi_2)
    z = max(y_1,y_2)

    for i in 1:lastindex(ranking)
        if ranking[i][2] >= z*(1-xi_2)
            push!(final_ranking,(ranking[i][1],ranking[i][2]))
        end
    end
    return final_ranking,start_time-time()
end





# = newone

function check_stopping_condition_topk(top_k::Array{Tuple{Int64,Float64,Float64,Float64}},epsilon::Float64)::Bool
    converged::Bool = true
    for i in 1:lastindex(top_k)
        converged = converged & ((top_k[i][2]/(1+epsilon) <= top_k[i][3])) & ((top_k[i][2]/(1-epsilon) >= top_k[i][4]))
    end
    return converged
end

function threaded_progressive_onbra_prefix_foremost_topk(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,tk::Int64,verbose_step::Int64)

    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    approx_top_k::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
    top_k::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
    b_est::Float64 = 0.0
    lb_v_k::Float64 = 0.0
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
        approx_top_k =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
        top_k =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
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
            _p_onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()],local_B[Base.Threads.threadid()],local_B_2[Base.Threads.threadid()])
            if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
                finish_partial = string(round(time() - start_time; digits=4))
                println("P-ONBRA-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial )
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
        println("P-ONBRA-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        betweenness = reduce(+, local_temporal_betweenness)
        for u in 1:tg.num_nodes
            b_est = betweenness[u]/sample_size_schedule[j]
            push!(approx_top_k,(u,b_est,b_est-xi,b_est+xi))
        end
        sort!(approx_top_k, by=approx_top_k->-approx_top_k[3])
        lb_v_k = approx_top_k[tk][3]
        println(lb_v_k)
        for i in 1:lastindex(approx_top_k)
            if approx_top_k[i][4] >= lb_v_k
                push!(top_k,approx_top_k[i])
            end
        end
        if check_stopping_condition_topk(top_k,epsilon)
            keep_sampling = false
        else
            j+=1
        end
        println("**************************++")

    end
    
    #betweenness = reduce(+, local_temporal_betweenness)
    #betweenness = betweenness .* [1/sample_size_schedule[j]]
    return top_k,sample_size_schedule,xi,time()-start_time
end