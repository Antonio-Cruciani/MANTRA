function temporal_node_index_srtp(tg::temporal_graph)::Dict{Tuple{Int64,Int64},Int64}
    d::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    t_z::Int64 = tg.temporal_edges[lastindex(tg.temporal_edges)][3]+1
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
    for z in 1:tg.num_nodes
        d[(z, t_z)] = current_index
        current_index = current_index + 1
    end
    return d
end


function threaded_trk(tg::temporal_graph,sample_size::Int64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}

    start_time = time()
    sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    z::Int64 = 0
    println("Using ",nthreads()," Trheads")

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i][1]
        z = sample[i][2]
        _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / processed_so_far )-(time() - start_time) ; digits=4))
            println("TRK-SH. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/sample_size]
    return betweenness,time()-start_time
end


function _trk_sh_accumulate!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,temporal_betweenness_centrality::Vector{Float64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SRTP_DS(tg.num_nodes, indexes)
    end
    t_z::Int64 = tg.temporal_edges[lastindex(tg.temporal_edges)][3]+1
    index_z::Int64 = tn_index[(z,t_z)]
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64,Int64} = (-1, -1,-1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
    cur_w::Tuple{Int64,Int64} = (-1,-1)
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
    enqueue!(bfs_ds.forward_queue, (s, 0,tni))
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
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
                        push!(bfs_ds.predecessors[index_z], (neig[1],neig[2],tni_w))
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
        tni = index_z
        while cur_w[1] != s
            if cur_w == (z,t_z)
                totalWeight = bfs_ds.sigma[z]
            else
                totalWeight = bfs_ds.sigma_t[tni]
            end
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
    return nothing

end






function threaded_progressive_trk(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
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
            _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
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
    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    while sampled_so_far < omega && !stop[1]
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            end
        end
        sampled_so_far += sample_step
        betweenness = reduce(+, local_temporal_betweenness)
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-TRK-SH. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
    end
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end
function commpute_f(btilde::Float64, iter_num::Int64,δ_l::Float64,ω::Int64)::Float64
    tmp::Float64 = ω/iter_num - 1.0/3.0
    err_chern::Float64 = (log(1.0/δ_l))*1.0/iter_num*(-tmp+sqrt(tmp*tmp +2*btilde * ω/(log(1.0/δ_l))))
    return min(err_chern,btilde)
end

function compute_g(btilde::Float64, iter_num::Int64,δ_u::Float64,ω::Int64)::Float64
    tmp::Float64 = ω/iter_num + 1.0/3.0
    err_chern::Float64 = (log(1.0/δ_u))*1.0/iter_num*(tmp+sqrt(tmp*tmp +2*btilde * ω/(log(1.0/δ_u))))
    return min(err_chern,1-btilde)
end
function _compute_finished!(stop::Array{Bool},omega::Int64,betweenness::Array{Float64},sampled_so_far::Int64,eps::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},delta_lb_guess::Array{Float64},delta_ub_guess::Array{Float64},delta_lb_min_guess::Float64,delta_ub_min_guess::Float64)
    #j::Int64 = 1
    n::Int64 = lastindex(betweenness)
    all_finished::Bool = true
    finished::Array{Bool} = falses(n)
    Base.Threads.@threads for i in 1:(n-1)
        betweenness[i] = betweenness[i] / sampled_so_far
        eps_lb[i] = commpute_f(betweenness[i],sampled_so_far,delta_lb_guess[i],omega)
        eps_ub[i] = compute_g(betweenness[i],sampled_so_far,delta_ub_guess[i],omega)
        #j+=1
    end
    betweenness[n] = betweenness[n] / sampled_so_far
    eps_lb[n] = commpute_f(betweenness[n],sampled_so_far,delta_lb_min_guess,omega)
    eps_ub[n] = compute_g(betweenness[n],sampled_so_far,delta_ub_min_guess,omega)
    for i in 1:n
        finished[i] = (eps_lb[i] < eps) & (eps_ub[i] < eps)
        all_finished = all_finished & finished[i]
    end
    stop[1] = all_finished
    return nothing
end


function compute_bet_err(eps::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},start_factor::Int64)::Tuple{Array{Float64},Array{Float64}}
    Base.Threads.@threads for i in 1:n
        eps_lb[i] = eps
        eps_ub[i] = eps
    end
    return eps_lb,eps_ub
end





function _compute_δ_guess!(betweenness::Array{Float64},eps::Float64,delta::Float64,balancing_factor::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},delta_lb_min_guess::Array{Float64},delta_ub_min_guess::Array{Float64},delta_lb_guess::Array{Float64},delta_ub_guess::Array{Float64}) 

    n::Int64 = lastindex(betweenness)
    a::Float64 = 0
    b::Float64 = 1.0 / eps / eps* log(n* 4* (1-balancing_factor)/delta)
    c::Float64 = (a+b)/2
    summation::Float64 = 0.0
    
    for i in 1:n
        eps_lb[i] = eps
        eps_ub[i] = eps
    end
    while (b-a > eps/10.0)
        c = (b+a)/2
        summation = 0
        for i in 1:n
            summation += exp(-c * eps_lb[i]*eps_lb[i] / betweenness[i] )
            summation += exp(-c * eps_ub[i]*eps_ub[i] / betweenness[i] )
        end
        summation += exp(-c * eps_lb[n-1]*eps_lb[n-1] / betweenness[n-1] ) * (n-n)
        summation += exp(-c * eps_ub[n-1]*eps_ub[n-1] / betweenness[n-1] ) * (n-n)
        if (summation >= delta/2.0 * (1-balancing_factor))
            a = c 
        else
            b = c
        end
    end
    delta_lb_min_guess[1] = exp(-b * eps_lb[n-1]* eps_lb[n-1] / betweenness[n-1]) + delta*balancing_factor/4.0 / n
    delta_ub_min_guess[1] = exp(-b * eps_ub[n-1]* eps_ub[n-1] / betweenness[n-1] ) + delta*balancing_factor/4.0 / n
    Base.Threads.@threads for i in 1:n
        delta_lb_guess[i] = exp(-b *  eps_lb[i]*eps_lb[i]/ betweenness[i])+ delta*balancing_factor/4.0 / n
        delta_ub_guess[i] = exp(-b *  eps_ub[i]*eps_ub[i] / betweenness[i]) + delta*balancing_factor/4.0 / n
    end
    return nothing
end




#------------------------------------------------------
# Progressive TRK using Bernstein Bound to compute ξ
#------------------------------------------------------

function threaded_progressive_trk_bernstein(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Array{Int64},Float64,Float64}

    start_time = time()
    print_algorithm_status("TRK","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)

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
            _p_trk_sh_accumulate_bernstein!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
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




function _p_trk_sh_accumulate_bernstein!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SRTP_DS(tg.num_nodes, indexes)
    end
    push!(B_1,zeros(tg.num_nodes))
    p::Int64 = lastindex(B_1)
    t_z::Int64 = tg.temporal_edges[lastindex(tg.temporal_edges)][3]+1
    index_z::Int64 = tn_index[(z,t_z)]
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64,Int64} = (-1, -1,-1)
    totalWeight,randomEdge,curWeight,curEdge = initialize_weights(bigint)
    cur_w::Tuple{Int64,Int64} = (-1,-1)
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
    enqueue!(bfs_ds.forward_queue, (s, 0,tni))
    d_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
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
                        push!(bfs_ds.predecessors[index_z], (neig[1],neig[2],tni_w))
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
        tni = index_z
        while cur_w[1] != s
            if cur_w == (z,t_z)
                totalWeight = bfs_ds.sigma[z]
            else
                totalWeight = bfs_ds.sigma_t[tni]
            end
            randomEdge = rand(0:totalWeight-1)
            curEdge = 0
            for pred in bfs_ds.predecessors[tni]
                pred_i = pred[3]
                curEdge += bfs_ds.sigma_t[pred_i]
                cur_w = (pred[1],pred[2])
                tni = pred_i
                if curEdge > randomEdge
                    if pred[1]!= s && pred[1] != z
                        B_1[p][pred[1]] += 1
                        B_2[pred[1]] += 1
                    end
                    break
                end
            end
            
        end
       
        
    end
    return nothing

end

