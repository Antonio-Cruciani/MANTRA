function threaded_progressive_wub_shortest_foremost(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,algo::String="trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end

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
            _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "rtb"
            _ssftp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
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
                _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "rtb"
                _ssftp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
            end        
        end
        sampled_so_far += sample_step
        betweenness = reduce(+, local_temporal_betweenness)
        if algo == "rtb"
            betweenness =  betweenness.*[1/(tg.num_nodes-1)]
        end
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
    end

    
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time
    


end



function threaded_progressive_wub_shortest_foremost_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,verbose_step::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    union_sample::Int64 = min(tg.num_nodes,max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/nthreads()),k+20))
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_foremost_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
    else
        omega = trunc(Int,(1.0/(2*eps^2))*log2(2*tg.num_nodes/delta))
    end
    println("Top-k algorithm: k =  "*string(k)*"  union sample = "*string(union_sample))
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
            _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "rtb"
            _ssftp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
    betweenness = betweenness .* [1/tau]
    for u in 1:tg.num_nodes
        if algo == "rtb"
            push!(approx_top_k,(u,betweenness[u]/(tg.num_nodes-1)))
        else
            push!(approx_top_k,(u,betweenness[u]))
        end
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    eps_lb::Array{Float64} = zeros(tg.num_nodes)
    eps_ub::Array{Float64} = zeros(tg.num_nodes)
    delta_lb_min_guess::Array{Float64} = [0.0]
    delta_ub_min_guess::Array{Float64} = [0.0]
    delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
    delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
    _compute_δ_guess_topk!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess,k,approx_top_k,start_factor,union_sample) 
    println("Bootstrap completed ")
    local_temporal_betweenness = [zeros(tg.num_nodes) for i in 1:nthreads()]
    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    while sampled_so_far < omega && !stop[1]
        approx_top_k = Array{Tuple{Int64,Float64}}([])
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_sfm_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "rtb"
                _ssftp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
            end
        end
        sampled_so_far += sample_step
        betweenness = reduce(+, local_temporal_betweenness)

        for u in 1:tg.num_nodes
            if algo == "rtb"
                push!(approx_top_k,(u,betweenness[u]/(tg.num_nodes-1)))
            else
                push!(approx_top_k,(u,betweenness[u]))
            end
        end
        sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
        _compute_finished_topk!(stop,omega,approx_top_k[begin:union_sample],sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1],union_sample)   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-SFM (TOP-K). Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
    end
    for i in 1:lastindex(approx_top_k)
        approx_top_k[i] = (approx_top_k[i][1],approx_top_k[i][2]/sampled_so_far)
    end
    return approx_top_k,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end