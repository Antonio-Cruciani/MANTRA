

function threaded_progressive_wub(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,bigint::Bool,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
        flush(stdout)
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+1+log(1/delta)))
    else
        omega = trunc(Int,(1.0/(2*eps^2))*log2(2*tg.num_nodes/delta))
    end
    println("Maximum sample size "*string(omega))
    tau::Int64 = trunc(Int64,omega/start_factor)
    s::Int64 = 0
    z::Int64 = 0
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
            elseif algo == "ob"
                _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
            elseif algo == "rtb"
                _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t])
            end
        end
    end
    #=
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "rtb"
            _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        end
    end
    =#
    betweenness = reduce(+, local_temporal_betweenness)
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
    betweenness = betweenness .* [1/tau]
    eps_lb::Array{Float64} = zeros(tg.num_nodes)
    eps_ub::Array{Float64} = zeros(tg.num_nodes)
    delta_lb_min_guess::Array{Float64} = [0.0]
    delta_ub_min_guess::Array{Float64} = [0.0]
    delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
    delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
    _compute_δ_guess!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess) 
    println("Bootstrap completed ")
    flush(stdout)
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    task_size = cld(ntasks, ntasks)
    vs_active = [i for i in 1:ntasks]
    while sampled_so_far < omega && !stop[1]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:ntasks, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t])
                end
            end
        end
        #=
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "rtb"
                _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
            end
        end
        =#
        sampled_so_far += ntasks
        betweenness = reduce(+, local_temporal_betweenness)
        if algo == "rtb"
            betweenness =  betweenness.*[1/(tg.num_nodes-1)]
        end
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-SH. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            flush(stdout)
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
  
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end








function threaded_progressive_wub_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,verbose_step::Int64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end

    balancing_factor::Float64 = 0.001

    mc_trials::Int64 = 25
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    wv::Array{Float64} = Array{Float64}([])
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]    
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    union_sample::Int64 = min(tg.num_nodes,max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/ntasks),k+20))
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
        flush(stdout)
    end
    omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
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
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    if vc_upper_bund == true
        println("Bootstrap phase "*string(tau)*" iterations")
        println("Bootstrap using VC dimension")
        flush(stdout)
        @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t])
                end
            end
        end
    else
        tau = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
        tau = trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
        task_size = cld(tau, ntasks)
        vs_active = [i for i in 1:tau]
        println("Bootstrap phase "*string(tau)*" iterations")
        println("Bootstrap using Variance")
        flush(stdout)
        @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "ob"
                    _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                elseif algo == "rtb"
                    _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
                end
    
            end
        end
    end
    #=
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "rtb"
            _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        end
    end
    =#
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
    if !vc_upper_bund 
        norm::Float64 = 1.0
        if algo == "rtb"
            norm = 1/(tg.num_nodes-1)
        end
        wv = reduce(+,local_wv)
        sp_lengths = reduce(+,local_sp_lengths) 
        max_tbc = maximum(betweenness)
        max_wv = maximum(wv)/tau
        avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,trunc(Int,diam),sp_lengths,tau,true,norm)
        top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv,delta/8,tau)
        println("AVERAGE DIAM UB "*string(avg_diam_ub))
        max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
        println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
        println("Sup tbc est "*string(max_tbc))
        println("Sup emp wimpy variance "*string(max_wv/tau))
        omega = trunc(Int,max_num_samples)
        flush(stdout)
        local_wv = [[]]
        wv = []
        local_sp_lengths = [[]]
        mcrade = [[]]
    end
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
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    task_size = cld(ntasks, ntasks)
    vs_active = [i for i in 1:ntasks]
    while sampled_so_far < omega && !stop[1]
        approx_top_k = Array{Tuple{Int64,Float64}}([])
        @sync for (t, task_range) in enumerate(Iterators.partition(1:ntasks, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[t])
                end
            end
        end
        #=
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_sh_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "rtb"
                _sstp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
            end
        end
        =#
        sampled_so_far += ntasks
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
            println("P-WUB-"*algo*"-SH (TOP-K). Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            flush(stdout)
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
    for i in 1:lastindex(approx_top_k)
        approx_top_k[i] = (approx_top_k[i][1],approx_top_k[i][2]/sampled_so_far)
    end
    return approx_top_k,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time

end





function compute_bet_err_topk(eps::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},start_factor::Int64,k::Int64,approx_top_k::Array{Tuple{Int64,Float64}},union_sample::Int64)::Tuple{Array{Float64},Array{Float64}}
    n::Int64 = lastindex(eps_lb)
    bet::Array{Float64} = zeros(n)
    max_error::Float64 = sqrt(start_factor) * eps/4
    Base.Threads.@threads for i in 1:union_sample 
        bet[i] = approx_top_k[i][2]
    end
    eps_ub[1] = max(eps,(bet[1]-bet[2])/2)
    eps_lb[1] = 10
    Base.Threads.@threads for i in 2:k
        eps_lb[i] = max(eps,(bet[i-1]-bet[i])/2)
        eps_ub[i] = max(eps,(bet[i]-bet[i+1])/2)
    end
    Base.Threads.@threads for i in (k+1):union_sample
        eps_lb[i] = 10
        eps_ub[i] = max(eps,bet[k-1]+(bet[k-1]-bet[k])/2 - bet[i])
    end
    for i in 1:(k-1)
        if bet[i] - bet[i + 1] < max_error
            eps_lb[i] = eps
            eps_ub[i] = eps
            eps_lb[i+1] = eps
            eps_ub[i+1] = eps
        end
    end
    for i in (k+1):union_sample
        if bet[k] - bet[i] < max_error
            eps_lb[k] = eps
            eps_ub[k] = eps
            eps_lb[i] = eps
            eps_ub[i] = eps
        end
    end
    return eps_lb,eps_ub
end

function _compute_δ_guess_topk!(betweenness::Array{Float64},eps::Float64,delta::Float64,balancing_factor::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},delta_lb_min_guess::Array{Float64},delta_ub_min_guess::Array{Float64},delta_lb_guess::Array{Float64},delta_ub_guess::Array{Float64},k::Int64,approx_top_k::Array{Tuple{Int64,Float64}},start_factor::Int64,union_sample::Int64) 

    n::Int64 = lastindex(betweenness)
    v::Int64 = -1
    a::Float64 = 0
    b::Float64 = 1.0 / eps / eps* log(n* 4* (1-balancing_factor)/delta)
    c::Float64 = (a+b)/2
    summation::Float64 = 0.0
    eps_lb,eps_ub = compute_bet_err_topk(eps,eps_lb,eps_ub,start_factor,k,approx_top_k,union_sample)
  
    while (b-a > eps/10.0)
        c = (b+a)/2
        summation = 0
        for i in 1:n
            summation += exp(-c * eps_lb[i]*eps_lb[i] / betweenness[i] )
            summation += exp(-c * eps_ub[i]*eps_ub[i] / betweenness[i] )
        end
        summation += exp(-c * eps_lb[union_sample-1]*eps_lb[union_sample-1] / betweenness[union_sample-1] ) * (n-union_sample)
        summation += exp(-c * eps_ub[union_sample-1]*eps_ub[union_sample-1] / betweenness[union_sample-1] ) * (n-union_sample)
        if (summation >= delta/2.0 * (1-balancing_factor))
            a = c 
        else
            b = c
        end
    end
    delta_lb_min_guess[1] = exp(-b * eps_lb[union_sample-1]* eps_lb[union_sample-1] / betweenness[union_sample-1]) + delta*balancing_factor/4.0 / n
    delta_ub_min_guess[1] = exp(-b * eps_ub[union_sample-1]* eps_ub[union_sample-1] / betweenness[union_sample-1] ) + delta*balancing_factor/4.0 / n
    Base.Threads.@threads for v in 1:n
        delta_lb_guess[v] = delta_lb_min_guess[1]
        delta_ub_guess[v] =  delta_ub_min_guess[1] 
    end

    Base.Threads.@threads for i in 1:union_sample
        v = approx_top_k[i][1]
        delta_lb_guess[v] = exp(-b *  eps_lb[i]*eps_lb[i]/ betweenness[i])+ delta*balancing_factor/4.0 / n
        delta_ub_guess[v] = exp(-b *  eps_ub[i]*eps_ub[i] / betweenness[i]) + delta*balancing_factor/4.0 / n
    end 

    return nothing
end


function _compute_finished_topk!(stop::Array{Bool},omega::Int64,top_k_approx::Array{Tuple{Int64,Float64}},sampled_so_far::Int64,eps::Float64,eps_lb::Array{Float64},eps_ub::Array{Float64},delta_lb_guess::Array{Float64},delta_ub_guess::Array{Float64},delta_lb_min_guess::Float64,delta_ub_min_guess::Float64,union_sample::Int64)
    #j::Int64 = 1
    k = lastindex(top_k_approx)
    all_finished::Bool = true
    finished::Array{Bool} = falses(union_sample)
    betweenness::Array{Float64} = zeros(union_sample)
    Base.Threads.@threads for i in 1:(union_sample-1)
        betweenness[i] = top_k_approx[i][2] / sampled_so_far
        eps_lb[i] = commpute_f(betweenness[i],sampled_so_far,delta_lb_guess[top_k_approx[i][1]],omega)
        eps_ub[i] = compute_g(betweenness[i],sampled_so_far,delta_ub_guess[top_k_approx[i][1]],omega)
        #j+=1
    end
    betweenness[union_sample] = top_k_approx[union_sample][2] / sampled_so_far
    eps_lb[union_sample] = commpute_f(betweenness[union_sample],sampled_so_far,delta_lb_min_guess,omega)
    eps_ub[union_sample] = compute_g(betweenness[union_sample],sampled_so_far,delta_ub_min_guess,omega)
    for i in 1:union_sample
        if i == 1
            finished[i] = (betweenness[i] - eps_lb[i] > betweenness[i+1] + eps_ub[i+1] )
        elseif i < k
            finished[i] = (betweenness[i-1] - eps_lb[i-1] > betweenness[i] + eps_ub[i] ) & (betweenness[i] - eps_lb[i] > betweenness[i+1] + eps_ub[i+1] )
        else
            finished[i] = (betweenness[k-1] - eps_ub[k-1] > betweenness[i] + eps_ub[i] )
        end
        
        all_finished = all_finished & finished[i] ||( (eps_lb[i] < eps) & (eps_ub[i] < eps))
    end
    stop[1] = all_finished

    return nothing
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
    
    Base.Threads.@threads for i in 1:n
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

