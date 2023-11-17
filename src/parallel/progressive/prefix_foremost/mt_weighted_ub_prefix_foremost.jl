
function threaded_progressive_wub_prefix_foremost_dep(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,algo::String="trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    balancing_factor::Float64 = 0.001

    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". VD = "*string(diam))
        flush(stdout)
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
        println("ω = ",omega)
        flush(stdout)
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
                _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
            elseif algo == "ob"
                _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
            elseif algo == "rtb"
                _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
            end
        end
    end
    #=
    Base.Threads.@threads for i in 1:tau
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
        s = sample[1][1]
        z = sample[1][2]
        if algo == "trk"
            _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "ob"
            _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
        elseif algo == "rtb"
            _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[Base.Threads.threadid()])
        end

    end
    =#
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
    flush(stdout)
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]

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
                    _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
                end
            end
        end
        #=
        Base.Threads.@threads for i in 1:sample_step 
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "ob"
                _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[Base.Threads.threadid()])
            elseif algo == "rtb"
                _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[Base.Threads.threadid()])
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
            println("P-WUB-"*algo*"-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            flush(stdout)
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
   
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end

function threaded_progressive_wub_prefix_foremost(tg::temporal_graph,eps::Float64,delta::Float64,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,geo::Float64 = 1.2,start_factor::Int64 = 100)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    balancing_factor::Float64 = 0.001
    omega::Int64 = 1000
    mc_trials::Int64 = 25
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    if (diam == -1) 
        println("Approximating (pfm)-Temporal Diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,0)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        flush(stdout)
    end
    omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-1)))+1+log(2/delta)))
    start_time_bootstrap = time()
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau = trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    tau = trunc(Int64,max(tau,omega/start_factor))
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    if vc_upper_bund == true
        println("Bootstrap using VC dimension")
    else
        println("Bootstrap using Variance")
    end
    flush(stdout)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _pfm_accumulate_trk!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            else
                _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
    betweenness = betweenness .* [1/tau]
    # getting maximum wimpy variance and tbc
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    max_tbc = maximum(betweenness)
    max_wv = maximum(wv)/tau
    if !vc_upper_bund 
        norm::Float64 = 1.0
        if algo == "rtb"
            norm = 1/(tg.num_nodes-1)
        end
        avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,trunc(Int,diam),sp_lengths,tau,true,norm)
        top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv,delta/8,tau)
        println("Average (pfm)-temporal path (upper bound) "*string(avg_diam_ub))
        max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
        println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
        println("Sup tbc estimation "*string(max_tbc))
        println("Sup empirical wimpy variance "*string(max_wv/tau))
        omega = trunc(Int,max_num_samples)
        flush(stdout)
    end
     # freeing memory
     local_wv = [[]]
     wv = []
     local_sp_lengths = [[]]
     mcrade = [[]]
     eps_lb::Array{Float64} = zeros(tg.num_nodes)
     eps_ub::Array{Float64} = zeros(tg.num_nodes)
     delta_lb_min_guess::Array{Float64} = [0.0]
     delta_ub_min_guess::Array{Float64} = [0.0]
     delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
     delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
     _compute_δ_guess!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess) 
     finish_bootstrap = string(round(time() - start_time_bootstrap; digits=4))
     println("Bootstrap completed in "*finish_bootstrap)
     println("Inferring initial sample size for the geometric sampler")
     flush(stdout)
     first_stopping_samples::Float64 = 0.0
     eps_guess::Float64 = 1.0
     first_sample_lower::Float64 = 1/eps *log(2/delta)
     first_sample_upper::Float64 = omega
     sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
     num_samples = 0
     while first_sample_upper - first_sample_lower> 10
         num_samples = (first_sample_upper+first_sample_lower)÷2
         eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples) + log(2/delta)/num_samples/3
         if eps_guess > eps
             first_sample_lower = num_samples
         else
             first_sample_upper = num_samples
         end
     end
     first_stopping_samples = num_samples
     last_stopping_samples = omega
     println("Maximum number of iterations "*string(omega))
     println("Initial sample size "*string(first_stopping_samples))
     if first_stopping_samples >= last_stopping_samples/4
         first_stopping_samples = last_stopping_samples/4
         println("Initial sample size dropped to "*string(first_stopping_samples))
     end
     flush(stdout)
     next_stopping_samples = trunc(Int,first_stopping_samples)
     prev_stopping_samples::Int64 = 0
     local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
     betweenness = zeros(tg.num_nodes)
     sampled_so_far::Int64 = 0
     stop::Array{Bool} = [false]
     sample_stop::Bool = false
     iteration::Int64 = 0
     while !sample_stop && !stop[1]
        sample_i = next_stopping_samples-prev_stopping_samples
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
                end
            end
        end
        sampled_so_far += sample_i
        if sampled_so_far >= omega
            sample_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(sampled_so_far)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end        
        betweenness = reduce(+, local_temporal_betweenness)
        if algo == "rtb"
            betweenness =  betweenness.*[1/(tg.num_nodes-1)]
        end
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1]) 
        if !sample_stop && !stop[1]
            prev_stopping_samples = next_stopping_samples
            next_stopping_samples = trunc(Int,prev_stopping_samples*geo)
            iteration+=1
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Increasing sample size to "*string(next_stopping_samples))
            flush(stdout)
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
 
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time

end



function threaded_progressive_wub_prefix_foremost_2(tg::temporal_graph,eps::Float64,delta::Float64,verbose_step::Int64,algo::String="trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    balancing_factor::Float64 = 0.001
    geo::Float64 = 1.2
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    if (diam == -1) && (!hb)
        println("Approximating diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". VD = "*string(diam))
        flush(stdout)
    end
    if !hb
        omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+log(1/delta)))
        println("ω = ",omega)
        flush(stdout)
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
                _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
            elseif algo == "ob"
                _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
            elseif algo == "rtb"
                _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
            end
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
    flush(stdout)
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]

    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    next_stopping_samples::Int64 = sample_step
    prev_stopping_samples::Int64 = 0
    sample_stop::Bool = false
    iteration::Int64 = 0
    num_samples::Int64 = 0
    while !sample_stop && !stop[1]
        sample_i = next_stopping_samples-prev_stopping_samples
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
                end
            end
        end

        sampled_so_far += sample_i
        num_samples = sampled_so_far
        if num_samples >= omega
            sample_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed sampling : "*string(num_samples)*"/"*string(omega)* " in "*finish_partial)
            flush(stdout)
        end
        betweenness = reduce(+, local_temporal_betweenness)
        if algo == "rtb"
            betweenness =  betweenness.*[1/(tg.num_nodes-1)]
        end
        _compute_finished!(stop,omega,betweenness,sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1])   
        if (verbose_step > 0 && sampled_so_far % verbose_step == 0)
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds ")
            flush(stdout)
        end
        if !sample_stop && !stop[1]
            prev_stopping_samples = next_stopping_samples
            next_stopping_samples = trunc(Int,prev_stopping_samples*geo)
            iteration+=1
        end
    end
    if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
   
    return betweenness,eps_lb,eps_ub,sampled_so_far,omega,time()-start_time


end

function threaded_progressive_wub_prefix_foremost_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,geo::Float64 = 1.2,start_factor::Int64 = 100,force_gc::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb" 
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
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
    if (diam == -1) 
        println("Approximating (pfm)-Temporal Diameter ")
        _,_,_,_,_,diam,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,0)
        println("Task completed in "*string(round(t_diam;digits = 4))*". VD = "*string(diam))
        flush(stdout)
    end
    omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-1)))+1+log(2/delta)))
    start_time_bootstrap = time()
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau = trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    tau = trunc(Int64,max(tau,omega/start_factor))
    println("Top-k algorithm: k =  "*string(k)*"  union sample = "*string(union_sample))
    println("Bootstrap phase "*string(tau)*" iterations")
    flush(stdout)
    task_size = cld(tau, ntasks)
    vs_active = [i for i in 1:tau]
    if vc_upper_bund == true
        println("Bootstrap using VC dimension")
    else
        println("Bootstrap using Variance")
    end
    flush(stdout)
    @sync for (t, task_range) in enumerate(Iterators.partition(1:tau, task_size))
        Threads.@spawn for _ in @view(vs_active[task_range])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _pfm_accumulate_trk!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _pfm_accumulate_onbra!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            else
                _pfm_accumulate_rtb!(tg,tal,s,z,mc_trials,true,local_temporal_betweenness[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/tau]
    if algo == "rtb"
        betweenness =  betweenness.*[1/(tg.num_nodes-1)]
    end
    # getting maximum wimpy variance and tbc
    wv = reduce(+,local_wv)
    sp_lengths = reduce(+,local_sp_lengths) 
    max_tbc = maximum(betweenness)
    max_wv = maximum(wv)/tau
    if !vc_upper_bund 
        norm::Float64 = 1.0
        if algo == "rtb"
            norm = 1/(tg.num_nodes-1)
        end
        avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,trunc(Int,diam),sp_lengths,tau,true,norm)
        top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
        wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv,delta/8,tau)
        println("Average (pfm)-temporal path (upper bound) "*string(avg_diam_ub))
        max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
        println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
        println("Sup tbc estimation "*string(max_tbc))
        println("Sup empirical wimpy variance "*string(max_wv/tau))
        omega = trunc(Int,max_num_samples)
        flush(stdout)
    end
    # freeing memory
    local_wv = [[]]
    wv = []
    local_sp_lengths = [[]]
    mcrade = [[]]
    for u in 1:tg.num_nodes
        push!(approx_top_k,(u,betweenness[u])) 
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    eps_lb::Array{Float64} = zeros(tg.num_nodes)
    eps_ub::Array{Float64} = zeros(tg.num_nodes)
    delta_lb_min_guess::Array{Float64} = [0.0]
    delta_ub_min_guess::Array{Float64} = [0.0]
    delta_lb_guess::Array{Float64} = zeros(tg.num_nodes)
    delta_ub_guess::Array{Float64} = zeros(tg.num_nodes)
    _compute_δ_guess_topk!(betweenness,eps,delta,balancing_factor,eps_lb,eps_ub,delta_lb_min_guess,delta_ub_min_guess,delta_lb_guess,delta_ub_guess,k,approx_top_k,start_factor,union_sample) 
    finish_bootstrap = string(round(time() - start_time_bootstrap; digits=4))
    println("Bootstrap completed in "*finish_bootstrap)
    println("Inferring initial sample size for the geometric sampler")
    flush(stdout)    
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    num_samples = 0
    while first_sample_upper - first_sample_lower> 10
        num_samples = (first_sample_upper+first_sample_lower)÷2
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2/delta) /num_samples) + log(2/delta)/num_samples/3
        if eps_guess > eps
            first_sample_lower = num_samples
        else
            first_sample_upper = num_samples
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("Maximum number of iterations "*string(omega))
    println("Initial sample size "*string(first_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("Initial sample size dropped to "*string(first_stopping_samples))
    end
    flush(stdout)
    next_stopping_samples = trunc(Int,first_stopping_samples)
    prev_stopping_samples::Int64 = 0
    local_temporal_betweenness = [zeros(tg.num_nodes) for _ in 1:ntasks]
    betweenness = zeros(tg.num_nodes)
    sampled_so_far::Int64 = 0
    stop::Array{Bool} = [false]
    sample_stop::Bool = false
    iteration::Int64 = 0
    while !sample_stop && !stop[1]
        approx_top_k = Array{Tuple{Int64,Float64}}([])
        sample_i = next_stopping_samples-prev_stopping_samples
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "trk"
                    _trk_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "ob"
                    _onbra_pfm_accumulate!(tg,tal,s,z,local_temporal_betweenness[t])
                elseif algo == "rtb"
                    _ssptp_accumulate!(tg,tal,s,local_temporal_betweenness[t])
                end
            end
        end
        sampled_so_far += sample_i
        if sampled_so_far >= omega
            sample_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(sampled_so_far)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end      
        betweenness = reduce(+, local_temporal_betweenness)

        if algo == "rtb"
            betweenness =  betweenness.*[1/(tg.num_nodes-1)]
        end
        for u in 1:tg.num_nodes
            push!(approx_top_k,(u,betweenness[u]))
        end
        sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
        _compute_finished_topk!(stop,omega,approx_top_k[begin:union_sample],sampled_so_far,eps,eps_lb,eps_ub,delta_lb_guess,delta_ub_guess,delta_lb_min_guess[1],delta_ub_min_guess[1],union_sample)   
        if !sample_stop && !stop[1]
            prev_stopping_samples = next_stopping_samples
            next_stopping_samples = trunc(Int,prev_stopping_samples*geo)
            iteration+=1
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-WUB-"*algo*"-PFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Increasing sample size to "*string(next_stopping_samples))
            flush(stdout)
            if (force_gc)
                clean_gc()
            end
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

