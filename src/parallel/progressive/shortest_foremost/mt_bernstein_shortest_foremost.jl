function threaded_progressive_bernstein_shortest_foremost(tg::temporal_graph,eps::Float64,delta::Float64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,geo::Float64 = 1.2,start_factor::Int64 = 100)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    local_temporal_betweenness_bootstrap::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for _ in 1:ntasks]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    xi::Float64 = 0
    balancing_factor::Float64 = 0.001
    omega::Int64 = 1000
    mc_trials::Int64 = 25
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    if (diam == -1)
        println("Approximating (sfm)-Temporal Diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_foremost_diameter(tg,64,0)
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
                _sfm_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _sfm_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "rtb"
                _sfm_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
        end
    end
    betweenness = reduce(+, local_temporal_betweenness_bootstrap)
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
        println("Average (sfm)-temporal path (upper bound) "*string(avg_diam_ub))
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
     local_temporal_betweenness_bootstrap = []
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
                if algo == "ob"
                    _p_onbra_sfm_bernstein_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "trk"
                    _p_trk_sfm_accumulate_bernstein!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "rtb"
                    _ssftp_accumulate_bernstein!(tg,tal,tn_index,bigint,s,local_temporal_betweenness[t],t_bc[t])
                end   
            end
        end
        sampled_so_far += sample_i
        #=
        if sampled_so_far >= omega
            sample_stop = true
            finish_partial = string(round(time() - start_time; digits=4))
            println("Completed, sampled "*string(sampled_so_far)*"/"*string(omega)* " couples in "*finish_partial)
            flush(stdout)
        end     
        =#
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        if algo == "rtb"
            xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc).*[1/(tg.num_nodes-1)],sampled_so_far,delta/2^iteration)
        else
            xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc),sampled_so_far,delta/2^iteration)
        end
        if xi <= eps
            stop[1] = true
        end
        if !sample_stop && !stop[1]
            prev_stopping_samples = next_stopping_samples
            next_stopping_samples = trunc(Int,prev_stopping_samples*geo)
            iteration+=1
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-Bernstein-"*algo*"-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = "*string(xi)*" | Increasing sample size to "*string(next_stopping_samples))
            flush(stdout)
        end
     end
     _reduce_list_of_arrays!(local_temporal_betweenness,betweenness,sampled_so_far,tg.num_nodes)
     betweenness = betweenness .* [1/sampled_so_far]
     return betweenness,[sampled_so_far],xi,time()-start_time 

end
function threaded_progressive_bernstein_shortest_foremost_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,bigint::Bool,algo::String = "trk",vc_upper_bund::Bool = true,diam::Int64 = -1,geo::Float64 = 1.2,start_factor::Int64 = 100)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    ntasks = nthreads()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    local_temporal_betweenness_bootstrap::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    reduced_betweenness::Vector{Float64} = Vector{Float64}([])
    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for _ in 1:ntasks]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    xi::Float64 = 0
    balancing_factor::Float64 = 0.001
    omega::Int64 = 1000
    mc_trials::Int64 = 25
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for _ in 1:ntasks]
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    if (diam == -1)
        println("Approximating (sfm)-Temporal Diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_foremost_diameter(tg,64,0)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        flush(stdout)
    end
    omega = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-1)))+1+log(2/delta)))
    union_sample::Int64 = min(tg.num_nodes,max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/ntasks),k+20))
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
                _sfm_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "ob"
                _sfm_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            elseif algo == "rtb"
                _sfm_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness_bootstrap[t],local_wv[t],mcrade[t],local_sp_lengths[t])
            end
        end
    end
    betweenness = reduce(+, local_temporal_betweenness_bootstrap)
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
        println("Average (sfm)-temporal path (upper bound) "*string(avg_diam_ub))
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
     local_temporal_betweenness_bootstrap = []
     betweenness = zeros(tg.num_nodes)
     sampled_so_far::Int64 = 0
     stop::Array{Bool} = [false]
     sample_stop::Bool = false
     iteration::Int64 = 0
     approx_top_k::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
     overall_converged::Int64 = 0
     while !sample_stop && !stop[1]
        sample_i = next_stopping_samples-prev_stopping_samples
        task_size = cld(sample_i, ntasks)
        vs_active = [i for i in 1:sample_i]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_i, task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "ob"
                    _p_onbra_sfm_bernstein_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "trk"
                    _p_trk_sfm_accumulate_bernstein!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "rtb"
                    _ssftp_accumulate_bernstein!(tg,tal,tn_index,bigint,s,local_temporal_betweenness[t],t_bc[t])
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
        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        betweenness = reduce(+,t_bc)
        if algo == "rtb"
            betweenness = betweenness.*[1/(tg.num_nodes-1)]
        end
        for u in 1:tg.num_nodes
            push!(approx_top_k,(u,betweenness[u]))
        end
        sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
        stop[1],overall_converged = theoretical_error_bound_topk(reduced_betweenness,approx_top_k[begin:union_sample],sampled_so_far,delta/2^iteration,eps,tg.num_nodes)
        if !sample_stop && !stop[1]
            prev_stopping_samples = next_stopping_samples
            next_stopping_samples = trunc(Int,prev_stopping_samples*geo)
            iteration+=1
            finish_partial = string(round(time() - start_time; digits=4))
            println("P-Bernstein-"*algo*"-SFM. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | top k converged "*string(overall_converged)*"/"*string(union_sample)*" | Increasing sample size to "*string(next_stopping_samples))
            flush(stdout)
        end
     end
     if stop[1]
        println("Progressive sampler converged at "*string(sampled_so_far)*"/"*string(omega)*" iterations")
        flush(stdout)
    end
    _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
    betweenness = reduce(+,t_bc)
    if algo == "rtb"
        betweenness = betweenness.*[1/(tg.num_nodes-1)]
    end
    approx_top_k = Array{Tuple{Int64,Float64}}([])
    for u in 1:tg.num_nodes
        push!(approx_top_k,(u,betweenness[u]/sampled_so_far))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])

    return approx_top_k,[sampled_so_far],time()-start_time 
end



function threaded_progressive_bernstein_shortest_foremost_dep(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool,algo::String = "trk")::Tuple{Array{Float64},Array{Int64},Float64,Float64}
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    start_time = time()
    print_algorithm_status("ONBRA","Bernstein",true)
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = Dict{Tuple{Int64,Int64},Int64}()
    if algo == "trk"
        tn_index = temporal_node_index_srtp(tg)
    else
        tn_index = temporal_node_index(tg)
    end
    ntasks = nthreads()

    local_temporal_betweenness::Vector{Vector{Vector{Float64}}} = [[] for _ in 1:ntasks]
    t_bc::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for _ in 1:ntasks]
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
    omega::Int64 = 1000
    t_diam::Float64 = 0.0
    
    println("Approximating diameter ")
    diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
    println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
    diam+=1
    
   
    omega = trunc(Int,(0.5/epsilon^2) * ((floor(log2(diam-2)))+log(1/delta)))
   
    println("Maximum sample size "*string(omega))
    println("Using ",nthreads()," Trheads")
    flush(stdout)

    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end
        task_size = cld(sample_size_schedule[j]-sample_size_schedule[j-1], ntasks)
        vs_active = [i for i in 1:sample_size_schedule[j]-sample_size_schedule[j-1]]
        @sync for (t, task_range) in enumerate(Iterators.partition(1:sample_size_schedule[j]-sample_size_schedule[j-1], task_size))
            Threads.@spawn for _ in @view(vs_active[task_range])
                sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
                s = sample[1][1]
                z = sample[1][2]
                if algo == "ob"
                    _p_onbra_sfm_bernstein_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "trk"
                    _p_trk_sfm_accumulate_bernstein!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[t],t_bc[t])
                elseif algo == "rtb"
                    _ssftp_accumulate_bernstein!(tg,tal,tn_index,bigint,s,local_temporal_betweenness[t],t_bc[t])
                end
            end
        end

        #=
        Base.Threads.@threads for i in 1:(sample_size_schedule[j]-sample_size_schedule[j-1])
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "ob"
                _p_onbra_sfm_bernstein_accumulate!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            elseif algo == "trk"
                _p_trk_sfm_accumulate_bernstein!(tg,tal,tn_index,bigint,s,z,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            elseif algo == "rtb"
                _ssftp_accumulate_bernstein!(tg,tal,tn_index,bigint,s,local_temporal_betweenness[Base.Threads.threadid()],t_bc[Base.Threads.threadid()])
            end
        end
        =#
        sampled_so_far+=sample_size_schedule[j]-sample_size_schedule[j-1]

        _reduce_arrays!(local_temporal_betweenness,reduced_betweenness)
        if algo == "rtb"
            xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc).*[1/(tg.num_nodes-1)],sample_size_schedule[j],delta/2^k)
        else
            xi = theoretical_error_bound(reduced_betweenness,reduce(+,t_bc),sample_size_schedule[j],delta/2^k)

        end


     
        finish_partial = string(round(time() - start_time; digits=4))
        println("P-"*algo*". Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        if xi <= epsilon || sampled_so_far >= omega
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
    if (bigint)
        bfs_ds = BI_BFS_SFM_ONBRA_DS(0,0)
    else
        bfs_ds = BFS_ONBRA_SFM_DS(0,0)
    end
    return nothing
end





function _p_trk_sfm_accumulate_bernstein!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    indexes::Int64 = length(keys(tn_index))
    if (bigint)
        bfs_ds = BI_BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
    else
        bfs_ds = BFS_SFM_SRTP_DS(tg.num_nodes, indexes)
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
        bfs_ds.t_hit[u] = -1
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
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.forward_queue, (s, 0,tni))
    d_z_min = Inf
    t_z_min = Inf
    while length(bfs_ds.forward_queue) != 0
        temporal_node = dequeue!(bfs_ds.forward_queue)
        u = temporal_node[1]
        t = temporal_node[2]
        tni = temporal_node[3]
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
                if bfs_ds.t_hit[w] == -1 || bfs_ds.t_hit[w] > t_w
                    bfs_ds.t_hit[w] = t_w
                end
            end
        end
    end
    if bfs_ds.dist[z] > 0
        
        totalWeight = 0
        randomEdge = 0
       
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
    if (bigint)
        bfs_ds = BI_BFS_SFM_SRTP_DS(0,0)
    else
        bfs_ds = BFS_SFM_SRTP_DS(0,0)
    end
    return nothing
end

function _ssftp_accumulate_bernstein!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},s::Int64,bigint::Bool,B_1::Vector{Vector{Float64}},B_2::Vector{Float64})
    if (bigint)
        bfs_ds = BI_BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_SFM_DS(tg.num_nodes, length(keys(tn_index)))
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
        bfs_ds.t_hit[u] = -1
        bfs_ds.dist[u] = -1
        bfs_ds.sigma[u] = 0
    end
    for tn in 1:(length(tn_index))
        bfs_ds.sigma_t[tn] = 0
        bfs_ds.delta_fm[tn] = 0
        bfs_ds.dist_t[tn] = -1
        bfs_ds.predecessors[tn] = Set{Tuple{Int64,Int64}}()
    end
    tni = tn_index[(s, 0)]
    bfs_ds.sigma[s] = 1
    bfs_ds.sigma_t[tni] = 1
    bfs_ds.dist[s] = 0
    bfs_ds.dist_t[tni] = 0
    bfs_ds.t_hit[s] = 0
    enqueue!(bfs_ds.queue, (s, 0))
    iter = 0
    while length(bfs_ds.queue) != 0
        iter += 1
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
            if bfs_ds.t_hit[w] == -1 || bfs_ds.t_hit[w] > t_w
                bfs_ds.t_hit[w] = t_w
            end
        end
    end
    sub = (count(x -> x >= 0, bfs_ds.dist) - 1) 
    B_1[p][s] -= sub *  (1/(tg.num_nodes-1))
    B_2[s] -= sub
    while length(bfs_ds.stack) != 0
        temporal_node = pop!(bfs_ds.stack)
        w = temporal_node[1]
        t_w = temporal_node[2]
        tni_w = tn_index[(w, t_w)]
        if bfs_ds.t_hit[w] == t_w
            bfs_ds.delta_fm[tni_w] += 1
        end
        for pred in bfs_ds.predecessors[tni_w]
            v = pred[1]
            t_v = pred[2]
            tni_v = tn_index[(v, t_v)]
            bfs_ds.delta_fm[tni_v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
            B_1[p][v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]*  (1/(tg.num_nodes-1))
            B_2[v] += (bfs_ds.sigma_t[tni_v] / bfs_ds.sigma_t[tni_w]) * bfs_ds.delta_fm[tni_w]
        end
    end
    if (bigint)
        bfs_ds = BI_BFS_SFM_DS(0,0)
    else
        bfs_ds = BFS_SFM_DS(0,0)
    end
    return nothing
end