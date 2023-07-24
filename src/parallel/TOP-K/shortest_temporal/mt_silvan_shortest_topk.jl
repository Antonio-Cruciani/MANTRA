function threaded_progressive_silvan_topk(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64,verbose_step::Int64,bigint::Bool,diam::Int64 = -1,algo::String = "trk",empirical_peeling_a::Float64 = 2.0,sample_step::Int64 = 10,hb::Bool = false)
    @assert (algo == "trk") || (algo == "ob") || (algo == "rtb") "Illegal algorithm, use: trk , ob , or rtb"
    norm::Float64 = 1.0
    if algo == "rtb"
        norm = 1/(tg.num_nodes-1)
    end
    start_time = time()
    union_sample::Int64 = min(tg.num_nodes,max(sqrt(lastindex(tg.temporal_edges))/nthreads(),k+20))
    mc_trials::Int64 = 25
    # Temoral Graph structures
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index_srtp(tg)
    # Wimpy variance
    local_wv::Array{Array{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    wv::Array{Float64} = Array{Float64}([])
    emp_w_node::Float64 = 0.0
    min_inv_w_node::Float64 = 0.0
    # Partitions
    number_of_non_empty_partitions::Int64 = 0.0
    non_empty_partitions::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partitions_ids_map::Dict{Int64,Int64} = Dict{Int64,Int64}()
    partition_index::Array{Int64} = zeros(tg.num_nodes)
    part_idx::Int64 = 1
    # TBC
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade::Array{Array{Float64}} = [zeros(tg.num_nodes*mc_trials) for i in 1:nthreads()]
    local_sp_lengths::Array{Array{Int64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega::Float64 = 1000
    t_diam::Float64 = 0.0
    max_num_samples::Float64 = 0.0
    top_k_result::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
    #diam::Float64 = 0.0
    if (diam == -1) 
        println("Approximating diameter ")
        diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step,0.9,false)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
        diam+=1
    end
    #max_num_samples = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    tau::Int64 = trunc(Int64,max(1. / eps * (log(1. / delta)) , 100.))
    tau =  trunc(Int64,max(tau,2*(diam -1) * (log(1. / delta))) )
    s::Int64 = 0
    z::Int64 = 0
    bootstrap_finished::Bool = false
    approx_to_add::Array{Int64} = zeros(tg.num_nodes)
    distinct_nodes_top_k::Array{Int64} = [0]
    tbc::Array{Float64} = zeros(tg.num_nodes)
    println("Bootstrap phase "*string(tau)*" iterations")
    n_pairs = 0
    while !bootstrap_finished
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk_topk!(tg,tal,tn_index,bigint,s,z,mc_trials,true,tbc,local_wv[Base.Threads.threadid()],distinct_nodes_top_k,approx_to_add)
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            else
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,true,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            end
        end
        #betweenness = reduce(+, local_temporal_betweenness)
        #betweenness = betweenness .* [1/tau]
        #
        n_pairs += sample_step
        bootstrap_finished = (distinct_nodes_top_k[1] >= 1.5*k )&&(n_pairs >= tau)
        println(" distinct_nodes_top_k ",distinct_nodes_top_k[1])
    end
    tau = n_pairs
    betweenness = tbc .*[1/tau]
    println("Empirical peeling phase:")
    max_tbc::Float64 = 0.0
    max_wv::Float64 = 0.0
    wv = reduce(+,local_wv)

    for i in 1:tg.num_nodes
        max_tbc = max(max_tbc,betweenness[i])
        max_wv = max(max_wv,wv[i])
        emp_w_node = wv[i] * 1. /tau
        min_inv_w_node = min(1. /emp_w_node,tau)
        node_partition_idx = trunc(Int,log(min_inv_w_node)/log(empirical_peeling_a)+1)
        partition_index[i] = node_partition_idx
        if haskey(non_empty_partitions,node_partition_idx)
            non_empty_partitions[node_partition_idx] += 1
        else
            non_empty_partitions[node_partition_idx] = 1
        end
    end
    number_of_non_empty_partitions = 0
    for key in keys(non_empty_partitions)
        partitions_ids_map[key] = part_idx
        part_idx+=1
        number_of_non_empty_partitions+=1
    end
    println("Number of non empty partitions "*string(number_of_non_empty_partitions))
    println("Bootstrap completed in "*string(round(time() - start_time; digits=4)))
    for key in keys(non_empty_partitions)
        println(" Part w. index "*string(key)*" has "*string(non_empty_partitions[key])*" elements, map to "*string(partitions_ids_map[key]))
    end
    #=
    # Upper bound on the average distance
    avg_diam_ub::Float64 = upper_bound_average_diameter(delta/8,diam,sp_lengths,tau,true,norm)
    # Upper bound on the top-1 temporal betweenness
    top1bc_upper_bound::Float64 = upper_bound_top_1_tbc(max_tbc,delta/8,tau)
    wimpy_var_upper_bound::Float64 = upper_bound_top_1_tbc(max_wv/tau,delta/8,tau)
    # define delta_for_progressive_bound
    println("AVERAGE DIAM UB "*string(avg_diam_ub))
    # Upper limit on number of samples
    max_num_samples = upper_bound_samples(top1bc_upper_bound,wimpy_var_upper_bound,avg_diam_ub,eps,delta/2 ,false)
    omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
    println("Maximum number of samples "*string(max_num_samples)*" VC Bound "*string(omega))
    println("Sup tbc est "*string(max_tbc))
    println("Sup emp wimpy variance "*string(max_wv/tau))
    =#
    iteration_index::Int64 =1 
    est_kth_bc::Float64 = 10.0
    for i in 1:tg.num_nodes
        push!(approx_top_k,(i,betweenness[i]))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    est_kth_bc = approx_top_k[k][2] +1/tau
    approx_top_k =  Array{Tuple{Int64,Float64}}([])
    println("Est k-th tbc "*string(est_kth_bc))
    local_wv = [zeros(tg.num_nodes) for i in 1:nthreads()]
    local_temporal_betweenness = [zeros(tg.num_nodes) for i in 1:nthreads()]
    mcrade = [zeros((tg.num_nodes+1)*mc_trials) for i in 1:nthreads()]
    local_sp_lengths = [zeros(tg.num_nodes) for i in 1:nthreads()]
    omega = 10^15
    if max_num_samples > 0
        omega = max_num_samples
    end
    first_stopping_samples::Float64 = 0.0
    eps_guess::Float64 = 1.0
    first_sample_lower::Float64 = 1/eps *log(2/delta)
    first_sample_upper::Float64 = omega
    println("First sample lower "*string(first_sample_lower)*"First sample upper "*string(first_sample_upper))
    sup_emp_wimpy_var_norm::Float64  = max_wv/tau +1/tau
    violated_condition_lower_1::Bool = false
    violated_condition_lower_2::Bool = false
    violated_condition_upper_1::Bool = false
    while first_sample_upper - first_sample_lower > 10
        num_samples = (first_sample_upper+first_sample_lower)/2
        eps_guess = sqrt(2*est_kth_bc*log(2/delta/est_kth_bc) /num_samples) + log(2/delta/est_kth_bc)/num_samples/3
        println("Eps guess "*string(eps_guess))
        violated_condition_lower_1 = (est_kth_bc + eps_guess > est_kth_bc/(1-eps))
        violated_condition_upper_1 = (est_kth_bc - eps_guess < est_kth_bc/(1+eps))
        violated_condition_lower_2 = (eps_guess >= (est_kth_bc *eps*((1-eps)/(1+eps))^2))
        println(violated_condition_lower_1, " ",violated_condition_upper_1," ",violated_condition_lower_2)
        if violated_condition_lower_1 || violated_condition_upper_1 || violated_condition_lower_2
            first_sample_lower = num_samples
            println(" (top-k) lower increased at: "*string(first_sample_lower))
        else
            first_sample_upper = num_samples
            println(" (top-k) upped increased at: "*string(first_sample_upper))
        end
    end
    first_stopping_samples = num_samples
    last_stopping_samples = omega
    println("First stopping samples "*string(first_stopping_samples)* " Last stopping samples "*string(last_stopping_samples))
    if first_stopping_samples >= last_stopping_samples/4
        first_stopping_samples = last_stopping_samples/4
        println("First stopping samples dropped to "*string(first_stopping_samples))
    end
    next_stopping_samples::Float64 = first_stopping_samples
    
    has_to_stop::Bool = false
    num_samples::Float64 = 0.0
    while !has_to_stop
        Base.Threads.@threads for i in 1:sample_step
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            if algo == "trk"
                _sh_accumulate_trk!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            elseif algo == "ob"
                _sh_accumulate_onbra!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            else
                _sh_accumulate_rtb!(tg,tal,tn_index,bigint,s,z,mc_trials,false,local_temporal_betweenness[Base.Threads.threadid()],local_wv[Base.Threads.threadid()],mcrade[Base.Threads.threadid()],local_sp_lengths[Base.Threads.threadid()])
            end
        end
        num_samples += sample_step
        #println("Num s "*string(num_samples)*" "*string(has_to_stop)* " "*string((num_samples < last_stopping_samples))*"  "*string((num_samples >= next_stopping_samples)))
        #println(" "*string(num_samples)*" "*string(last_stopping_samples)*" FRST = "*string(first_stopping_samples))
       if num_samples >= omega
            println("Num samples/ω : "*string(num_samples)*"/"*string(omega))
            has_to_stop = true
        end
        #println("Checking stopping condition")
        #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
        #println(" num_samples ",num_samples,"  ",next_stopping_samples)
        # & (num_samples >= next_stopping_samples)
        if !has_to_stop & (num_samples < last_stopping_samples)&(num_samples >= next_stopping_samples)
            betweenness = reduce(+, local_temporal_betweenness)
            approx_top_k =  Array{Tuple{Int64,Float64}}([])

            wv = reduce(+,local_wv)
            sp_lengths = reduce(+,local_sp_lengths) 
            r_mcrade = reduce(+,mcrade)
            println("Checking stopping condition")
            for i in 1:tg.num_nodes
                push!(approx_top_k,(i,betweenness[i]/num_samples))
            end
            sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
            #println(" num_samples ",num_samples," last_stopping_samples ",last_stopping_samples)
            #println(" num_samples ",num_samples,"  ",next_stopping_samples)
            #(approx_top_k,union_sample,num_samples,eps,delta,iteration_index,partition_index,partitions_ids_map,number_of_non_empty_partitions,mc_trials,sup_bcest_partition,sup_empwvar_partition,norm)
            has_to_stop,top_k_result  = check_stopping_condition_topk(betweenness,k,union_sample,num_samples,eps,delta,iteration_index,mc_trials,r_mcrade,norm,partition_index, number_of_non_empty_partitions, non_empty_partitions,partitions_ids_map,wv,diam,sp_lengths)

            if has_to_stop
                println("Progressive sampler converged!")
            else
                next_stopping_samples,iteration_index = get_next_stopping_sample(next_stopping_samples,iteration_index )
                println("Increasing sample size to "*string(next_stopping_samples))
            end
                    
        end


    end

    println("(SH)-Temporal Betweenness estimated in "*string(round(time() - start_time; digits=4)))

    return top_k_result,num_samples,time()-start_time


end



function check_stopping_condition_topk(betweenness::Array{Float64},k::Int64,union_sample::Int64,num_samples::Float64,eps::Float64,delta::Float64,iteration::Int64,mc_trials::Int64,mcrade::Array{Float64},norm::Float64,partition_index::Array{Int64},number_of_non_empty_partitions::Int64, non_empty_partitions::Dict{Int64,Int64},partitions_ids_map::Dict{Int64,Int64},emp_wimpy_vars::Array{Float64},diam,sp_lengths::Array{Float64})
    n::Int64 = lastindex(betweenness)
    approx_top_k::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    for i in 1:tg.num_nodes
        push!(approx_top_k,(i,betweenness[i]))
    end
    sort!(approx_top_k, by=approx_top_k->-approx_top_k[2])
    num_samples_d::Float64 = num_samples
    delta_for_progressive_bound::Float64 = delta/2^iteration
    println("Checking stopping condition at iteration "*string(iteration)*" sample size "*string(num_samples)*" δ = "*string(delta_for_progressive_bound))

    sup_bcest_partition::Array{Float64} = zeros(n)
    sup_empwvar_partition::Array{Float64} = zeros(n)
    epsilon_partition::Array{Float64} = ones(n)
    max_mcera_partition::Array{Float64} = [-num_samples for i in 1:(mc_trials*n) ]
    # Update MCERA
    for j in 1:k
        i = approx_top_k[j][1]
        v_rade_idx = i*mc_trials
        node_partition_idx = partition_index[i]
        mapped_partition_index = partitions_ids_map[node_partition_idx]
        sup_bcest_partition[mapped_partition_index] = max(sup_bcest_partition[mapped_partition_index],betweenness[i] )
        sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index],emp_wimpy_vars[i])
        mcera_partition_index = mc_trials*mapped_partition_index
        for j in 1:mc_trials
            max_mcera_partition[j+mcera_partition_index] = max(max_mcera_partition[j+mcera_partition_index] ,mcrade[v_rade_idx+j])
        end
    end
    mcera_partition_avg::Array{Float64} = zeros(number_of_non_empty_partitions)
    mcera_avg::Float64 = 0.0
    mcera_partition_index::Int64 = 0.0 
    delta_each_partition::Float64 = delta_for_progressive_bound/number_of_non_empty_partitions
    for i in 1:number_of_non_empty_partitions
        mcera_avg = 0.0
        mcera_partition_index = mc_trials *i
        for j in 1:mc_trials
            mcera_avg+=max_mcera_partition[mcera_partition_index+j]/mc_trials
        end
        mcera_avg = mcera_avg/num_samples_d
        mcera_partition_avg[i] = mcera_avg
        sup_emp_wimpy_var = sup_empwvar_partition[i]/num_samples_d
        current_eps = epsilon_mcrade(sup_emp_wimpy_var,mcera_avg,delta_each_partition,num_samples_d,mc_trials)
        epsilon_partition[i] = current_eps
    end
    top_k_converged,top_k_result = check_top_k_convergence(approx_top_k,union_sample,num_samples,eps,delta,iteration,partition_index,partitions_ids_map,number_of_non_empty_partitions,mc_trials,sup_bcest_partition,sup_empwvar_partition,norm,emp_wimpy_vars,max_mcera_partition,mcrade,diam,sp_lengths,k)
    if top_k_converged
        println("MCRADE STOPS at iteration : "*string(iteration))
        println("MCRADE STOPS at sample size : "*string(num_samples))
    else
     println("MCRADE DOEST NOT STOP" )
    end
    return top_k_converged,top_k_result
end

function compute_relative_bound(n::Int64,bc_est::Float64,delta::Float64,avg_dist::Float64,m_samples::Float64,upper_b::Bool)::Float64
    quant_condition::Float64 = 0.0
    upper::Float64 = 1.0
    lower::Float64 = 0.0
    v_test::Float64 = 0.0
    log_term::Float64 = 0.0
    sqrt_term::Float64 = 0.0
    if upper_b
        lower = bc_est
    else
        upper = bc_est
    end
    test::Float64 = 1.0
    i::Int64 = 0
    while i < 30
        test = (upper+lower)÷2
        v_test = test*(1-test)
        log_term = log(2/delta*min(n),avg_dist/test)/m_samples
        sqrt_term = sqrt(2*log_term*v_test)
        if upper_b
            quant_condition = bc_est + log_term /3 + sqrt_term
            if quant_condition >= test
                lower = test
            else
                upper = test
            end
        else
            quant_condition = bc_est - log_term /3 - sqrt_term
            if quant_condition <= test
                upper = test
            else
                lower = test
            end
        end
        i+=1
    end
    if upper_b
        return upper
    end
    return lower
end

function check_top_k_convergence(approx_top_k::Array{Tuple{Int64,Float64}},union_sample::Int64,num_samples::Float64,err::Float64,delta::Float64,iteration::Int64,partition_index::Array{Int64},partitions_ids_map::Dict{Int64,Int64},number_of_non_empty_partitions::Int64,mc_trials::Int64,sup_bcest_partition::Array{Float64},sup_empwvar_partition::Array{Float64},norm::Float64,emp_wimpy_vars::Array{Float64},max_mcera_partition,mcrade::Array{Float64},diam,sp_lengths::Array{Float64},k::Int64)
    num_samples_d::Float64 = num_samples
    n::Int64 = lastindex(approx_top_k)
    println("Evaluating stopping condition for TOP-K relative approximation")
    delta_for_progressive_bound::Float64 = delta/2^iteration
    eps_final_topk::Float64 = 1.0
    top_k_tbc_lb::Float64 = 0.0
    top_k_tbc_ub::Float64 = 1.0
    top_k_result::Array{Tuple{Int64,Float64,Float64,Float64}} =  Array{Tuple{Int64,Float64,Float64,Float64}}([])
    epsilon_partition::Array{Float64} = zeros(union_sample)
    approx_top_k_lb::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    approx_top_k_ub::Array{Tuple{Int64,Float64}} =  Array{Tuple{Int64,Float64}}([])
    num_inserted::Int64 = 0
    v::Int64 = 0
    is_relative_approx::Bool = true
    avg_diam_upperbound::Float64 = 0.0
    approx_v::Float64 = 0.0

    for i in 1:number_of_non_empty_partitions
        v = approx_top_k[i][1]
        v_rade_idx = v*mc_trials
        node_partition_idx = partition_index[v]
        mapped_partition_index = partitions_ids_map[node_partition_idx]
        sup_bcest_partition[mapped_partition_index] = max(sup_bcest_partition[mapped_partition_index],approx_top_k[i][2])
        sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index],emp_wimpy_vars[i])
        mcera_partition_index = mc_trials*mapped_partition_index
        for j in 1:mc_trials
            max_mcera_partition[j+mcera_partition_index] = max(max_mcera_partition[j+mcera_partition_index] ,mcrade[v_rade_idx+j])
        end
    end
    mcera_partition_avg::Array{Float64} = zeros(number_of_non_empty_partitions)
    mcera_avg::Float64 = 0.0
    mcera_partition_index::Int64 = 0.0 
    delta_each_partition::Float64 = delta_for_progressive_bound/number_of_non_empty_partitions
    for i in 1:number_of_non_empty_partitions
        mcera_avg = 0.0
        mcera_partition_index = mc_trials *i
        for j in 1:mc_trials
            mcera_avg+=max_mcera_partition[mcera_partition_index+j]/mc_trials
        end
        mcera_avg = mcera_avg/num_samples_d
        mcera_partition_avg[i] = mcera_avg
        sup_emp_wimpy_var = sup_empwvar_partition[i]/num_samples_d
        current_eps = epsilon_mcrade(sup_emp_wimpy_var,mcera_avg,delta_each_partition,num_samples_d,mc_trials)
        epsilon_partition[i] = current_eps
    end
    for i in 1:union_sample
        v = approx_top_k[i][1]
        approx_v = approx_top_k[i][2]/num_samples_d
        node_partition_idx = partition_index[v]
        map_node_partition_index = partitions_ids_map[node_partition_idx]
        eps_current_node = epsilon_partition[map_node_partition_index]
        lowerbound_bc = approx_v-eps_current_node
        upperbound_bc = approx_v+eps_current_node
        avg_diam_upperbound = upper_bound_average_diameter(delta_for_progressive_bound,trunc(Int,diam),sp_lengths,trunc(Int,num_samples_d),true,norm)
        lb_tbc_rel = compute_relative_bound(n,approx_v,delta_for_progressive_bound,avg_diam_upperbound,num_samples_d,false)
        ub_tbc_rel = compute_relative_bound(n,approx_v,delta_for_progressive_bound,avg_diam_upperbound,num_samples_d,true)
        lowerbound_bc = max(lowerbound_bc,lb_tbc_rel)
        upperbound_bc = min(upperbound_bc,ub_tbc_rel)
        push!(approx_top_k_lb,(v,lowerbound_bc))
        push!(approx_top_k_ub,(v,upperbound_bc))
        num_inserted += 1
        if num_inserted >= k
            top_k_tbc_ub = min(top_k_tbc_ub,approx_top_k_ub[k][2])
            top_k_tbc_lb = max(top_k_tbc_lb,approx_top_k_lb[k][2])
        end
        if upperbound_bc > top_k_tbc_lb
            eps_final_topk = min(eps_final_topk,eps_current_node)
            l_condition = (lowerbound_bc >= approx_v/(1+err))
            u_condition = (upperbound_bc <= approx_v/(1-err))
            l_condition_rel = (lb_tbc_rel >= approx_v/(1+err))
            u_condition_rel = (ub_tbc_rel <= approx_v/(1-err))
            if !l_condition
                println(" Not stopping as current node has to improve its lower bound")
                println(" Lower bound Temporal BC "*string(lowerbound_bc))
                println(" Lower bound Temporal BC (relative) "*string(lb_tbc_rel))
                println(" Relative condition "*string(l_condition_rel))
                println(" approximation_node / (1+ε) "*string(approx_v/(1+err)))
            end
            if !u_condition
                println(" Not stopping as current node has to improve its upper bound")
                println(" Upper bound Temporal BC "*string(upperbound_bc))
                println(" Upper bound Temporal BC (relative) "*string(ub_tbc_rel))
                println(" Relative condition "*string(u_condition_rel))
                println(" approximation_node / (1-ε) "*string(approx_v/(1-err)))
            end
            if !l_condition || !u_condition
                println("Approx. node "*string(approx_v))
                println("Node partition index "*string(node_partition_idx))
            end
            is_relative_approx = is_relative_approx & l_condition & u_condition
        else
            eps_final_topk = min(eps_final_topk,top_k_tbc_lb-approx_v)
        end
    end
    println(" Upper bound top-k TBC "*string(top_k_tbc_ub))
    println(" Lower bound top-k TBC "*string(top_k_tbc_lb))
    if is_relative_approx
        println("Relative Approximation for the TOP-k Nodes achieved")
        for i in 1:num_inserted
            v = approx_top_k[i][1]
            approx_v = approx_top_k[i][2]/num_samples_d
            node_partition_idx = partition_index[v]
            map_node_partition_index = partitions_ids_map[node_partition_idx]
            eps_current_node = epsilon_partition[map_node_partition_index]
            lowerbound_bc = approx_v-eps_current_node
            upperbound_bc = approx_v+eps_current_node
            if upperbound_bc > lowerbound_bc
                push!(top_k_result,(v,approx_v,lowerbound_bc,upperbound_bc))
            end
        end
    end

    return is_relative_approx,top_k_result


end



function _sh_accumulate_trk_topk!(tg::temporal_graph,tal::Array{Array{Tuple{Int64,Int64}}},tn_index::Dict{Tuple{Int64,Int64},Int64},bigint::Bool,s::Int64,z::Int64,mc_trials::Int64,boostrap_phase::Bool,temporal_betweenness_centrality::Vector{Float64},wimpy_variance::Vector{Float64},distinct_nodes_top_k::Vector{Int64},approx_to_add::Array{Int64})
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
                        wimpy_variance[pred[1]] += 1
                        lk = ReentrantLock()
                        lock(lk)
                        try 
                            temporal_betweenness_centrality[pred[1]] += 1
                       
                            if (approx_to_add[pred[1]] == 0) && (temporal_betweenness_centrality[pred[1]] >= 3)
                                approx_to_add[pred[1]] = 1
                                distinct_nodes_top_k[1] +=1
                            end
                        finally 
                            unlock(lk)
                        end
                        
                    end
                    break
                end
            end
            
        end
        
    end
    return nothing

end