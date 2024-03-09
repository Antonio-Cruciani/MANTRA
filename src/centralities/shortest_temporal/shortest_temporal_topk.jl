


function shortest_temporal_betweenness_topk(tg::temporal_graph,eps::Float64,delta_1::Float64,delta_2::Float64,k::Int64,method::String, verbose_step::Int64, bigint::Bool,diam::Int64 = -1)
    start_time::Float64 = time()
    betweenness::Array{Float64}  = zeros(tg.num_nodes)
    ranking::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
    final_ranking::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
    eps_1::Float64 = 0.0
    eps_2::Float64 = 0.0
    b_k_1::Float64 = 0.0
    b_k_2::Float64 = 0.0
    if diam == -1
        println("Approximating diameter ")
        diam,_,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,verbose_step)
        println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
        diam+=1
    end
    diam = 18

    
    omega_1::Int64 = trunc(Int,(0.5/eps^2) * ((floor(log2(diam-2)))+1+log(1/delta_1)))

    println("Sample size for the first phase "*string(omega_1))
    if method == "rtb"
        betweenness,_ = rtb(tg,omega_1,verbose_step, bigint)
    elseif method == "trk"
        betweenness,_ = trk(tg,omega_1,verbose_step, bigint)

    elseif method == "ob"
        betweenness,_ = onbra_prefix_foremost(tg,omega_1,verbose_step)
    else
        println("Error choose one among: rtb , trk , ob")
        exit(1)
    end
    println("First phase completed in "*string(round(time() - start_time; digits=4)))
    for u in 1:tg.num_nodes
        push!(ranking,(u,betweenness[u]/omega_1))
    end
    sort!(ranking,by=ranking -> -ranking[2])
    println(ranking)
    b_k_1 = ranking[k][2]
    eps_1 = b_k_1 - eps
    println("EPS 1 ",eps_1," BK 1 ",b_k_1)
    omega_2::Int64 = trunc(Int,(0.5/(eps^2*eps_1)) * (((floor(log2(diam-2)))+1)*log(1/eps_1)+log(1/delta_2)))

    println("Sample size for the first phase "*string(omega_2))
    if method == "rtb"
        betweenness,_ = rtb(tg,omega_2,verbose_step, bigint)
    elseif method == "trk"
        betweenness,_ = trk(tg,omega_2,verbose_step, bigint)

    elseif method == "ob"
        betweenness,_ = onbra_prefix_foremost(tg,omega_2,verbose_step)
    end
    println("Second phase completed in "*string(round(time() - start_time; digits=4)))
    ranking = Array{Tuple{Int64,Float64}}([])
    for u in 1:tg.num_nodes
        push!(ranking,(u,betweenness[u]))
    end
    sort!(ranking,by=ranking -> -ranking[2])
    b_k_2 = ranking[k][2]/omega_2
    eps_2 = b_k_2/(1+eps)
    for u in 1:tg.num_nodes
        if ranking[u][2]/(1-eps) >= eps_2
            push!(final_ranking,ranking[u]/omega_2)
        end
    end
    term_time::Float64 = time()-start_time
    return final_ranking,time()-start_time
end