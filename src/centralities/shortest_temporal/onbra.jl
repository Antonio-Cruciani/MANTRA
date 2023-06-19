using DataStructures
using StatsBase
using NLopt
using JuMP

struct BFS_ONBRA_DS
    sigma::Array{UInt128}
    dist::Array{Int64}
    sigma_t::Array{UInt128}
    sigma_z::Array{UInt128}
    dist_t::Array{Int64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    boolean_matrix::Array{Bool}
    forward_queue::Queue{Tuple{Int64,Int64}}
    backward_queue::Queue{Tuple{Int64,Int64}}
    function BFS_ONBRA_DS(nn::Int64, ntn::Int64)
        return new(Array{UInt128}(undef, nn), Array{Int64}(undef, nn), Array{UInt128}(undef, ntn), zeros(Int64, ntn), Array{Int64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), falses(ntn), Queue{Tuple{Int64,Int64}}(), Queue{Tuple{Int64,Int64}}())
    end
end

struct BI_BFS_ONBRA_DS
    sigma::Array{BigInt}
    dist::Array{Int64}
    sigma_t::Array{BigInt}
    sigma_z::Array{BigInt}
    dist_t::Array{Int64}
    predecessors::Array{Set{Tuple{Int64,Int64}}}
    boolean_matrix::Array{Bool}
    forward_queue::Queue{Tuple{Int64,Int64}}
    backward_queue::Queue{Tuple{Int64,Int64}}
    function BI_BFS_ONBRA_DS(nn::Int64, ntn::Int64)
        return new(Array{BigInt}(undef, nn), Array{Int64}(undef, nn), Array{BigInt}(undef, ntn), zeros(Int64, ntn), Array{Int64}(undef, ntn), Array{Set{Tuple{Int64,Int64}}}(undef, ntn), falses(ntn), Queue{Tuple{Int64,Int64}}(), Queue{Tuple{Int64,Int64}}())
    end
end

function onbra_sample(tg::temporal_graph, sample_size::Int64)::Array{Tuple{Int64,Int64}}
    sample_pairs::Array{Tuple{Int64,Int64}} = []
    s::Int64 = 0
    z::Int64 = 0
    while length(sample_pairs) < sample_size
        s, z = sample(1:tg.num_nodes, 2, replace=false)
        push!(sample_pairs, (s, z))
    end
    return sample_pairs
end
#=
function empirical_variance(tilde_b::Array{Float64}, sample_size::Int64, v::Int64)::Float64
    n::Int64 = div(length(tilde_b), sample_size)
    variance::Float64 = 0
    for i in 1:sample_size
        for j in (i+1):sample_size
            variance += (((tilde_b[(i-1)*n+v] - tilde_b[(j-1)*n+v]) / sample_size)^2)
        end
    end
    return variance / (sample_size * (sample_size - 1))
end

function theoretical_error_bound(tilde_b::Array{Float64}, sample_size::Int64, eta::Float64)::Float64
    n::Int64 = div(length(tilde_b), sample_size)
    error::Float64 = 0.0
    max_error::Float64 = 0.0
    errors::Array{Float64} = zeros(n)
    for u in 1:n
        variance::Float64 = empirical_variance(tilde_b, sample_size, u)
        errors[u] = sqrt(2 * variance * log(4 * n / eta) / sample_size) + 7 * log(4 * n / eta) / (3 * (sample_size - 1))

    end
    return maximum(errors)
end
=#
function onbra(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Float64}
    start_time = time()
    sample = test_sample
    if (length(sample) == 0 || length(sample) != sample_size)
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    end
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    end
    tilde_b::Array{Float64} = zeros(tg.num_nodes)
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    processed_so_far::Int64 = 0
    exec_time::Array{Float64} = zeros(sample_size)
    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i][1]
        z = sample[i][2]
        for u in 1:tg.num_nodes
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
        enqueue!(bfs_ds.forward_queue, (s, 0))
        d_z_min = Inf
        while length(bfs_ds.forward_queue) != 0
            temporal_node = dequeue!(bfs_ds.forward_queue)
            u = temporal_node[1]
            t = temporal_node[2]
            tni = tn_index[(u, t)]
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
        end
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
                tilde_b[temporal_node[1]] += (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]))
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
        exec_time[i] = time() - exec_time[i]
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("ONBRA. Processed " * string(processed_so_far) * "/" * string(sample_size) * " pairs in " * finish_partial * " seconds")
        end
    end
    return tilde_b, time() - start_time
end


function def_summand(bigint::Bool)
    if bigint
        a::Float64 = 0
        b::BigInt = 0
        c::BigInt = 0
        return a,b,c
    end
    f::UInt128 = 0
    g::UInt128 = 0
    d::UInt128 = 0
    return f,g,d
end
function initialize_structures(bigInt::Bool,n::Int64)
    if bigInt
        a::Dict{Float64,Float64} = Dict()
        b::Array{Float64} =  zeros(n)
        return a,b
    end
    c::Dict{Float64,Float64} = Dict()
    d::Array{Float64} =  zeros(n)
    return c,d
end

function progressive_onbra(tg::temporal_graph, initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64 ,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Array{Int64},Float64,Float64}
    start_time = time()
    B,B_2 = initialize_structures(bigint,tg.num_nodes)
    B_1::Array{Float64} = zeros(tg.num_nodes)
    k::Int64 = 0
    j::Int64 = 2
    keep_sampling::Bool = true
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64} = temporal_node_index(tg)
    summand,b ,b_1 = def_summand(bigint)
    if (bigint)
        bfs_ds = BI_BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    else
        bfs_ds = BFS_ONBRA_DS(tg.num_nodes, length(keys(tn_index)))
    end
    u::Int64 = -1
    w::Int64 = -1
    t::Int64 = -1
    t_w::Int64 = -1
    tni::Int64 = -1
    tni_w::Int64 = -1
    temporal_node::Tuple{Int64,Int64} = (-1, -1)
    sampled_so_far::Int64 = 0
    reached_so_far::Int64 = 0
    new_sample::Int64 = 0
    sample_size_schedule::Array{Int64} = [0,initial_sample]
    xi::Float64 = 0
    exec_time::Array{Float64} =Array{Float64}([])
    exect_start::Float64 = -1.0
    while keep_sampling
        k+=1
        if (k >= 2)
            new_sample = trunc(Int,geo^k*sample_size_schedule[2])
            push!(sample_size_schedule,new_sample)
        end
        for i in 1:(sample_size_schedule[j]-sample_size_schedule[j-1])
            exect_start = time()
            sampled_so_far+=1
            sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, 1)
            s = sample[1][1]
            z = sample[1][2]
            for u in 1:tg.num_nodes
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
            enqueue!(bfs_ds.forward_queue, (s, 0))
            d_z_min = Inf
            while length(bfs_ds.forward_queue) != 0
                temporal_node = dequeue!(bfs_ds.forward_queue)
                u = temporal_node[1]
                t = temporal_node[2]
                tni = tn_index[(u, t)]
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
                    end
                end
            end
            if bfs_ds.dist[z] > 0
                reached_so_far+=1
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
                            summand = Float64(bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z])    , RoundUp)
                        else
                            summand = (bfs_ds.sigma_z[tni] * (bfs_ds.sigma_t[tni] / bfs_ds.sigma[z]) )
                        end

                        # Updating phase
                        b = B_2[temporal_node[1]]
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
                        B_1[temporal_node[1]] += summand
                        B_2[temporal_node[1]] += summand^2

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
            push!(exec_time, time() - exect_start)
        end
        B_vectorized = collect(keys(B))
        if length(B_vectorized)>0
            omega = compute_xi(B_vectorized,sample_size_schedule[j])
            
            delta_i = delta/2^k
            xi = 2*omega + (log(3/delta_i)+sqrt((log(3/delta_i)+4*sample_size_schedule[j]*omega)*log(3/delta_i)))/sample_size_schedule[j]  + sqrt(log(3/delta_i)/(2*sample_size_schedule[j]))
        else
            xi = Inf
        end
        if (verbose_step > 0 && j % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            println("P-ONBRA. Processed " * string(sampled_so_far) * " pairs in " * finish_partial * " seconds | Est. ξ = ",xi)
        end
        if xi <= epsilon
            keep_sampling = false
        else
            j+=1
        end

    end
    
    finish::Float64 = round(time() - start_time; digits=4)
    println("ONBRA-SH. Total number of samples ",sampled_so_far," ξ = ",xi, " Time ",string(finish))

    #betweenness::Array{Float64} = [B_1[i] / sample_size_schedule[j] for i in 1:lastindex(B_1)]
    return B_1,sample_size_schedule,xi,time()-start_time
end

function compute_xi(B,r)
    myf(x)= 1/x * log(sum([exp(x^2 * B[i]/(2*r^2)) for i in 1:lastindex(B)])) 
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    local_optimizer = NLopt.Opt(:LD_LBFGS, 1)
    set_optimizer_attribute(model, "local_optimizer", :LD_LBFGS)
    local_optimizer.lower_bounds = [floatmin(Float64)]
    local_optimizer.xtol_abs = 1e-4
    local_optimizer.ftol_abs =1e-6
    set_optimizer_attribute(model, "local_optimizer", local_optimizer)
    #set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    register(model,:myf, 1, myf, autodiff=true)
    @variable(model, x >= 0)
    @constraint(model,x >= 0)
    @NLobjective(model, Min,myf(x))
    JuMP.optimize!(model)
    return(objective_value(model))
end


