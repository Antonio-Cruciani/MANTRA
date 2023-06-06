using DataStructures
using StatsBase

struct BFS_ONBRA_PFM_DS
    t_hit::Array{Int64}
    sigma::Array{UInt128}
    sigma_z::Array{UInt128}
    predecessors::Array{Set{Int64}}
    boolean_matrix::Array{Bool}
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    backward_queue::Queue{Int64}
    function BFS_ONBRA_PFM_DS(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{UInt128}(undef,nn),Array{UInt128}(undef,nn),Array{Set{Int64}}(undef,nn),falses(nn), PriorityQueue{Tuple{Int64,Int64,Int64},Int64}(),Queue{Int64}() )
    end
end


struct BFS_ONBRA_PFM_DS_BI
    t_hit::Array{Int64}
    sigma::Array{BigInt}
    sigma_z::Array{BigInt}
    predecessors::Array{Set{Int64}}
    boolean_matrix::Array{Bool}
    priority_queue::PriorityQueue{Tuple{Int64,Int64,Int64},Int64}
    backward_queue::Queue{Int64}
    function BFS_ONBRA_PFM_DS_BI(nn::Int64)
        return new(Array{Int64}(undef,nn),Array{BigInt}(undef,nn),Array{BigInt}(undef,nn),Array{Set{Int64}}(undef,nn),falses(nn),  PriorityQueue{Tuple{Int64,Int64,Int64},Int64}(),Queue{Int64}() )
    end
end


function onbra_prefix_foremost(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Float64}
    start_time = time()
    sample = test_sample
    if (length(sample) == 0 || length(sample) != sample_size)
        sample::Array{Tuple{Int64,Int64}} = onbra_sample(tg, sample_size)
    end
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    if (bigint)
        bfs_ds = BFS_ONBRA_PFM_DS_BI(tg.num_nodes)
    else
        bfs_ds = BFS_ONBRA_PFM_DS(tg.num_nodes)
    end
    
    tilde_b::Array{Float64} = zeros(tg.num_nodes)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    processed_so_far::Int64 = 0
    exec_time::Array{Float64} = zeros(sample_size)
    for i in 1:sample_size
        exec_time[i] = time()
        s = sample[i][1]
        z = sample[i][2]
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
        end
        for pred in bfs_ds.predecessors[z]
            bfs_ds.sigma_z[pred] =1
            if !bfs_ds.boolean_matrix[pred]
                enqueue!(bfs_ds.backward_queue,pred)
                bfs_ds.boolean_matrix[pred] = true
            end
        end
        while length(bfs_ds.backward_queue) > 0
            v = dequeue!(bfs_ds.backward_queue)
            if v != s
                tilde_b[v] += (bfs_ds.sigma[v] * (bfs_ds.sigma_z[v] / bfs_ds.sigma[z]))
                for pred in bfs_ds.predecessors[v]
                    if (!bigint && bfs_ds.sigma_z[pred] > typemax(UInt128) - bfs_ds.sigma_z[pred])
                        println("Overflow occurred with sample (", s, ",", z, ")")
                        return [], 0.0
                    end
                    bfs_ds.sigma_z[pred] += bfs_ds.sigma_z[v]
                    if !bfs_ds.boolean_matrix[pred]
                        enqueue!(bfs_ds.backward_queue, pred) 
                        bfs_ds.boolean_matrix[pred] = true
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
    return tilde_b,  time() - start_time
end



function progressive_onbra_prefix_foremost(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64 ,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])::Tuple{Array{Float64},Array{Int64},Float64, Float64}
    start_time = time()
    B,B_2 = initialize_structures(bigint,tg.num_nodes)
    B_1::Array{Float64} = zeros(tg.num_nodes)
    k::Int64 = 0
    j::Int64 = 2
    keep_sampling::Bool = true
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    if (bigint)
        bfs_ds = BFS_ONBRA_PFM_DS_BI(tg.num_nodes)
    else
        bfs_ds = BFS_ONBRA_PFM_DS(tg.num_nodes)
    end

    summand,b ,b_1 = def_summand(bigint)
    w::Int64 = -1
    v::Int64 = -1
    t_w::Int64 = -1
    temporal_edge::Tuple{Int64,Int64,Int64} = (-1,-1,-1)
    exec_time::Array{Float64} = Array{Float64}([])
    sampled_so_far::Int64 = 0
    new_sample::Int64 = 0
    sample_size_schedule::Array{Int64} = [0,initial_sample]
    xi::Float64 = 0
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
    println("P-ONBRA. Total number of samples ",sampled_so_far," ξ = ",xi, " Time ",string(finish))

    #betweenness::Array{Float64} = [B_1[i] / sample_size_schedule[j] for i in 1:lastindex(B_1)]

    return B_1,sample_size_schedule,xi, time() - start_time


end