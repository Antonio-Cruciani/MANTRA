function threaded_sstp_shortest_foremost(tg::temporal_graph,sample_size::Int64,verbose_step::Int64, bigint::Bool)::Tuple{Array{Float64},Float64}
    start_time::Float64 = time()
    tal::Array{Array{Tuple{Int64,Int64}}} = temporal_adjacency_list(tg)
    tn_index::Dict{Tuple{Int64,Int64},Int64}  = temporal_node_index(tg)
    local_temporal_betweenness::Vector{Vector{Float64}} = [zeros(tg.num_nodes) for i in 1:nthreads()]
    processed_so_far::Int64 = 0
    s::Int64 = 0
    sample::Array{Int64} = sstp_sample(tg, sample_size)
    println("Using ",nthreads()," Trheads")

 

    Base.Threads.@threads for i in 1:sample_size
        s = sample[i]
        _ssftp_accumulate!(tg,tal,tn_index,s,bigint,local_temporal_betweenness[Base.Threads.threadid()])
        processed_so_far = processed_so_far + 1
        if (verbose_step > 0 && processed_so_far % verbose_step == 0)
            finish_partial::String = string(round(time() - start_time; digits=4))
            time_to_finish::String = string(round((sample_size*(time() - start_time) / i )-(time() - start_time) ; digits=4))
            println("RTB-SFM. Processed " * string(processed_so_far) * "/" * string(sample_size) * " samples in " * finish_partial * " seconds | Est. remaining time : "*time_to_finish)
        end
    end
    betweenness = reduce(+, local_temporal_betweenness)
    betweenness = betweenness .* [1/(sample_size*(tg.num_nodes-1))]
    return betweenness,time()-start_time
end
