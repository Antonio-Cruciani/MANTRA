function read_centrality_values(file_name::String)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist "*file_name
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = []
    value::Float64 = 0.0
    for l in eachline(f)
        value = parse(Float64, l)
        if (value < -0.1)
            println("ERROR. There are negative values with absolute big values")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, value)
    end
    close(f)
    return centrality
end

function read_centrality_values_topk(file_name::String)::Array{Tuple{Int64,Float64}}
    @assert isfile(file_name) "The centrality value file does not exist"
    f::IOStream = open(file_name, "r")
    centrality::Array{Tuple{Int64,Float64}} = []
    value::Float64 = 0.0
    for l in eachline(f)
        split_line::Vector{String} = split(l, " ")
        value = parse(Float64, split_line[2])
        node = parse(Float64, split_line[1])

        if (value < -0.1)
            println("ERROR. There are negative values with absolute big values")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, (node,value))
    end
    close(f)
    return centrality
end

function _reduce_data_a!(u::Int64,tn::Int64,src_x::Vector{Vector{Float64}},src_y::Vector{Array{Float64}},src_z::Vector{Array{Int64}},dst_x::Array{Float64},dst_y::Array{Float64},dst_z::Array{Int64})
    for t in 1:tn
        dst_x[u]+=src_x[t][u]
        dst_y[u]+=src_y[t][u]
        dst_z[u]+=src_z[t][u]
    end
    return nothing
end

function _reduce_data_b!(u::Int64,tn::Int64,src_x::Vector{Vector{Float64}},src_y::Vector{Array{Float64}},src_z::Array{Array{Float64}},dst_x::Array{Float64},dst_y::Array{Float64},dst_z::Array{Float64})
    for t in 1:tn
        dst_x[u]+=src_x[t][u]
        dst_y[u]+=src_y[t][u]
        dst_z[u]+=src_z[t][u]
    end
    return nothing
end




function save_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "w")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end

function append_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "a")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end

function append_centrality_values_topk(file_name::String, centrality::Array{Tuple{Int64,Float64}})::Nothing
    f::IOStream = open(file_name, "a")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u][1])*" "*string(centrality[u][2]) * "\n")
    end
    close(f)
end

function read_time(file_name::String)
    @assert isfile(file_name) "The time value file does not exist: "*string(file_name)
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[1]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end

function read_sample_size(file_name::String)
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[2]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end
function read_xi(file_name::String)
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Array{Float64} = []
    for line in eachline(f)
        split_line::Vector{String} = split(line, " ")
        push!(t,parse(Float64, split_line[3]))
    end
    close(f)
    if length(t) == 1
        return t[1]
    else
        return t
    end 
end
function save_time(nn::String,algo::String,tt::Float64)::nothing
    mkpath("times/" * nn * "/")
    f = open("times/" * nn *"/"*algo*".txt","a")
    write(f, string(round(tt; digits=4)))
    close(f)
end


function save_xi(nn::String,algo::String,xi::Float64)::nothing
    mkpath("xis/" * nn * "/")
    f = open("xis/" * nn *"/"*algo*".txt","a")
    write(f, string(round(xi; digits=4)))
    close(f)
end

function save_sample_size(nn::String,algo::String,ss::Int64)::nothing
    mkpath("sample_sizes/" * nn * "/")
    f = open("sample_sizes/" * nn *"/"*algo*".txt","a")
    write(f, string(ss))
    close(f)
end

function save_results(nn::String, cn::String, c::Array{Float64}, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/" * cn * ".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, string(t))
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, "-1.0\n")
        close(f)
    end
end

function save_results_sampling(nn::String, cn::String,c::Array{Float64}, ss::Int64, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/"*cn*"_"*string(ss)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(ss)*".txt", "w")
        write(f, string(round(t; digits=4)) * " " * string(ss) * "\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(ss)*".txt", "w")
        write(f, "-1.0 -1.0 -1.0 -1.0")
        close(f)
    end
end

function save_results_samplings(nn::String, cn::String,c::Array{Float64}, ss::Int64, t::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        append_centrality_values("scores/" * nn * "/"*cn*"_"*string(ss)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(ss)*".txt", "a")
        write(f, string(round(t; digits=4)) * " " * string(ss) * "\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(ss)*".txt", "a")
        write(f, "-1.0 -1.0 -1.0 -1.0")
        close(f)
    end
end


function save_results_progressive_sampling(nn::String, cn::String,c::Array{Float64}, ss::Int64, t::Float64,starting_ss::Int64,xi::Float64 = -1.0)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        append_centrality_values("scores/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "a")
        write(f, string(round(t; digits=4)) * " " * string(ss) *" "*string(xi) *"\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "a")
        write(f, "-1.0 -1.0 -1.0 -1.0,-1.0")
        close(f)
    end
end

function save_results_topk(nn::String,cn::String,c::Array{Tuple{Int64,Float64}},k::Int64,ss::Int64,t::Float64,starting_ss::Int64,xi::Float64 = -1.0)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        append_centrality_values_topk("scores/" * nn * "/"*cn*"_"*string(starting_ss)*"_top_"*string(k)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*"_top_k.txt", "a")
        write(f, string(round(t; digits=4)) * " " * string(ss) *" "*string(xi) *" "*string(k)*"\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*"_top_k.txt", "a")
        write(f, "-1.0 -1.0 -1.0 -1.0,-1.0,-1.0")
        close(f)
    end
end
function save_results_diameter(nn::String,diameter::Int64,vertex_diameter::Int64,avg_distance::Float64,number_of_pairs::Float64,eff_diam::Float64,zeta::Float64,tt::Float64,path_optimality::String,diam_path_size::Float64 = 0.0,eff_diam_path_size::Float64 = 0.0,avg_path_size::Float64 = 0.0,total_couples_path_size::Float64 = 0.0)
    mkpath("diameter/" * nn * "/")
    f::IOStream = open("diameter/" * nn * "/"*path_optimality*".txt", "a")
    if (path_optimality != "pfm") && (!occursin("pfm",path_optimality))
        write(f, string(diameter)*","*string(vertex_diameter)*","*string(avg_distance)*","*string(number_of_pairs)*","*string(eff_diam)*","*string(zeta)*","*string(tt)* "\n")
    else
        write(f, string(diameter)*","*string(vertex_diameter)*","*string(avg_distance)*","*string(number_of_pairs)*","*string(eff_diam)*","*string(zeta)*","*string(tt)*","*string(diam_path_size)*","*string(eff_diam_path_size)*","*string(avg_path_size)*","*string(total_couples_path_size)* "\n")
    end
    close(f)
end

function read_vertex_diameter(nn::String,path_optimality::String)
    file_name = "diameter/"*nn*"/"*path_optimality*".txt"
    @assert isfile(file_name) "The vertex diameter does not exists"
    f = open(file_name, "r")
    t::Int64 = 0
    t = parse(Int64, split(readline(f),",")[2])
    close(f)
    return t  
end


function read_distance_measures(nn::String,path_optimality::String,apx::Bool = false,seed::String = " ")
    
    res = []
    if !apx
        file_name = "diameter/"*nn*"/"*path_optimality*".txt"
        @assert isfile(file_name) "The vertex diameter does not exists"
        f = open(file_name, "r")
        res = [parse(Float64,x) for x in split(readline(f),",")]
        close(f)
    else
        file_name = "diameter/"*nn*"/apx_"*path_optimality*"_"*seed*".txt"
        @assert isfile(file_name) "The vertex diameter does not exists"
        f = open(file_name, "r")
        for line in eachline(f)
            push!(res,[parse(Float64,x) for x in split(line,",")])
        end
        close(f)
    end
    return res
end



function print_algorithm_status(algo::String,bound::String,progressive::Bool)
    println("=======================================")
    println("****************Status*****************")
    println("Algorithm : "*algo)
    println("Using "*bound*" bound ")
    if progressive
        println("Running progressive sampling version")
    else
        println("Running fixed sample size version")
    end
    println("=======================================")
    flush(stdout)
end


function compute_vapnik_chervonenkis_bound(vd::Int64,epsilon::Float64,delta::Float64,c::Float64 = 0.5)::Int64
    return trunc(Int,(c/epsilon^2) * ((floor(log2(vd-2)))+1+log(1/delta)))
end

function compute_hoeffding_bound(n::Int64,epsilon::Float64,delta::Float64)::Int64
    return trunc(Int,(1.0/(2*epsilon^2))*log2(2*n/delta))
end

function compute_starting_sample(epsilon::Float64,delta::Float64)::Int64
    return trunc(Int,((1+8*epsilon+sqrt(1+16*epsilon))*log(6/delta)/(4*epsilon^2)))+1
end

# ================================
#       |Exact Algorithms|
# ================================

function shortest_temporal_betweenness(tg::temporal_graph,verbose_step::Int64, bigint::Bool)
    if nthreads() > 1
        return threaded_temporal_shortest_betweenness(tg,verbose_step,bigint)
    else
        return temporal_shortest_betweenness(tg,verbose_step,bigint)
    end
end


function shortest_foremost_temporal_betweenness(tg::temporal_graph,verbose_step::Int64, bigint::Bool)
    if nthreads() > 1
        return threaded_temporal_shortest_foremost_betweenness(tg,verbose_step,bigint)
    else
        return temporal_shortest_foremost_betweenness(tg,verbose_step,bigint)
    end
end


function prefix_foremost_temporal_betweenness(tg::temporal_graph,verbose_step::Int64)
    if nthreads() > 1
        return threaded_temporal_prefix_foremost_betweenness(tg,verbose_step)
    else
        return temporal_prefix_foremost_betweenness(tg,verbose_step)
    end
end



# ================================
#       |Progressive Sampling|
# ================================



# TO CHANGE
#-------| KADABRA's progressive sampling |-------

function progressive_wub(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64 = 0,bigint::Bool = false,algo::String = "trk",topt::String = "sh",vc_upper_bund::Bool = true,geo::Float64 = 1.2,diam::Int64 = -1,start_factor::Int64 = 100,force_gc::Bool = false)
    @assert (topt == "sh") || (topt == "sfm") || (topt == "pfm") "Illegal temporal-path optimality, use: sh for shortest , sfm for shortest foremost , or pfm for prefix foremost"
    if nthreads() >= 1
        println("Algorithm "*algo* " Temporal path optimality "*topt)
        flush(stdout)
        if k > 0
            if topt == "sh"
                return threaded_progressive_wub_topk(tg,eps,delta,k,bigint,algo,vc_upper_bund,diam,geo,start_factor,force_gc)
            elseif topt == "sfm"
                return threaded_progressive_wub_shortest_foremost_topk(tg,eps,delta,k,bigint,algo,vc_upper_bund,diam,geo,start_factor,force_gc)
            else
                return threaded_progressive_wub_prefix_foremost_topk(tg,eps,delta,k,algo,vc_upper_bund,diam,geo,start_factor,force_gc)
            end
        else
            if topt == "sh"
                return threaded_progressive_wub(tg,eps,delta,bigint,algo,vc_upper_bund,diam,geo,start_factor)
            elseif topt == "sfm"
                return threaded_progressive_wub_shortest_foremost(tg,eps,delta,bigint,algo,vc_upper_bund,diam,geo,start_factor)
            else
                return threaded_progressive_wub_prefix_foremost(tg,eps,delta,algo,vc_upper_bund,diam,geo,start_factor)
            end
        end
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end
end


#-------| Bernstein |-------

# ONBRA


function progressive_bernstein(tg::temporal_graph,epsilon::Float64,delta::Float64,geo::Float64, bigint::Bool = false,algo::String = "trk",topt::String = "sh",vc_upper_bound::Bool = true,force_gc::Bool = false)
    @assert (topt == "sh") || (topt == "sfm") || (topt == "pfm") "Illegal temporal-path optimality, use: sh for shortest , sfm for shortest foremost , or pfm for prefix foremost"
    if nthreads() > 1
        println("Algorithm "*algo* " Temporal path optimality "*topt)
        flush(stdout)
        
        if topt == "sh"
            return  threaded_progressive_bernstein(tg,epsilon,delta,bigint, algo,vc_upper_bound,-1,geo,100,force_gc)
        elseif topt == "sfm"
            return threaded_progressive_bernstein_shortest_foremost(tg,epsilon,delta,bigint, algo,vc_upper_bound,-1,geo,100,force_gc)
        else
            return threaded_progressive_onbra_prefix_foremost_bernstein(tg,epsilon,delta, algo,vc_upper_bound,-1,geo,100,force_gc)
        end
        
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end
end

function progressive_mantra(tg::temporal_graph,eps::Float64,delta::Float64,bigint::Bool,algo::String = "trk",topt::String = "sh",vc_upper_bound::Bool = false,geo::Float64 = 1.2,diam::Int64 = -1,empirical_peeling_a::Float64 = 2.0,force_gc::Bool=false)
    @assert (topt == "sh") || (topt == "sfm") || (topt == "pfm") "Illegal temporal-path optimality, use: sh for shortest , sfm for shortest foremost , or pfm for prefix foremost"
    if nthreads() > 1
        println("Algorithm "*algo* " Temporal path optimality "*topt)
        flush(stdout)
        
        if topt == "sh"
            return  threaded_progressive_cmcera(tg,eps,delta,bigint,algo,vc_upper_bound,diam,empirical_peeling_a,force_gc)
        elseif topt == "sfm"
            return threaded_progressive_cmcera_shortest_foremost(tg,eps,delta,bigint,algo,vc_upper_bound,diam,empirical_peeling_a,force_gc)
        else
            return threaded_progressive_cmcera_prefix_foremost(tg,eps,delta,algo,vc_upper_bound,diam,empirical_peeling_a,force_gc)
        end

    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end
end


function temporal_distance_based_metrics(tg::temporal_graph,topt::String = "sh",sample_size::Int64 = 0,verbose_step::Int64 = 0,threshold::Float64 = 0.9)
    @assert (topt == "sh") || (topt == "sfm") || (topt == "pfm") "Illegal temporal-path optimality, use: sh for shortest , sfm for shortest foremost , or pfm for prefix foremost"
    if nthreads() > 1
        println("Computing temporal distance-based metrics | Temporal path optimality "*topt)
        flush(stdout)
        
        if topt == "sh"
            return  threaded_temporal_shortest_diameter(tg,sample_size,verbose_step,threshold)
        elseif topt == "sfm"
            return threaded_temporal_shortest_foremost_diameter(tg,sample_size,verbose_step,threshold)
        else
            return threaded_temporal_prefix_foremost_diameter(tg,sample_size,verbose_step,threshold)
        end

    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end

end