function read_centrality_values(file_name::String)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist"
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



function read_time(file_name::String)::Float64
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Float64 = 0.0
    t = parse(Float64, split(readline(f)," ")[1])
    close(f)
    return t  
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



function save_results_diameter(nn::String,diameter::Int64,vertex_diameter::Int64,avg_distance::Float64,number_of_pairs::Float64,eff_diam::Float64,zeta::Float64,tt::Float64,path_optimality::String)
    mkpath("diameter/" * nn * "/")
    f::IOStream = open("diameter/" * nn * "/"*path_optimality*".txt", "a")
    write(f, string(diameter)*","*string(vertex_diameter)*","*string(avg_distance)*","*string(number_of_pairs)*","*string(eff_diam)*","*string(zeta)*","*string(tt)* "\n")
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
    return trunc(Int,(c/epsilon^2) * ((floor(log2(vd-2)))+log(1/delta)))
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
        return threaded_temporal_shortest_betweennes(tg,verbose_step,bigint)
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
#       |Fixed Sample Size|
# ================================

# Random Temporal Betweenness

function rtb_shortest_temporal_betweenness(tg::temporal_graph, sample_size::Int64,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1

        return threaded_rtb(tg,sample_size,verbose_step,bigint)

    else
        return rtb(tg,sample_size,verbose_step,bigint)
    end
end


function rtb_shortest_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64,verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1

        return threaded_rtb_shortest_foremost(tg,sample_size,verbose_step,bigint)
    
    else
        return rtb_shortest_foremost(tg,sample_size,verbose_step,bigint)
    end
end

function rtb_prefix_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64,verbose_step::Int64; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1

        return threaded_rtb_prefix_foremost(tg,sample_size,verbose_step)
    
    else
        return rtb_prefix_foremost(tg,sample_size,verbose_step)
    end

end

# ONBRA

function onbra_shortest_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_onbra(tg,sample_size,verbose_step,bigint)
    else
        return onbra(tg,sample_size,verbose_step,bigint)
    end
end


function onbra_shortest_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_onbra_shortest_foremost(tg,sample_size,verbose_step,bigint)
    else
        return onbra_shortest_foremost(tg,sample_size,verbose_step,bigint)
    end
end

function onbra_prefix_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_onbra_prefix_foremost(tg,sample_size,verbose_step)
    else
        return onbra_prefix_foremost(tg,sample_size,verbose_step)
    end
end

# Temporal R.K.

function trk_shortest_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_trk(tg,sample_size,verbose_step,bigint)
    else
        return trk(tg,sample_size,verbose_step,bigint)
    end
end


function trk_shortest_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64, bigint::Bool; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_trk_shortest_foremost(tg,sample_size,verbose_step,bigint)
    else
        return trk_shortest_foremost(tg,sample_size,verbose_step,bigint)
    end
end

function trk_prefix_foremost_temporal_betweenness(tg::temporal_graph, sample_size::Int64, verbose_step::Int64; test_sample=Array{Tuple{Int64,Int64}}[])
    if nthreads() > 1
        return threaded_trk_prefix_foremost(tg,sample_size,verbose_step)
    else
        return trk_prefix_foremost(tg,sample_size,verbose_step)
    end
end


# ================================
#       |Progressive Sampling|
# ================================

#-------| Rademacher Averages |-------

# ONBRA

function progressive_onbra_shortest_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads() > 1
        return threaded_progressive_onbra(tg,initial_sample,epsilon,delta,geo,verbose_step,bigint)
    else
        return progressive_onbra(tg,initial_sample,epsilon,delta,geo,verbose_step,bigint)
    end
end

function progressive_onbra_shortest_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads() > 1
        return threaded_progressive_onbra_shortest_foremost(tg,initial_sample,epsilon,delta,geo,verbose_step,bigint)
    else
        return progressive_onbra_shortest_foremost(tg,initial_sample,epsilon,delta,geo,verbose_step,bigint)
    end
end


function progressive_onbra_prefix_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)
    if nthreads() > 1
        return threaded_progressive_onbra_prefix_foremost(tg,initial_sample,epsilon,delta,geo,verbose_step)
    else
        return progressive_onbra_prefix_foremost(tg,initial_sample,epsilon,delta,geo,verbose_step)
    end
end

#-------| KADABRA's progressive sampling |-------

function progressive_trk_shortest_temporal_betweenness(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64 = 0,verbose_step::Int64 = 0,bigint::Bool = false,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    if nthreads()>1
        if k > 0
            return threaded_progressive_trk_topk(tg,eps,delta,k,verbose_step,bigint,algo,diam,start_factor,sample_step,hb)
        else
            return threaded_progressive_trk(tg,eps,delta,verbose_step,bigint,algo,diam,start_factor,sample_step,hb)
        end
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end
end


function progressive_trk_shortest_foremost_temporal_betweenness(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64 = 0,verbose_step::Int64 = 0,bigint::Bool = false,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)    
    if nthreads()>1
        if k > 0
            return threaded_progressive_trk_shortest_foremost_topk(tg,eps,delta,k,verbose_step,bigint,algo,diam,start_factor,sample_step,hb)   
        else
            return threaded_progressive_trk_shortest_foremost(tg,eps,delta,verbose_step,bigint,algo,diam,start_factor,sample_step,hb)      
        end
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end 
end


function progressive_trk_prefix_foremost_temporal_betweenness(tg::temporal_graph,eps::Float64,delta::Float64,k::Int64 = 0,verbose_step::Int64= 0,algo::String = "trk",diam::Int64 = -1,start_factor::Int64 = 100,sample_step::Int64 = 10,hb::Bool = false)
    if nthreads()>1
        if k > 0
            return threaded_progressive_trk_prefix_foremost_topk(tg,eps,delta,k,verbose_step,algo,diam,start_factor,sample_step,hb)
        else
            return threaded_progressive_trk_prefix_foremost(tg,eps,delta,verbose_step,algo,diam,start_factor,sample_step,hb)
        end
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end 
end


#-------| Bernstein |-------

# ONBRA

function progressive_onbra_bernstein_shortest_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_onbra_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end 
end


function progressive_onbra_bernstein_shortest_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_onbra_shortest_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end 
end

function progressive_onbra_bernstein_prefix_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)
    if nthreads()>1
        return threaded_progressive_onbra_prefix_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step )
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end 
end

# TRK

function progressive_trk_bernstein_shortest_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_trk_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end


function progressive_trk_bernstein_shortest_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_trk_shortest_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end

function progressive_trk_bernstein_prefix_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)
    if nthreads()>1
        return threaded_progressive_trk_prefix_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step )
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end

# RTB 

function progressive_rtb_bernstein_shortest_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_rtb_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end


function progressive_rtb_bernstein_shortest_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64, bigint::Bool)
    if nthreads()>1
        return threaded_progressive_rtb_shortest_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step, bigint)
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end

function progressive_rtb_bernstein_prefix_foremost_temporal_betweenness(tg::temporal_graph,initial_sample::Int64,epsilon::Float64,delta::Float64,geo::Float64,verbose_step::Int64)
    if nthreads()>1
        return threaded_progressive_rtb_prefix_foremost_bernstein(tg,initial_sample,epsilon,delta,geo,verbose_step )
    else
        println("Error: set the number of threads > 1 julia --threads <thread_number>")
    end  
end
#=

if nthreads() > 1
    return
else
    return
end

=#