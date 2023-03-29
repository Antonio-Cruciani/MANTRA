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


function save_results_progressive_sampling(nn::String, cn::String,c::Array{Float64}, ss::Int64, t::Float64,starting_ss::Int64,xi::Float64 = -1.0)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "w")
        write(f, string(round(t; digits=4)) * " " * string(ss) *" "*string(xi) *"\n")
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/"*cn*"_"*string(starting_ss)*".txt", "w")
        write(f, "-1.0 -1.0 -1.0 -1.0,-1.0")
        close(f)
    end
end