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

function read_onbra_centrality_values(file_name::String, ss::Int64, nn::Int64)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist"
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = zeros(nn)
    value::Float64 = 0.0
    for _ in 1:ss
        for n in 1:nn
            value = parse(Float64, readline(f))
            if (value < -0.1)
                println("ERROR. There are negative values with absolute big values")
                return Array{Float64}([])
            end
            if (value < 0)
                value = 0
            end
            centrality[n] += value / ss
        end
    end
    close(f)
    return centrality
end


function read_time(file_name::String,exact::Bool,pos::Int64 = 1)::Float64
    @assert isfile(file_name) "The time value file does not exist"
    f = open(file_name, "r")
    t::Float64 = 0.0
    if exact
        t = parse(Float64, split(readline(f)," ")[2])
    else
        t = parse(Float64, split(readline(f)," ")[pos])
    end
    close(f)
    return t
    
end