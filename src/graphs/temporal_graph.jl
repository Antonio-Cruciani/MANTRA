struct temporal_graph
    num_nodes::Int64
    temporal_edges::Array{Tuple{Int64,Int64,Int64}}
    file_id::Array{String}
    file_time::Array{Int64}
end


struct static_graph
    num_nodes::Int64
    edges::Array{Tuple{Int64,Int64}}
end

 
function load_temporal_graph(file_name::String, sep::String)
    @assert isfile(file_name) "The temporal edge list file " * file_name * " does not exist"
    current_node_id::Int64 = 1
    file_id::Vector{String} = []
    file_id_to_graph_id::Dict{String,Int64} = Dict{String,Int64}()
    current_time::Int64 = 1
    file_time::Vector{Int64} = []
    file_time_to_graph_time::Dict{Int64,Int64} = Dict{Int64,Int64}()
    temporal_edges::Set{Tuple{Int64,Int64,Int64}} = Set{Tuple{Int64,Int64,Int64}}()
    t::Int64 = 0
    f::IOStream = open(file_name, "r")
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        @assert length(split_line) == 3 "Bad line format: " * line
        t = parse(Int64, split_line[3])
        if split_line[1] != split_line[2]
            if (!haskey(file_id_to_graph_id, split_line[1]))
                file_id_to_graph_id[split_line[1]] = current_node_id
                push!(file_id, split_line[1])
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_id_to_graph_id, split_line[2]))
                file_id_to_graph_id[split_line[2]] = current_node_id
                push!(file_id, split_line[2])
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_time_to_graph_time, t))
                file_time_to_graph_time[t] = current_time
                push!(file_time, t)
                current_time = current_time + 1
            end
        end
    end
    close(f)
    sort!(file_time)
    for t in 1:lastindex(file_time)
        file_time_to_graph_time[file_time[t]] = t
    end
    f = open(file_name, "r")
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        if split_line[1] != split_line[2]
            t = parse(Int64, split_line[3])
            push!(temporal_edges, (file_id_to_graph_id[split_line[1]], file_id_to_graph_id[split_line[2]], file_time_to_graph_time[t]))
        end
    end
    temporal_edges_array::Array{Tuple{Int64,Int64,Int64}} = sort!(collect(temporal_edges), by=te -> te[3])
    return temporal_graph(length(file_id_to_graph_id), temporal_edges_array, file_id, file_time)
end

function print_stats(tg::temporal_graph; graph_name="anonymous")
    println("====================================================")
    println("Temporal network: " * graph_name)
    println("====================================================")
    println("Number of nodes " * string(tg.num_nodes))
    println("Number temporal of edges " * string(length(tg.temporal_edges)))
    println("Number of unique time stamps " * string(length(tg.file_time)))
    println("====================================================")
end

function temporal_adjacency_list(tg::temporal_graph)::Array{Array{Tuple{Int64,Int64}}}
    tal::Array{Array{Tuple{Int64,Int64}}} = Array{Array{Tuple{Int64,Int64}}}(undef, tg.num_nodes)
    for u in 1:tg.num_nodes
        tal[u] = Tuple{Int64,Int64,Int64}[]
    end
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        push!(tal[te[1]], (te[2], te[3]))
    end
    return tal
end

function temporal_incidency_list(tg::temporal_graph)::Array{Array{Tuple{Int64,Int64}}}
    tal::Array{Array{Tuple{Int64,Int64}}} = Array{Array{Tuple{Int64,Int64}}}(undef, tg.num_nodes)
    for u in 1:tg.num_nodes
        tal[u] = Tuple{Int64,Int64,Int64}[]
    end
    te::Tuple{Int64,Int64,Int64} = (0, 0, 0)
    for i in 1:lastindex(tg.temporal_edges)
        te = tg.temporal_edges[i]
        push!(tal[te[2]], (te[1], te[3]))
    end
    return tal
end

function adjacency_list(sg::static_graph)::Array{Array{Int64}}
    al::Array{Int64} = Array{Array{Int64}}(undef, sg.num_nodes)
    for u in 1:sg.num_nodes
        al[u] = Array{Int64}([])
    end
    for i in 1:lastindex(sg.edges)
        push!(al[sg.edges[i][1]], sg.edges[i][2])
    end
    return al
end

function incidency_list(sg::static_graph)::Array{Array{Int64}}
    al::Array{Int64} = Array{Array{Int64}}(undef, sg.num_nodes)
    for u in 1:sg.num_nodes
        al[u] = Array{Int64}([])
    end
    for i in 1:lastindex(sg.edges)
        push!(al[sg.edges[i][2]], sg.edges[i][1])
    end
    return al
end

function underlying_graph(tg::temporal_graph)::static_graph
    edge_set::Set{Tuple{Int64,Int64}} = Set{Tuple{Int64,Int64}}()
    for te in tg.temporal_edges
        push!(edge_set,(te[1],te[2]))
    end  
    return static_graph(tg.num_nodes,sort!(collect(edge_set), by=e -> (e[1],e[2])))
end