
using DataStructures
using StatsBase
using Base.Threads
using Distributed

include("graphs/temporal_graph.jl")
# Shortest
include("centralities/shortest_temporal/temporal_shortest_betweenness.jl")
include("centralities/shortest_temporal/onbra.jl")
include("centralities/shortest_temporal/rtb.jl")
include("centralities/shortest_temporal/trk.jl")

# Shortest Foremost
include("centralities/shortest_foremost/temporal_shortest_foremost_betweenness.jl")
include("centralities/shortest_foremost/onbra_shortest_foremost.jl")
include("centralities/shortest_foremost/rtb_shortest_foremost.jl")
include("centralities/shortest_foremost/trk_shortest_foremost.jl")

# Prefix Foremost
include("centralities/prefix_foremost/prefix_foremost_betweenness.jl")
include("centralities/prefix_foremost/onbra_prefix_foremost.jl")
include("centralities/prefix_foremost/rtb_prefix_foremost.jl")
include("centralities/prefix_foremost/trk_prefix_foremost.jl")
# Distance function
include("temporal_neighborhood_function/mt_diameter.jl")
include("centralities/utilities.jl")
include("statistics/rankings.jl")
#include("statistics/correlations_and_error.jl")

#PARALLEL 
# EXACT
include("parallel/shortest_temporal/mt_temporal_shortest_betweenness.jl")
include("parallel/shortest_foremost/mt_temporal_shortest_foremost_betweenness.jl")
include("parallel/prefix_foremost/mt_prefix_foremost_betweenness.jl")


include("parallel/progressive/shortest_temporal/mt_bernstein.jl")
include("parallel/progressive/shortest_temporal/mt_cmcera.jl")
include("parallel/progressive/shortest_temporal/mt_weighted_ub.jl")

include("parallel/progressive/shortest_foremost/mt_bernstein_shortest_foremost.jl")
include("parallel/progressive/shortest_foremost/mt_cmcera_shortest_foremost.jl")
include("parallel/progressive/shortest_foremost/mt_weighted_ub_shortest_foremost.jl")

include("parallel/progressive/prefix_foremost/mt_bernstein_prefix_foremost.jl")
include("parallel/progressive/prefix_foremost/mt_cmcera_prefix_foremost.jl")
include("parallel/progressive/prefix_foremost/mt_weighted_ub_prefix_foremost.jl")


#STATS 
#include("statistics/correlations_and_error.jl")

# SOME PRINTS 
packet_name::String =raw"
         __  __          _   _ _______ _____            
        |  \/  |   /\   | \ | |__   __|  __ \     /\    
        | \  / |  /  \  |  \| |  | |  | |__) |   /  \   
        | |\/| | / /\ \ | . ` |  | |  |  _  /   / /\ \  
        | |  | |/ ____ \| |\  |  | |  | | \ \  / ____ \ 
        |_|  |_/_/    \_\_| \_|  |_|  |_|  \_\/_/    \_\
       "

println(packet_name)
flush(stdout)