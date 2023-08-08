include("src/APXTBC.jl")
path = "graphs/"


datasets = ["20_askubuntu.txt","22_superuser.txt","23_flickr_grow.txt"]


println("Computing Ground Truth values for the prefix-foremost temporal betweenness")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_prefix_foremost_betweenness(tg,1000)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])
end
#=

println("Computing Groun Truth values for the shortest temporal betweenness")

datasets = ["21_mathoverflow.txt"]

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
end

println("Computing Groun Truth values for the shortest-foremost temporal betweenness")


for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_foremost_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
end
=#
#=
roads = ["18_venice.txt","19_bordeaux.txt"]

println("Computing Ground Truth values for the prefix-foremost temporal betweenness")

for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_prefix_foremost_betweenness(tg,0)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])
end


println("Computing Groun Truth values for the shortest temporal betweenness")


for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_betweenness(tg,0,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
end

println("Computing Groun Truth values for the shortest-foremost temporal betweenness")


for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_foremost_betweenness(tg,0,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
end

=#