include("src/APXTBC.jl")
using JSON
function clean_gc()
    GC.gc()
end
path = "graphs/"
big_int = false
f = JSON.parsefile("exp.json")
gn = f["ds"]
topt = f["topt"]
nn = String(split(gn, ".t")[1])
tg = load_temporal_graph(path*gn," ")
print_stats(tg, graph_name= gn)
flush(stdout)

if topt == "pfm"
    result = threaded_temporal_prefix_foremost_betweenness(tg,1000)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])

elseif topt == "sh"
    result = threaded_temporal_shortest_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
else
    result =threaded_temporal_shortest_foremost_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
end
clean_gc()
