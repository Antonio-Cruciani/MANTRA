#@everywhere include("src/APXTBC.jl")
include("src/APXTBC.jl")
path = "graphs/"



#=
function clean_gc()
    GC.gc()
end

=#
println("Computing Groun Truth values for the shortest-foremost temporal betweenness")
## "26_delicius.txt",
## "24_flickr_grow.txt"

datasets = [
"25_stackoverflow.txt"

]



for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
   
    result = threaded_temporal_prefix_foremost_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])
    clean_gc()
    
   #=
    result = threaded_temporal_shortest_betweenness(tg,1000,false,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
    clean_gc()
   
    
    result =threaded_temporal_shortest_foremost_betweenness(tg,1000,false,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
    clean_gc()
    =#
    
end
