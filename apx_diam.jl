include("src/APXTBC.jl")
datasets = ["17_brain_100206_70.txt"]
#"23_flickr_grow.txt" 
path = "graphs/"
seeds = 256
println("Computing APX values for the prefix-foremost temporal diameter")
trials = 10
verb = 0
for gn in datasets
    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        start_time = time()
        println("Experiment "*string(i)*"/"*string(trials))
        flush(stdout)
        result = threaded_temporal_prefix_foremost_diameter(tg,seeds,verb)
        nn = String(split(gn, ".t")[1])
        finish_time = string(round(time() - start_time; digits=4))
        println("Iteration "*string(i)*"/"*string(trials)*" completed in "*finish_time*" seconds")
        flush(stdout)
        save_results_diameter(nn,result[1],trunc(Int,result[6]),result[2],result[4],result[3],result[5],result[7],"pfm")
    end
end
