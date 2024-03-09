include("src/MANTRA.jl")
path = "graphs/"





datasets = [
    "01_college_msg.txt",
    "02_digg_reply.txt",
    "03_slashdot_reply.txt",
    "04_facebook_wall.txt",
    "05_topology.txt",
    "06_bordeaux.txt",
    "07_mathoverflow.txt",
    "08_SMS.txt",
    "09_askubuntu.txt",
    "10_superuser.txt",
    "11_wiki_talk.txt"
]
bint = false
println("Computing Ground Truth values for the (*)-temporal betweenness ")
println("Suggestion : Go and grab a coffee ;)")
flush(stdout)


println("Computing Ground Truth values for the prefix-foremost temporal betweenness")
for gn in datasets
    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_prefix_foremost_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])
    save_results(nn,"pfm",result[1],result[2])
    clean_gc()   
end
println("Computing Ground Truth values for the shortest temporal betweenness")
for gn in datasets
    tg = load_temporal_graph(path*gn," ")
    if gn == "06_bordeaux.txt"
        bint = true
    else
        bint = false
    end
    result = threaded_temporal_shortest_betweenness(tg,1000,bint,false)
    nn = String(split(gn, ".t")[1])
    save_results(nn,"sh",result[1],result[2])
    clean_gc()
end
println("Computing Ground Truth values for the shoretst-foremost temporal betweenness")
for gn in datasets
    tg = load_temporal_graph(path*gn," ")
    if gn == "06_bordeaux.txt"
        bint = true
    else
        bint = false
    end
    result =threaded_temporal_shortest_foremost_betweenness(tg,1000,bint,false)
    nn = String(split(gn, ".t")[1])
    save_results(nn,"sfm",result[1],result[2])
    clean_gc()
end