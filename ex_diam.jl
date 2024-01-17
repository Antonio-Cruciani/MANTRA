include("src/APXTBC.jl")


datasets = ["16_brain_100206_90.txt",
"17_brain_100206_70.txt","01_hypertext.txt",
"02_highschool.txt",
"03_hospital_ward.txt",
"04_college_msg.txt",
"05_wiki_elections.txt",
"06_highschool.txt",
"07_digg_reply.txt",
"08_infectious.txt",
"09_primary_school.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"12_highschool.txt",
"13_topology.txt",
"14_SMS.txt",
"18_venice.txt",
"19_bordeaux.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]

path = "graphs/"
trials = 5
#=
println("Computing Ground Truth values for the prefix-foremost temporal diameter")

seed = 0
for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_prefix_foremost_diameter(tg,seed,1000)
    println(result)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],trunc(Int,result[6]),result[2],result[4],result[3],result[5],result[11],"pfm",result[7],result[8],result[9],result[10])

end

seed = 64
for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_prefix_foremost_diameter(tg,seed,0)
        println(result)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],trunc(Int,result[6]),result[2],result[4],result[3],result[5],result[11],"apx_pfm_"*string(seed),result[7],result[8],result[9],result[10])
    end
end

seed = 256
for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_prefix_foremost_diameter(tg,seed,0)
        println(result)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],trunc(Int,result[6]),result[2],result[4],result[3],result[5],result[11],"apx_pfm_"*string(seed),result[7],result[8],result[9],result[10])
    end
end


seed = 0
println("Computing Ground Truth values for the shortest temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_diameter(tg,seed,1000)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"sh")
end

seed = 64
println("Computing Ground Truth values for the shortest temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_shortest_diameter(tg,seed,0)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"apx_sh_"*string(seed))
    end
end

seed = 256
println("Computing Ground Truth values for the shortest temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_shortest_diameter(tg,seed,0)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"apx_sh_"*string(seed))
    end
end

datasets = [
"22_superuser.txt"
]


datasets = ["16_brain_100206_90.txt",
"17_brain_100206_70.txt","01_hypertext.txt",
"02_highschool.txt",
"03_hospital_ward.txt",
"04_college_msg.txt",
"05_wiki_elections.txt",
"06_highschool.txt",
"07_digg_reply.txt",
"08_infectious.txt",
"09_primary_school.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"12_highschool.txt",
"13_topology.txt",
"14_SMS.txt",
"18_venice.txt",
"19_bordeaux.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]
=#



seed = 64

println("Computing Ground Truth values for the shortest-foremost temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_shortest_foremost_diameter(tg,seed,0)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"apx_sfm_"*string(seed))
    end
end
seed = 256

println("Computing Ground Truth values for the shortest-foremost temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    for i in 1:trials
        result = threaded_temporal_shortest_foremost_diameter(tg,seed,0)
        nn = String(split(gn, ".t")[1])
        
        save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"apx_sfm_"*string(seed))
    end
end


datasets = [
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]

seed = 0
println("Computing Ground Truth values for the shortest-foremost temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_foremost_diameter(tg,seed,1000)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"sfm")
end
