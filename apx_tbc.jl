@everywhere include("src/APXTBC.jl")
path = "graphs/"


datasets = [
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]

#=
println("Computing Ground Truth values for the prefix-foremost temporal betweenness")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_prefix_foremost_betweenness(tg,1000)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])
end
println("Computing Groun Truth values for the shortest temporal betweenness")
=#

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_shortest_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
end

println("Computing Groun Truth values for the shortest-foremost temporal betweenness")

datasets = [
    "16_brain_100206_90.txt",
    "17_brain_100206_70.txt",
    "01_hypertext.txt",
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
    "21_mathoverflow.txt",
    "20_askubuntu.txt",
    "22_superuser.txt"
]

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_shortest_foremost_betweenness(tg,1000,false)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
end


roads = ["18_venice.txt","19_bordeaux.txt"]

println("Computing Ground Truth values for the prefix-foremost temporal betweenness")

for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_prefix_foremost_betweenness(tg,1000)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"pfm",result[1],result[2])
end


println("Computing Groun Truth values for the shortest temporal betweenness")


for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_shortest_betweenness(tg,1000,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sh",result[1],result[2])
end

println("Computing Groun Truth values for the shortest-foremost temporal betweenness")


for gn in roads

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = distributed_temporal_shortest_foremost_betweenness(tg,1000,true)
    nn = String(split(gn, ".t")[1])

    save_results(nn,"sfm",result[1],result[2])
end