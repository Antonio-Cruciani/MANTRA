include("src/APXTBC.jl")

#=
datasets = ["01_hypertext.txt",
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
"16_brain_100206_90.txt",
"17_brain_100206_70.txt",
"18_venice.txt",
"19_bordeaux.txt"
]
=#
datasets = [
"23_wiki_talk.txt"
]
path = "graphs/"


println("Computing Ground Truth values for the prefix-foremost temporal diameter")


for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_prefix_foremost_diameter(tg,0,1000)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],trunc(Int,result[6]),result[2],result[4],result[3],result[5],result[7],"pfm")

end

#=
println("Computing Ground Truth values for the shortest temporal diameter")

for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_diameter(tg,0,0)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"sh")
end


println("Computing Ground Truth values for the shortest-foremost temporal diameter")


for gn in datasets

    tg = load_temporal_graph(path*gn," ")
    print_stats(tg, graph_name= gn)
    flush(stdout)
    result = threaded_temporal_shortest_foremost_diameter(tg,0,0)
    nn = String(split(gn, ".t")[1])
    
    save_results_diameter(nn,result[1],result[1]+1,result[2],result[4],result[3],result[5],result[6],"sfm")
end

=#