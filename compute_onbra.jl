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
epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,1350,1500,2000]
trials = 10
big_int = false
geo = 1.2
algo = "ob"
vc_upperbound = true
println("Computing absolute (δ,ε)-Approximation of the (*)-temporal betweenness using ONBRA")
println("Suggestion : Go and grab a coffee ;)")
flush(stdout)

println("Computing (δ,ε)-Approximation for the prefix-foremost temporal betweenness")
flush(stdout)
topt = "pfm"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Progressive ONBRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upperbound)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end

println("Computing (δ,ε)-Approximation for the shortest temporal betweenness")
flush(stdout)
topt = "sh"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        if gn == "06_bordeaux.txt"
            big_int = true
        else
            big_int = false
        end
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Progressive ONBRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upperbound)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end

println("Computing (δ,ε)-Approximation for the shortest-foremost temporal betweenness")
flush(stdout)
topt = "sfm"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        if gn == "06_bordeaux.txt"
            big_int = true
        else
            big_int = false
        end
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Progressive ONBRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upperbound)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end