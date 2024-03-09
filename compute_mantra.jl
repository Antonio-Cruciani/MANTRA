include("src/MANTRA.jl")
path = "graphs/"
function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
#= test
datasets = ["00_workplace.txt"]
epsilon_list = [0.1,0.07]
sample_list = [100,350]
=#
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
delta =0.1
trials = 10
global big_int = false
geo = 1.2
algo = "ob"
println("Computing absolute (ε,δ)-Approximation of the (*)-temporal betweenness using MANTRA")
println("Suggestion : Go and grab a (quick) coffee ;)")
flush(stdout)
println("Computing (ε,δ)-Approximation for the prefix-foremost temporal betweenness")
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
        println("Running Progressive MANTRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta, big_int,algo,topt)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end

println("Computing (ε,δ)-Approximation for the shortest temporal betweenness")
flush(stdout)
topt = "sh"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        if gn == "06_bordeaux.txt"
            global big_int = true
        else
            global big_int = false
        end
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Progressive MANTRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta, big_int,algo,topt)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end

println("Computing (ε,δ)-Approximation for the shortest-foremost temporal betweenness")
flush(stdout)
topt = "sfm"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        if gn == "06_bordeaux.txt"
            global big_int = true
        else
            global big_int = false
        end
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Progressive MANTRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta, big_int,algo,topt)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end