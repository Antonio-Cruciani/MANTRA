include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
function clean_gc()
    GC.gc()
end
# PFM

epsilon_list = [0.1,0.07,0.05,0.01]
sample_list = [100,350,750,1000]
path = "graphs/"
sample_step = 32
delta = 0.1
trials = 5
geo = 1.2
big_int = false
#topt = "pfm"
algo = "ob"
k = 0
topt = "pfm"

datasets = [
    "16_brain_100206_90.txt",
    "17_brain_100206_70.txt",

]



for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Bernstein")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end





#SH
topt = "sh"
trials = 5





for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Bernstein")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end

        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end





# SFM

big_int = false
epsilon_list = [0.1,0.07,0.05,0.01]
sample_list = [100,350,750,1000]
topt = "sfm"


for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Bernstein")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end

datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]
big_int = true

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Bernstein")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end



datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]
big_int = true
topt = "sh"

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running Bernstein")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end


