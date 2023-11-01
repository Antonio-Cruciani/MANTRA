include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
function clean_gc()
    GC.gc()
end
# PFM
#epsilon_list = [0.005]
#sample_list = [1500]

path = "graphs/"
sample_step = 32
delta = 0.1
trials = 5
geo = 1.2
big_int = false
vc_upper_bound = true
#topt = "pfm"
algo = "ob"
k = 0
topt = "sh"
upperbound_sample = "vc"
#=
epsilon_list = [0.005]
sample_list = [1500]
datasets = [

    "22_superuser.txt"
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
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upper_bound)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,big_int,algo,topt,vc_upper_bound,geo)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Running c-MC ERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
    end
end
=#
#=
upperbound_sample = "rho"
vc_upper_bound = false

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
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upper_bound)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,big_int,algo,topt,vc_upper_bound,geo)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
    end
end



datasets = [

    "22_superuser.txt"
]
epsilon_list = [0.01]
sample_list = [1000]

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
            result = progressive_bernstein(tg,epsilon,delta,geo, big_int,algo,topt,vc_upper_bound)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end

        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,big_int,algo,topt,vc_upper_bound,geo)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[4],result[6],starting_ss,epsilon)
            clean_gc()
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end

    end
end
clean_gc()
=#
# SFM
epsilon_list = [0.01]
sample_list = [1000]
topt = "sfm"

trials = 5

datasets = [
    "22_superuser.txt"
]

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,big_int,algo,topt,vc_upper_bound,geo)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[4],result[6],starting_ss,epsilon)
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


