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
vc_upper_bound = false
#topt = "pfm"
algo = "ob"
k = 0
topt = "pfm"
upperbound_sample = "vc"


#=

trials = 1



datasets = [

"03_hospital_ward.txt",
"05_wiki_elections.txt",
"06_highschool.txt",
"08_infectious.txt",
"09_primary_school.txt",
"12_highschool.s for help. Thanks for your effort.
txt"
]



epsilon_list = [0.001]
sample_list = [2000]




topt = "sh"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)

        
        println("Running c-MC ERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,big_int,algo,topt,false,1.2,-1,2.0,true)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
        
        
        
    end
end


topt = "sfm"
for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        println("Running c-MC ERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,big_int,algo,topt,false,1.2,-1,2.0,true)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
        
               
    end
end

=#

trials = 1
epsilon_list = [0.007,0.005,0.001]
sample_list = [1350,1500,2000]
datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]



topt = "sh"
big_int = true

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)

        
        println("Running c-MC ERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,big_int,algo,topt,false,1.2,-1,2.0,true)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
        
        
        
    end
end



topt = "sfm"

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)

        
        println("Running c-MC ERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,big_int,algo,topt,false,1.2,-1,2.0,true)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
            clean_gc()
        end
        
        
        
    end
end