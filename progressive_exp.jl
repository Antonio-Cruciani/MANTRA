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
epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,0.007,0.005,0.001]
#epsilon_list = [0.001]
#sample_list = [2000]
#epsilon_list = [0.007]
#sample_list = [1350]
#epsilon_list = [0.005]
#sample_list = [1500]


datasets = [
    "18_venice.txt",
    "19_bordeaux.txt",
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

#=



datasets = [
"21_mathoverflow.txt",
"13_topology.txt"
]

datasets = [
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",
"23_wiki_talk.txt"
]

datasets = [
"23_wiki_talk.txt"
]

datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]
=#
#datasets = ["21_mathoverflow.txt","22_superuser.txt"]
#=
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

trials = 1

datasets = [
"04_college_msg.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"07_digg_reply.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]

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
trials = 5
#=
datasets = [
"04_college_msg.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"07_digg_reply.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",
"18_venice.txt",
"19_bordeaux.txt"
]

"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",

datasets = [
"18_venice.txt",
"19_bordeaux.txt"
]
=#
datasets = [
    "24_flickr_grow.txt",
    "25_stackoverflow.txt",
    "26_delicius.txt"
    ]

topt = "pfm"
big_int = false
epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,1350,1500,2000]
#epsilon_list = [0.007,0.001]
#sample_list = [1350,2000]
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
#=
trials = 1

epsilon_list = [0.01]
sample_list = [1000]


datasets = [
    "22_superuser.txt"

]

topt = "sfm"
big_int = false

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




trials = 1
epsilon_list = [0.001]
sample_list = [2000]
datasets = [
    "13_topology.txt"
]



topt = "sfm"
big_int = false

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


epsilon_list = [0.007,0.005,0.001]
sample_list = [1350,1500,2000]

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
epsilon_list = [0.001]
sample_list = [2000]
trials = 4
datasets = [
"04_college_msg.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"07_digg_reply.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]

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
epsilon_list = [0.007,0.005,0.001]
sample_list = [1350,1500,2000]
datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]


trials = 4
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


=#