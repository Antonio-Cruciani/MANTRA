include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
function clean_gc()
    GC.gc()
end

path = "graphs/"

delta = 0.1
trials = 5
geo = 1.2
big_int = false
vc_upper_bound = false
algo = "rtb"
k = 0
upperbound_sample = "vc"

topt = "pfm"
epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,1350,1500,2000]


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
"22_superuser.txt",
"18_venice.txt",
"19_bordeaux.txt"
]

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
