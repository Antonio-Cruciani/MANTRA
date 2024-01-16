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
trials = 1
geo = 1.2
big_int = false
vc_upper_bound = true
#topt = "pfm"
algo = "ob"
k = 0
topt = "sh"
upperbound_sample = "vc"

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
"18_venice.txt",
"19_bordeaux.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"
]
trials = 5
epsilon_list = [0.01,0.007,0.005,0.001]
sample_list = [1000,1350,1500,2000]

topt = "pfm"

for i in 1:lastindex(epsilon_list)
    epsilon = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(epsilon,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)

        
        println("Running ONBRA")
        flush(stdout)
        for i in 1:trials
            result = progressive_bernstein(tg,epsilon,delta,big_int,algo,topt,true,true)
            save_results_progressive_sampling(nn,"b_"*algo*"_"*topt*"_"*upperbound_sample,result[1],result[2][end],result[4],starting_ss,result[3])
            clean_gc()
        end
        
        
        
    end
end
