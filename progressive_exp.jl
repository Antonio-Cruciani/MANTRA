include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end

# PFM

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
"19_bordeaux.txt"
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt"

]

#=
datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]
=#
path = "graphs/"
sample_step = 96
epsilon_list = [0.05,0.01,0.005]
#epsilon = 0.05
delta = 0.1
trials = 10
sample_list = [750,1000,1500]
#starting_ss = 750
geo = 1.2
big_int = false
topt = "pfm"
algo = "ob"
k = 0
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
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
        end

    end
end

#SH

topt = "sh"
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
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
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
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
        end

    end
end


# SFM

big_int = false

topt = "sfm"
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
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
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
        end
        println("Running WUB")
        flush(stdout)
        for i in 1:trials
            result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
            save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
        end
        println("Runnin CMCERA")
        flush(stdout)
        for i in 1:trials
            result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
            save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
        end

    end
end

big_int = false
