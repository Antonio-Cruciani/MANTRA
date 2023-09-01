include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
#=
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
=#
datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]

path = "graphs/"
sample_step = 32
epsilon = 0.05
delta = 0.1
trials = 5
starting_ss = 750
geo = 1.2
big_int = true
topt = "sh"
algo = "ob"
k = 0
for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,starting_ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein")
    flush(stdout)
    for i in 1:trials
        result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
        save_results_progressive_sampling(nn,"b_onbra_sh",result[1],result[2][end],result[4],starting_ss,result[3])
    end
    println("Running WUB")
    flush(stdout)
    for i in 1:trials
        result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
        save_results_progressive_sampling(nn,"wub_onbra_sh",result[1],result[4],result[6],result[5],epsilon)
    end
    println("Runnin CMCERA")
    flush(stdout)
    for i in 1:trials
        result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
        save_results_progressive_sampling(nn,"sil_onbra_sh",result[1],result[2],result[4],starting_ss,epsilon)
    end

end