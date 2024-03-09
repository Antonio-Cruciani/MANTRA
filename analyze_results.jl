include("src/MANTRA.jl")
#STATS 
include("src/statistics/correlations_and_error.jl")
#=test
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


path_opt = ["pfm","sh","sfm"]

algo = "ob"

prog_samplers = ["cm","b"]

trials = 10


kappas = [10,25,50]

j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for p in path_opt
        for s in ["cm"]
            for ub in ["var"]
                println("ANALIZYING EPS: "*string(eps)*" PATH OPT "*string(p)*" PROG SAMPLER "*s* " VC")
                get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
                get_errors(p,ss,eps,datasets,algo,s,trials,ub)
                get_times(p,ss,eps,datasets,algo,s,ub)
                for kappa in kappas
                    get_ranking_intersections(p,kappa,ss,eps,datasets,algo,s,trials,ub)
                end
            end
        end
        
    end
end
global j = 1
for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for p in path_opt
        for s in ["b"]
            for ub in ["vc"]
                println("ANALIZYING EPS: "*string(eps)*" PATH OPT "*string(p)*" PROG SAMPLER "*s* " VC")
                get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
                get_errors(p,ss,eps,datasets,algo,s,trials,ub)
                get_times(p,ss,eps,datasets,algo,s,ub)
                for kappa in kappas
                    get_ranking_intersections(p,kappa,ss,eps,datasets,algo,s,trials,ub)
                end
            end
        end
        
    end
end
