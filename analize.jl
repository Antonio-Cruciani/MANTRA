include("src/APXTBC.jl")

# mathoverflow 0.007 must be done
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

]

path_opt = ["sh"]

algo = "ob"

#prog_samplers = ["b","wub","cm"]
prog_samplers = ["cm"]

trials = 5

#upper_bound_sampless = ["vc","var"]
upper_bound_sampless = ["vc"]

kappas = [10,25,50]

epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,1350,1500,2000]
#epsilon_list = [0.007]
#sample_list = [1350]
j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for p in path_opt
        for s in ["cm"]
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
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",

]

path_opt = ["sfm"]

epsilon_list = [0.1,0.07,0.05,0.01,0.007,0.005,0.001]
sample_list = [100,350,750,1000,1350,1500,2000]
#epsilon_list = [0.007]
#sample_list = [1350]
j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for p in path_opt
        for s in ["cm"]
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

datasets = [
"18_venice.txt",
"19_bordeaux.txt",
"03_hospital_ward.txt",
"04_college_msg.txt",
"07_digg_reply.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"13_topology.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",

]


datasets = [
"18_venice.txt",
"19_bordeaux.txt",
"04_college_msg.txt",
"05_wiki_elections.txt",
"07_digg_reply.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"13_topology.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",

]

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
"14_SMS.txt",
"21_mathoverflow.txt",

]
trials = 5
#=
datasets = [

"16_brain_100206_90.txt",
"17_brain_100206_70.txt",
"01_hypertext.txt",
"02_highschool.txt",
"03_hospital_ward.txt",
"06_highschool.txt",
"08_infectious.txt",
"09_primary_school.txt",
"12_highschool.txt",

]


path_opt = ["pfm"]

epsilon_list = [0.001]
sample_list = [2000]
#epsilon_list = [0.007]
#sample_list = [1350]
j = 1

for eps in epsilon_list
    ss = sample_list[j]
    global j += 1
    for p in path_opt
        for s in ["cm"]
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
=#

# FIX THE VAR AND VC STUFF in files
#=
datasets = [
    "18_venice.txt",
    "19_bordeaux.txt"
]
trials = 1
j = 1
for eps in epsilon_list
    ss = sample_list[j]
    global j+=1
    for p in path_opt
        for s in ["b","wub","cm"]
            for ub in ["vc"]
                get_correlations(p,ss,eps,datasets,algo,s,trials,ub)
                get_errors(p,ss,eps,datasets,algo,s,trials,ub)
                get_times(p,ss,eps,datasets,algo,s,ub)
                for kappa in kappas
                    get_ranking_intersections(p,kappa,ss,eps,datasets,algo,s,trials,ub)
                end
            end
        end
        for s in ["b","wub","cm"]
            for ub in ["var"]
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

=#