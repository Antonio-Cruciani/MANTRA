include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
function clean_gc()
    GC.gc()
end
# PFM


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




trials = 5



datasets = [
"16_brain_100206_90.txt",
"17_brain_100206_70.txt",
"01_hypertext.txt",
"02_highschool.txt",
"03_hospital_ward.txt",
"05_wiki_elections.txt",
"06_highschool.txt",
"08_infectious.txt",
"09_primary_school.txt",
"12_highschool.txt",
"04_college_msg.txt",
"10_facebook_wall.txt",
"11_slashdot_reply.txt",
"07_digg_reply.txt",
"13_topology.txt",
"14_SMS.txt",
"21_mathoverflow.txt",
"20_askubuntu.txt",
"22_superuser.txt",
"18_venice.txt",
"19_bordeaux.txt"
]


epsilon_list = [0.007,0.005,0.001]
sample_list = [1350,1500,2000]
results = []
# pfm
for i in 1:lastindex(epsilon_list)
    eps = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(eps,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        for _ in 1:trials
            println("Approximating (pfm)-Temporal Diameter ")
            _,_,_,_,_,diam_v,_,_,_,_,t_diam = threaded_temporal_prefix_foremost_diameter(tg,64,0,0.9)
            diam_v = diam_v -1
            diam = trunc(Int,diam_v)
            println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam))
            flush(stdout)
            omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
            push!(results,[gn,eps,"pfm",omega])
        end
    end
end


# write results
f::IOStream = open("samplesizes.txt", "a")
for u in 1:lastindex(results)
    for i in 1:lastindex(results[u])
        if i<length(results[u])
            write(f, string(results[u][i]) * " ")
        else
            write(f, string(results[u][i]) * "\n")
        end
    end
end
close(f)
# sh
results = []

for i in 1:lastindex(epsilon_list)
    eps = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(eps,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        for _ in 1:trials
            println("Approximating (sh)-Temporal Diameter ")
            diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_diameter(tg,64,0,0.9,false)
            println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
            flush(stdout)
            omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
            push!(results,[gn,eps,"sh",omega])
        end
    end
end

# write results
f::IOStream = open("samplesizes.txt", "a")
for u in 1:lastindex(results)
    for i in 1:lastindex(results[u])
        if i<length(results[u])
            write(f, string(results[u][i]) * " ")
        else
            write(f, string(results[u][i]) * "\n")
        end
    end
end
close(f)
# sfm



results = []

for i in 1:lastindex(epsilon_list)
    eps = epsilon_list[i]
    starting_ss = sample_list[i]
    for gn in datasets
        nn = String(split(gn, ".t")[1])
        tg = load_temporal_graph(path*gn," ")
        print_samplig_stats(eps,delta,trials,starting_ss)
        print_stats(tg, graph_name= gn)
        for _ in 1:trials
            println("Approximating (sfm)-Temporal Diameter ")
            diam,avg_dist,_,_,_,t_diam = threaded_temporal_shortest_foremost_diameter(tg,64,0,0.9)
            println("Task completed in "*string(round(t_diam;digits = 4))*". Δ = "*string(diam)*" ρ = "*string(avg_dist))
            flush(stdout)
            omega = 0.5/eps/eps * (log2(diam-1)+1+log(2/delta))
            push!(results,[gn,eps,"sfm",omega])
        end
    end
end

# write results
f::IOStream = open("samplesizes.txt", "a")
for u in 1:lastindex(results)
    for i in 1:lastindex(results[u])
        if i<length(results[u])
            write(f, string(results[u][i]) * " ")
        else
            write(f, string(results[u][i]) * "\n")
        end
    end
end
close(f)