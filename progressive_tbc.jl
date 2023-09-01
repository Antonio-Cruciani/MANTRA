include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
#"20_askubuntu.txt",

datasets = [
"22_superuser.txt"
]

path = "graphs/"

epsilon = 0.01
delta = 0.1
trials = 5
ss = 1000
geo = 1.2
big_int = false

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_bernstein(tg,ss,epsilon,delta,geo,10000,false,"ob","pfm")
        save_results_progressive_sampling(nn,"b_onbra_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    #=
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_trk_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_rtb_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end

#=

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_onbra_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,1000,big_int)
        save_results_progressive_sampling(nn,"b_onbra_sh",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
#=
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_trk_sh",result[1],result[2][end],result[4],ss,result[3])
    end
    
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_rtb_sh",result[1],result[2][end],result[4],ss,result[3])
    end
   =# 
#end
#=
datasets = ["18_venice.txt","19_bordeaux.txt"]

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_onbra_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,big_int)
        save_results_progressive_sampling(nn,"b_onbra_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    #=
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_trk_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_rtb_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end


=#