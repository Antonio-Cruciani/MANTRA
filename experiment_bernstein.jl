include("src/APXTBC.jl")
using JSON

function clean_gc()
    GC.gc()
end

path = "graphs/"
epsilon_list = [0.1,0.07,0.05,0.01]
sample_list = [100,350,750,1000]
path = "graphs/"
sample_step = 32
delta = 0.1
trials = 5
geo = 1.2
big_int = false
#topt = "pfm"
algo = "ob"
k = 0
topt = "pfm"
epsilon = 0.1
delta = 0.1
f = JSON.parsefile("exp.json")
gn = f["ds"]
epsilon = f["eps"]
starting_ss = f["ss"]
topt = f["topt"]
algo = f["algo"]
sampler = f["sampler"]
nn = String(split(gn, ".t")[1])
tg = load_temporal_graph(path*gn," ")
print_stats(tg, graph_name= gn)
if (sampler == "bernstein")
    println("Running Bernstein")
    flush(stdout)
    result = progressive_bernstein(tg,starting_ss,epsilon,delta,geo,1000, big_int,algo,topt)
    save_results_progressive_sampling(nn,"b_"*algo*"_"*topt,result[1],result[2][end],result[4],starting_ss,result[3])
elseif sampler == "wub"
    println("Running WUB")
    result = progressive_wub(tg,epsilon,delta,k,10000,big_int,algo,topt,-1,100,sample_step)
    save_results_progressive_sampling(nn,"wub_"*algo*"_"*topt,result[1],result[4],result[6],starting_ss,epsilon)
else
    println("Running CMCERA")
    result = progressive_cmcera(tg,epsilon,delta,0,big_int,algo,topt,2.0,sample_step)
    save_results_progressive_sampling(nn,"cm_"*algo*"_"*topt,result[1],result[2],result[4],starting_ss,epsilon)
end
clean_gc()
