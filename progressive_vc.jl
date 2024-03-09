include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
datasets = [
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
"16_brain_100206_90.txt",
"17_brain_100206_70.txt",
"18_venice.txt",
"19_bordeaux.txt",
"20_askubuntu.txt",
"21_mathoverflow.txt",
"22_superuser.txt"
]

path = "graphs/"

epsilon = 0.005
delta = 0.1
trials = 5
k = 0
sample_step = 32

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    #vd = read_vertex_diameter(nn,"sh")
    #vc = compute_vapnik_chervonenkis_bound(vd,epsilon,delta)
    #print_samplig_stats(epsilon,delta,trials,vc)
    print_stats(tg, graph_name= gn)
    println("Running W.UB. ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_prefix_foremost_temporal_betweenness(tg,epsilon,delta,k,10000,"ob",-1,100,sample_step)
        save_results_progressive_sampling(nn,"wub_onbra_pfm",result[1],result[4],result[6],result[5],epsilon)
    end
    #=
    println("Running W.UB. TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_prefix_foremost_temporal_betweenness(tg,epsilon,delta,k,0,"trk",vd)
        save_results_progressive_sampling(nn,"wub_trk_pfm",result[1],result[4],result[6],vc,epsilon)
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
    vd = read_vertex_diameter(nn,"sh")
    vc = compute_vapnik_chervonenkis_bound(vd,epsilon,delta)
    print_samplig_stats(epsilon,delta,trials,vc)
    print_stats(tg, graph_name= gn)
    println("Running W.UB. ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_shortest_temporal_betweenness(tg,epsilon,delta,k,0,false,"ob",vd)
        save_results_progressive_sampling(nn,"wub_onbra_sh",result[1],result[4],result[6],vc,epsilon)
    end
    #=
    println("Running W.UB. TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_shortest_temporal_betweenness(tg,epsilon,delta,k,0,false,"trk",vd)
        save_results_progressive_sampling(nn,"wub_trk_sh",result[1],result[4],result[6],vc,epsilon)
    end
    
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_rtb_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end



for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    vd = read_vertex_diameter(nn,"sh")
    vc = compute_vapnik_chervonenkis_bound(vd,epsilon,delta)
    print_samplig_stats(epsilon,delta,trials,vc)
    print_stats(tg, graph_name= gn)
    println("Running W.UB. ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_shortest_foremost_temporal_betweenness(tg,epsilon,delta,k,0,false,"ob",vd)
        save_results_progressive_sampling(nn,"wub_onbra_sfm",result[1],result[4],result[6],vc,epsilon)
    end
    #=
    println("Running W.UB. TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_shortest_foremost_temporal_betweenness(tg,epsilon,delta,k,0,false,"trk",vd)
        save_results_progressive_sampling(nn,"wub_trk_sfm",result[1],result[4],result[6],vc,epsilon)
    end
    
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_rtb_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end
=#
