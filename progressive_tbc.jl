include("src/APXTBC.jl")

function print_samplig_stats(epsilon,delta,trials,ss)
    println(" ε = "*string(epsilon)*" δ = "*string(delta)*" #trials = "*string(trials)*" starting sample size/ub sample size "*string(ss))
    flush(stdout)
end
datasets = ["01_hypertext.txt",
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
"17_brain_100206_70.txt"]

path = "graphs/"

epsilon = 0.05
delta = 0.1
trials = 10
ss = 750
geo = 1.2
for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_onbra_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_onbra_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_prefix_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0)
        save_results_progressive_sampling(nn,"b_trk_pfm",result[1],result[2][end],result[4],ss,result[3])
    end
    #=
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
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_onbra_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_onbra_sh",result[1],result[2][end],result[4],ss,result[3])
    end
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_trk_sh",result[1],result[2][end],result[4],ss,result[3])
    end
    #=
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_shortest_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_rtb_sh",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    print_samplig_stats(epsilon,delta,trials,ss)
    print_stats(tg, graph_name= gn)
    println("Running Bernstein ONBRA")
    flush(stdout)
    for i in 1:trials
        result = progressive_onbra_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_onbra_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    println("Running Bernstein TRK")
    flush(stdout)
    for i in 1:trials
        result = progressive_trk_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_trk_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    #=
    println("Running Bernstein RTB")
    flush(stdout)
    for i in 1:trials
        result = progressive_rtb_bernstein_shortest_foremost_temporal_betweenness(tg,ss,epsilon,delta,geo,0,false)
        save_results_progressive_sampling(nn,"b_rtb_sfm",result[1],result[2][end],result[4],ss,result[3])
    end
    =#
end


#=

for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    vd = read_vertex_diameter(nn,"sh")
    vc = compute_vapnik_chervonenkis_bound(vd,epsilon,delta)
    print_stats(tg, graph_name= gn)
    println("Running ONBRA")
    flush(stdout)
    for i in 1:trials
        result = onbra_shortest_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"onbra_sh",result[1],vc,result[2])
    end
    println("Running TRK")
    flush(stdout)
    for i in 1:trials
        result = trk_shortest_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"trk_sh",result[1],vc,result[2])
    end
    println("Running RTB")
    flush(stdout)
    for i in 1:trials
        result = rtb_shortest_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"rtb_sh",result[1],vc,result[2])
    end
end


for gn in datasets
    nn = String(split(gn, ".t")[1])
    tg = load_temporal_graph(path*gn," ")
    vd = read_vertex_diameter(nn,"sfm")
    vc = compute_vapnik_chervonenkis_bound(vd,epsilon,delta)
    print_stats(tg, graph_name= gn)
    println("Running ONBRA")
    flush(stdout)
    for i in 1:trials
        result = onbra_shortest_foremost_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"onbra_sfm",result[1],vc,result[2])
    end
    println("Running TRK")
    flush(stdout)
    for i in 1:trials
        result = trk_shortest_foremost_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"trk_sfm",result[1],vc,result[2])
    end
    println("Running RTB")
    flush(stdout)
    for i in 1:trials
        result = rtb_shortest_foremost_temporal_betweenness(tg,vc,0,false)
        save_results_samplings(nn,"rtb_sfm",result[1],vc,result[2])
    end
end

=#