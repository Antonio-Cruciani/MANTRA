
#using Conda
using PyCall
#Conda.add("scipy")

py"""
from scipy import stats
def ktau(x, y):
    return stats.kendalltau(x, y)
def wktau(x, y):
    return stats.weightedtau(x, y)
def spear(x,y):
    return stats.spearmanr(x, y)
def pear(x,y):
    return stats.pearsonr(x, y)
"""

function compute_correlations(x, y, verbose::Bool)
    sp = py"spear"(x, y)
    if (verbose)
        log("    Spearman computed")
    end
    kt = py"ktau"(x, y)
    if (verbose)
        log("    Kendall tau computed")
    end
    wkt = py"wktau"(x, y)
    if (verbose)
        log("    Weighted Kendall tau computed")
    end
    p =  py"pear"(x, y)
    if (verbose)
        log("    Pearson computed")
    end
    return sp, kt, wkt,p
end

function normalize_centrality(x::Array{Float64})::Array{Float64}
    n::Int64 = length(x)
    for i in 1:n
        x[i] = x[i]/(n*(n-1))
    end
    return x
end

function get_correlations(method::String,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="ob",prog_sampler::String = "wub",trials::Int64 = 5,upper_bound_samples::String = "vc")
    results = []
    
    for graph in datasets
        tg = load_temporal_graph("graphs/"*graph, " ")
        gn = split(graph,".txt")[1]
        #println("Nnodes "*string(tg.num_nodes))
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/"*method*".txt"))
        #if prog_sampler == "cm"
        apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*string(starting_sample)*".txt")
        #else
        #    apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt")
        #end
        apx_bc = zeros(tg.num_nodes)

        for k in 1:trials
            for i in 1:tg.num_nodes
                apx_bc[i] += apx_cc[i+(tg.num_nodes*(k-1))]
            end
        end
        apx_bc = apx_bc .* [1/trials]
       
        for i in 1:tg.num_nodes
            if apx_bc[i]> 1
                apx_bc[i] = 0.0
            end
        end
        
        corr = compute_correlations(exact,apx_bc,false)
        push!(results,[gn,method,algo,prog_sampler,corr[1],corr[2],corr[3],corr[4]])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/correlations.txt")
        header = true
    end
    open("analysis/correlations.txt","a") do file
        if header
            write(file,"Graph,Method,Epsilon,Algorithm,Sampler,Spearman,KendallTau,weightedKT,Pearson,UpperBoundSampleSize\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*prog_sampler*","*res[3]*","*string(res[5][1])*","*string(res[6][1])*","*string(res[7][1])*","*string(res[8][1])*","*upper_bound_samples*"\n")
        end
    end
end

function supremum_deviation(x,y)
    SD = 0
    for i in 1:lastindex(x)
        if abs(x[i]-y[i]) > SD
            SD = abs(x[i]-y[i]) 
        end
    end
    return SD
end

function mean_sqaured_error(x,y)
    MSE = 0
    for i in 1:lastindex(x)
        MSE += (x[i]-y[i])^2
    end
    return (1/lastindex(x)) * MSE
end

function get_max_temporal_bc(method::String,datasets::Array{String})
    results = []
    
    for graph in datasets
        gn = split(graph,".txt")[1]
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/"*method*".txt"))
        push!(results,[gn,method,maximum(exact)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/max_tbc.txt")
        header = true
    end
    open("analysis/max_tbc.txt","a") do file
        if header
            write(file,"Graph,Method,MaxTBC\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(res[3])*"\n")
        end
    end
end

function get_errors(method::String,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="ob",prog_sampler::String = "wub",trials::Int64 = 5,upper_bound_samples::String = "vc")
    results = []
    
    for graph in datasets
        tg = load_temporal_graph("graphs/"*graph, " ")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/"*method*".txt"))
        #if prog_sampler == "cm"
        apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*string(starting_sample)*".txt")
        #else
        #apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt")
        #end
        apx_bc = zeros(tg.num_nodes)
        mse_list = []
        sd_list = []
        for k in 1:trials
            tmp_scores = []
            for i in 1:tg.num_nodes
                push!(tmp_scores,apx_cc[i+(tg.num_nodes*(k-1))])
            end
            push!(mse_list,mean_sqaured_error(exact,tmp_scores))
            push!(sd_list,supremum_deviation(exact,tmp_scores))
        end
        for k in 1:trials
            for i in 1:tg.num_nodes
                apx_bc[i] += apx_cc[i+(tg.num_nodes*(k-1))]
            end
        end
        apx_bc = apx_bc .* [1/trials]
        for i in 1:tg.num_nodes
            if apx_bc[i]> 1
                apx_bc[i] = 0.0
            end
        end
        SD = supremum_deviation(exact,apx_bc)
        MSE = mean_sqaured_error(exact,apx_bc)
        push!(results,[gn,method,algo,prog_sampler,SD,MSE,mean(sd_list),std(sd_list),mean(mse_list),std(mse_list)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/errors.txt")
        header = true
    end
    open("analysis/errors.txt","a") do file
        if header
            write(file,"Graph,Method,Epsilon,Algorithm,Sampler,SD,MSE,SDavg,SDstd,MSEavg,MSEstd,UpperBoundSampleSize\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*prog_sampler*","*res[3]*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*","*string(res[10])*","*upper_bound_samples*"\n")
        end
    end
end

function get_times(method::String,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="ob",prog_sampler::String = "wub",upper_bound_samples::String = "vc")
    results = []
    
    for graph in datasets
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = read_time("times/"*gn*"/time_"*method*".txt")
        #if prog_sampler == "cm"
        apx_path = "times/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*string(starting_sample)*".txt"
        #else
        #apx_path = "times/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt"
        #end
        times = read_time(apx_path)
        samples =read_sample_size(apx_path)
        xi = read_xi(apx_path)
        push!(results,[gn,method,algo,prog_sampler,exact,mean(times),std(times),mean(samples),std(samples),mean(xi),std(xi)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/times.txt")
        header = true
    end
    open("analysis/times.txt","a") do file
        if header
            write(file,"Graph,Method,Epsilon,Algorithm,Sampler,ExactTime,ApxTimeMean,ApxTimeStd,SampleSizeMean,SampleSizeStd,XiMean,XiStd,UpperBoundSampleSize\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*res[3]*","*prog_sampler*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*","*string(res[10])*","*string(res[11])*","*upper_bound_samples*"\n")
        end
    end
end

function get_stats_topk(method::String,tk::Int64,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="ob",prog_sampler::String = "wub",trials::Int64 = 5,upper_bound_samples::String = "vc",ntasks::Int64 = 10)
    results = []
    for graph in datasets
        tg = load_temporal_graph("graphs/"*graph, " ")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/"*method*".txt"))
        apx_path = gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*"_top_"*string(tk)*".txt"
        apx_path_t = gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*"_top_k.txt"
        us = min(tg.num_nodes, max(trunc(Int,sqrt(lastindex(tg.temporal_edges))/ntasks),tk+20))
        apx_cc = read_centrality_values_topk("scores/"*apx_path)
        times = read_time("times/"*apx_path_t)
        samples =read_sample_size("times/"*apx_path_t)
        top_k_exact::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
        for u in 1:tg.num_nodes
            push!(top_k_exact,(u,exact[u]))
        end
        sort!(top_k_exact, by=top_k_exact->(-top_k_exact[2],top_k_exact[1]))
        exact_rank = [u[1] for u in top_k_exact]
        intersection_list = []
        for k in 1:trials
            tmp_top_k_apx = []
            for i in 1:tg.num_nodes
                push!(tmp_top_k_apx,(apx_cc[i+(tg.num_nodes*(k-1))][1],apx_cc[i+(tg.num_nodes*(k-1))][2]))
            end
            tmp_top_k_rank = [u[1] for u in tmp_top_k_apx]
            push!(intersection_list,length(intersect(Set(exact_rank[1:tk]), Set(tmp_top_k_rank[1:us]))))
        end
        push!(results,[gn,method,algo,tk,mean(intersection_list),std(intersection_list),mean(times),std(times),mean(samples),std(samples)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/intersections_topk.txt")
        header = true
    end
    open("analysis/intersections_topk.txt","a") do file
        if header
            write(file,"Graph,Method,Epsilon,Algorithm,Sampler,k,IntersectionMean,IntersectionStd,ApxTimeMean,ApxTimeStd,SampleSizeMean,SampleSizeStd,UpperBoundSampleSize\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*res[3]*","*prog_sampler*","*string(res[4])*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*","*string(res[10])*","*upper_bound_samples*"\n")
        end
    end

end
function get_ranking_intersections(method::String,tk::Int64,starting_sample::Int64,target_epsilon::Float64,datasets::Array{String},algo::String ="ob",prog_sampler::String = "wub",trials::Int64 = 5,upper_bound_samples::String = "vc")
    results = []
    for graph in datasets
        tg = load_temporal_graph("graphs/"*graph, " ")
        gn = split(graph,".txt")[1]
        #println("Analyzing "*gn*" ε ",target_epsilon, " Algorithm ",algo," Progressive Sampler ",prog_sampler, " Path Optimality ",method, " Upper bound samples via ",upper_bound_samples)
        exact = normalize_centrality(read_centrality_values("scores/"*gn*"/"*method*".txt"))
        #if prog_sampler == "cm"
        apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*string(starting_sample)*".txt")
        #else
        #apx_cc = read_centrality_values("scores/"*gn*"/"*prog_sampler*"_"*algo*"_"*method*"_"*upper_bound_samples*"_"*string(starting_sample)*".txt")
        #end
        top_k_exact::Array{Tuple{Int64,Float64}} = Array{Tuple{Int64,Float64}}([])
        for u in 1:tg.num_nodes
            push!(top_k_exact,(u,exact[u]))
        end
        sort!(top_k_exact, by=top_k_exact->(-top_k_exact[2],top_k_exact[1]))
        exact_rank = [u[1] for u in top_k_exact]
        apx_bc = zeros(tg.num_nodes)
        intersection_list = []
        min_h_list = []
        jaccard_list = []
        for k in 1:trials
            tmp_top_k_apx = []
            tmp_tbc = Array{Float64}([])
            for i in 1:tg.num_nodes
                push!(tmp_top_k_apx,(i,apx_cc[i+(tg.num_nodes*(k-1))]))
                push!(tmp_tbc,apx_cc[i+(tg.num_nodes*(k-1))])
            end
            sort!(tmp_top_k_apx, by=tmp_top_k_apx->(-tmp_top_k_apx[2],tmp_top_k_apx[1]))
            tmp_top_k_rank = [u[1] for u in tmp_top_k_apx]
            push!(intersection_list,length(intersect(Set(exact_rank[1:tk]), Set(tmp_top_k_rank[1:tk]))))
            push!(min_h_list,min_h_k(exact,tmp_tbc,tk))
            push!(jaccard_list,jaccard(exact,tmp_tbc,tk)[end])
        end
        for k in 1:trials
            for i in 1:tg.num_nodes
                apx_bc[i] += apx_cc[i+(tg.num_nodes*(k-1))]
            end
        end
        apx_bc = apx_bc .* [1/trials]
        tmp_top_k_apx = []
        for i in 1:tg.num_nodes
            push!(tmp_top_k_apx,(i,apx_bc[i]))
        end
        sort!(tmp_top_k_apx, by=tmp_top_k_apx->(-tmp_top_k_apx[2],tmp_top_k_apx[1]))
        tmp_top_k_rank = [u[1] for u in tmp_top_k_apx]
        push!(results,[gn,method,algo,tk,length(intersect(Set(exact_rank[1:tk]), Set(tmp_top_k_rank[1:tk]))),mean(intersection_list),std(intersection_list),min_h_k(exact,apx_bc,tk),mean(min_h_list),std(min_h_list),jaccard(exact,apx_bc,tk)[end],mean(jaccard_list),std(jaccard_list)])
    end
    mkpath("analysis/")
    header = false
    if !isfile("analysis/intersections.txt")
        header = true
    end
    open("analysis/intersections.txt","a") do file
        if header
            write(file,"Graph,Method,Epsilon,Algorithm,Sampler,k,Intersection,IntersectionMean,IntersectionStd,h,hMean,hStd,Jaccard,JaccardMean,JaccardStd,UpperBoundSampleSize\n")            
        end
        for res in results
            write(file,res[1]*","*res[2]*","*string(target_epsilon)*","*res[3]*","*prog_sampler*","*string(res[4])*","*string(res[5])*","*string(res[6])*","*string(res[7])*","*string(res[8])*","*string(res[9])*","*string(res[10])*","*string(res[11])*","*string(res[12])*","*string(res[13])*","*upper_bound_samples*"\n")
        end
    end
end



function compute_temporal_distance_residuals(datasets,path_optimality)
 
    all_results = []
    for graph in datasets
        gn = split(graph,".txt")[1]
        exact = read_distance_measures(gn,path_optimality)
        apx = read_distance_measures(gn,path_optimality,true)
        res = compute_residuals(exact,apx,path_optimality)
        # write the procedure
        push!(all_results,res)
    end
    return all_results
end









