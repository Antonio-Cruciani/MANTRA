
using Conda
using PyCall
Conda.add("scipy")

py"""
from scipy import stats
def ktau(x, y):
    return stats.kendalltau(x, y)
def wktau(x, y):
    return stats.weightedtau(x, y)
def spear(x,y):
    return stats.spearmanr(x, y)
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
    return sp, kt, wkt
end

function normalize_centrality(x::Array{Float64})::Array{Float64}
    n::Int64 = length(x)
    for i in 1:n
        x[i] = x[i]/(n*(n-1))
    end
    return x
end
function correlations_fixed_ss(method::String,ss_array::Array{String},datasets::Array{String},algorithms::Array{String} = ["rtb","ob","trk"])
    i= 1
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        println("Loaded")
        ss = ss_array[i]
        i+=1
        for algo in algorithms
            println("Reading : "*"scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            println("Loaded")
            corr = compute_correlations(exact,apx_cc,false)
            push!(results,[graph,method,ss,algo,corr[1][1],corr[2][1],corr[3][1]])
        end
    end
    mkpath("correlations/")

    open("correlations/fixed_ss.txt","a") do file
        write(file,"Graph \t Method \t SampleSize \t Algorithm\t Spearman\t KendallTau\t weightedKT\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\t"*string(res[7])*"\n")
        end
    end
end

function correlations_progressive_ss(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})
    i= 1
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        println("Loaded")
        i+=1
        for algo in algo_starting_ss
            println("Reading : "*"scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            println("Loaded")
            corr = compute_correlations(exact,apx_cc,false)
            push!(results,[graph,method,algo[2],algo[1],corr[1][1],corr[2][1],corr[3][1]])
        end
    end
    mkpath("correlations/")

    open("correlations/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tStartingSampleSize\tAlgorithm\tSpearman\tKendallTau\tweightedKT\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\t"*string(res[7])*"\n")
        end
    end
end

function compute_errors(x,y)
    @assert length(x) == length(y) "The two rankings have different number of elements"
    sup = 0.0
    mse = 0.0
    for i in 1:lastindex(x)
        ae = abs(x[i]-y[i])
        mse += (x[i]-y[i])^2
        if ae > sup 
            sup = ae
        end
    end
    mse = mse * 1/length(x)
    return sup,mse

end
function abs_error_fixed_ss(method::String,ss_array::Array{String},datasets::Array{String},algorithms::Array{String} = ["rtb","ob","trk"])
    i = 1
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        ss = ss_array[i]
        i+=1
        for algo in algorithms
            println("Reading : "*"scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            errors = compute_errors(exact,apx_cc)
            push!(results,[graph,method,ss,algo,errors[1],errors[2]])
        end
    end
    mkpath("errors/")

    open("errors/fixed_ss.txt","a") do file
        write(file,"Graph\tMethod\tSampleSize\tAlgorithm\tSD\tMSE\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end
end
function abs_error_progressive_ss(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        for algo in algo_starting_ss
            println("Reading : "*"scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            errors = compute_errors(exact,apx_cc)
            push!(results,[graph,method,algo[2],algo[1],errors[1],errors[2]])
        end
    end
    mkpath("errors/")

    open("errors/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tStartingSampleSize\tAlgorithm\tSD\tMSE\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end
end

function times_fixed_ss(method::String,ss_array::Array{String},datasets::Array{String},algorithms::Array{String} = ["rtb","ob","trk"],pos::Int64 = 1)

    i = 1
    results = []

    for graph in datasets
        println("Reading : "*"times/" * graph * "/exact_"*method*"/time_exact_"*method*"_0.txt")
        exact = read_time("times/" * graph * "/exact_"*method*"/time_exact_"*method*"_0.txt",true)
        ss = ss_array[i]
        i+=1
        for algo in algorithms
            println("Reading : "*"times/" * graph * "/"*algo*"_"*method*"/time_"*algo*"_"*method*"_"*ss*".txt")
            apx_tt = read_time("times/" * graph * "/"*algo*"_"*method*"/time_"*algo*"_"*method*"_"*ss*".txt",false,pos)
            push!(results,[graph,method,ss,algo,exact,apx_tt])
        end
    end

    mkpath("speed/")

    open("speed/fixed_ss.txt","a") do file
        write(file,"Graph\tMethod\tSampleSize\tAlgorithm\tExact\tApproximation\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end

end



function times_progressive_ss(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})

   
    results = []

    for graph in datasets
        println("Reading : "*"times/" * graph * "/exact_"*method*"/time_exact_"*method*"_0.txt")
        exact = read_time("times/" * graph * "/exact_"*method*"/time_exact_"*method*"_0.txt",true)
        for algo in algo_starting_ss
            println("Reading : "*"times/" * graph * "/"*algo[1]*"_"*method*"/time_"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_tt = read_time("times/" * graph * "/"*algo[1]*"_"*method*"/time_"*algo[1]*"_"*method*"_"*algo[2]*".txt",false)
            push!(results,[graph,method,algo[2],algo[1],exact,apx_tt])
        end
    end

    mkpath("speed/")

    open("speed/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tStartingSampleSize\tAlgorithm\tExact\tApproximation\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end

end



function xi_progressive_ss(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})

   
    results = []

    for graph in datasets
        for algo in algo_starting_ss
            println("Reading : "*"epsilons/" * graph * "/"*algo[1]*"_"*method*"/eps_"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_tt = read_time("epsilons/" * graph * "/"*algo[1]*"_"*method*"/eps_"*algo[1]*"_"*method*"_"*algo[2]*".txt",false)
            push!(results,[graph,method,algo[2],algo[1],apx_tt])
        end
    end

    mkpath("empiricalerror/")

    open("empiricalerror/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tStartingSampleSize\tAlgorithm\tXi\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\n")
        end
    end

end

function kappas_progressive_ss(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})

   
    results = []

    for graph in datasets
        for algo in algo_starting_ss
            println("Reading : "*"kappas/" * graph * "/"*algo[1]*"_"*method*"/k_"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_tt = read_time("kappas/" * graph * "/"*algo[1]*"_"*method*"/k_"*algo[1]*"_"*method*"_"*algo[2]*".txt",false)
            push!(results,[graph,method,algo[2],algo[1],apx_tt])
        end
    end

    mkpath("empiricalerror/")

    open("empiricalerror/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tStartingSampleSize\tAlgorithm\tXi\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\n")
        end
    end

end


function intersection_progressive(method::String,algo_starting_ss::Array{Tuple{String,String}},datasets::Array{String})
    i= 1
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        println("Loaded")
        i+=1
        for algo in algo_starting_ss
            println("Reading : "*"scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo[1]*"_"*method*"/"*algo[1]*"_"*method*"_"*algo[2]*".txt")
            println("Loaded")
            inter = intersection(exact,apx_cc,50)
            push!(results,[graph,method,algo[2],algo[1],inter[10],inter[50]])
        end
    end

    mkpath("intersections/")

    open("intersections/progressive_ss.txt","a") do file
        write(file,"Graph\tMethod\tSampleSize\tAlgorithm\tIntersection10\tIntersection50\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end
end


function intersection_fixed_ss(method::String,ss_array::Array{String},datasets::Array{String},algorithms::Array{String} = ["rtb","ob","trk"])
    i= 1
    results = []

    for graph in datasets
        #@assert isfile("scores/" * nn * "/" * cn1 * ".txt") "The values of " * cn1 * " are not available for the network " * nn
        println("Reading : "*"scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt")
        exact = normalize_centrality(read_centrality_values("scores/" * graph * "/exact_"*method*"/exact_"*method*"_0.txt"))
        println("Loaded")
        ss = ss_array[i]
        i+=1
        for algo in algorithms
            println("Reading : "*"scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            apx_cc = read_centrality_values("scores/" * graph * "/"*algo*"_"*method*"/"*algo*"_"*method*"_"*ss*".txt")
            println("Loaded")
            inter = intersection(exact,apx_cc,50)
            push!(results,[graph,method,ss,algo,inter[10],inter[50]])
        end
    end
    mkpath("intersections/")

    open("intersections/fixed_ss.txt","a") do file
        write(file,"Graph\tMethod\tSampleSize\tAlgorithm\tIntersection10\tIntersection50\n")
        for res in results
            write(file,res[1]*"\t"*res[2]*"\t"*res[3]*"\t"*res[4]*"\t"*string(res[5])*"\t"*string(res[6])*"\n")
        end
    end
end
