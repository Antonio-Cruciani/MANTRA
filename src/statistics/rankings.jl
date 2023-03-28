function min_h_k(a::Array{Float64}, b::Array{Float64}, k::Int64)::Int64
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (k > length(a))
        k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    min_h_b_k_a::Int64 = 0
    for a_j in 1:k
        b_j::Int64 = 1
        while (bi[b_j] != ai[a_j])
            b_j = b_j + 1
        end
        if (b_j > min_h_b_k_a)
            min_h_b_k_a = b_j
        end
    end
    return min_h_b_k_a
end

function intersection(a::Array{Float64}, b::Array{Float64}, max_k::Int64)::Array{Float64}
    @assert length(a) == length(b) "The two rankings have different number of elements"
    if (max_k > length(a))
        max_k = length(a)
    end
    ai::Vector{Int64} = sortperm(a, rev=true)
    bi::Vector{Int64} = sortperm(b, rev=true)
    inter::Array{Float64} = zeros(max_k)
    for k in 1:max_k
        inter[k] = length(intersect(Set(ai[1:k]), Set(bi[1:k])))
    end
    return inter
end

