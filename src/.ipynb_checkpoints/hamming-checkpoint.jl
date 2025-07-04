# hamming_distance_calculations.jl
# Functions to calculate Hamming distance and related metrics.

"""
    ham_dist(vec1::Array{Int8,1}, vec2::Array{Int8,1})

Calculate the Hamming distance between two sequences.

# Arguments
- `vec1::Array{Int8,1}`: First sequence.
- `vec2::Array{Int8,1}`: Second sequence.

# Returns
- `distance::Int`: Hamming distance between `vec1` and `vec2`.
"""
function ham_dist(vec1::Array{Int8,1}, vec2::Array{Int8,1})
    return sum(vec1 .!= vec2)
end

function ham_dist(vec1::Array{Int,1}, vec2::Array{Int,1})
    return sum(vec1 .!= vec2)
end

function ham_dist(vec1::Array{String,1}, vec2::Array{String,1})
    return sum(vec1 .!= vec2)
end

# Assuming all_codons is already defined as a 61-element Vector{String}

# Function to compute Hamming distance between two codons
function cod_dist(c1::String, c2::String)
    sum(c1[i] != c2[i] for i in 1:3)
end

function create_codon_dist_dict()
    all_codons = vcat([amino2cod[i] for i in 1:20]...)
    codon_distance_dict = Dict{Tuple{String, String}, Int}()

    for c1 in all_codons, c2 in all_codons
        codon_distance_dict[(c1, c2)] = codon_distance(c1, c2)
    end
    
    return codon_distance_dict
end

    
function ham_dist_deg(a::Array{String,1}, b::Array{String,1})
    return sum([cod_dist(a[i],b[i]) for i in 1:length(a)])
end

"""
    ham_dist(vec::Array{Int8,1}, msa::Array{Int8,2})

Calculate the Hamming distance between a sequence and each column in a multiple sequence alignment (MSA).

# Arguments
- `vec::Array{Int8,1}`: Sequence to compare.
- `msa::Array{Int8,2}`: MSA, where each column is a sequence.

# Returns
- `distances::Array{Int}`: Array of Hamming distances for each sequence in the MSA.
"""
function ham_dist(vec::Array{Int8,1}, msa::Array{Int8,2})
    return [ham_dist(vec, msa[:, i]) for i in 1:size(msa, 2)]
end

function ham_dist(vec::Array{Int,1}, msa::Array{Int,2})
    return [ham_dist(vec, msa[:, i]) for i in 1:size(msa, 2)]
end

function ham_dist(vec::Array{String,1}, msa::Array{String,2})
    return [ham_dist(vec, msa[:, i]) for i in 1:size(msa, 2)]
end


function ham_dist_deg(vec::Array{String,1}, msa::Array{String,2})
    return [ham_dist_deg(vec, msa[:, i]) for i in 1:size(msa, 2)]
end

"""
    ham_dist(msa1::Array{Int8,2}, msa2::Array{Int8,2})

Calculate pairwise Hamming distances between corresponding columns in two MSAs.

# Arguments
- `msa1::Array{Int8,2}`: First MSA.
- `msa2::Array{Int8,2}`: Second MSA.

# Returns
- `distances::Array{Int}`: Array of Hamming distances for each pair of sequences.
"""
function ham_dist(msa1::Array{Int8,2}, msa2::Array{Int8,2})
    return [ham_dist(msa1[:, i], msa2[:, i]) for i in 1:size(msa1, 2)]
end

function ham_dist(msa1::Array{Int,2}, msa2::Array{Int,2})
    return [ham_dist(msa1[:, i], msa2[:, i]) for i in 1:size(msa1, 2)]
end

function ham_dist_deg(msa1::Array{String,2}, msa2::Array{String,2})
    return [ham_dist_deg(msa1[:, i], msa2[:, i]) for i in 1:size(msa1, 2)]
end

"""
    ham_dist(step_msa::Array{Array{Int8,2},1})

Calculate Hamming distances for each step in a sequence alignment path relative to the initial step.

# Arguments
- `step_msa::Array{Array{Int8,2},1}`: Array of MSAs, each representing a step in the path.

# Returns
- `distances::Array{Float64,2}`: Matrix of Hamming distances for each chain at each step.
"""
function ham_dist(step_msa::Array{Array{Int8,2},1})
    N_steps = length(step_msa)
    N_chains = size(step_msa[1], 2)
    res = zeros(N_steps, N_chains)
    for n in 1:N_steps
        res[n, :] .= ham_dist(step_msa[1], step_msa[n])
    end
    return res
end

function ham_dist(step_msa::Array{Array{Int,2},1})
    N_steps = length(step_msa)
    N_chains = size(step_msa[1], 2)
    res = zeros(N_steps, N_chains)
    for n in 1:N_steps
        res[n, :] .= ham_dist(step_msa[1], step_msa[n])
    end
    return res
end

"""
    ham_dist_AB_rp2(step_msa::Array{Int8,3}, step_msa_B::Array{Int8,3}, num::Int)

Calculate mean and mean squared Hamming distances between two sequence alignments.

# Arguments
- `step_msa::Array{Int8,3}`: First set of MSAs.
- `step_msa_B::Array{Int8,3}`: Second set of MSAs.
- `num::Int`: Number of chains to consider.

# Returns
- `res::Array{Float64}`: Array of mean Hamming distances.
- `res_sq::Array{Float64}`: Array of mean squared Hamming distances.
"""
function ham_dist_AB_rp2(step_msa::Array{Int8,3}, step_msa_B::Array{Int8,3}, num::Int)
    L, N_points, N_chains = size(step_msa)
    res = zeros(N_points)
    res_sq = zeros(N_points)
    for n in 1:N_points
        a, a_sq = [], []
        for r in 1:num
            for p in r+1:num
                hh = ham_dist(step_msa[:, n, r], step_msa_B[:, n, p])
                push!(a, hh)
                push!(a_sq, hh^2)
            end
        end
        res[n] = mean(a)
        res_sq[n] = mean(a_sq)
    end
    return res, res_sq
end

"""
    ham_dist_AB_rp(step_msa::Array{Array{Int,2},1}, step_msa_B::Array{Array{Int,2},1}, num::Int)

Calculate the average Hamming distances between two sets of MSAs across multiple steps.

# Arguments
- `step_msa::Array{Array{Int,2},1}`: First set of MSAs.
- `step_msa_B::Array{Array{Int,2},1}`: Second set of MSAs.
- `num::Int`: Number of chains to compare.

# Returns
- `res::Array{Float64}`: Mean Hamming distance for each step.
- `res_sq::Array{Float64}`: Mean squared Hamming distance for each step.
"""
function ham_dist_AB_rp(step_msa::Array{Array{Int,2},1}, step_msa_B::Array{Array{Int,2},1}, num::Int)
    N_steps = length(step_msa)
    res, res_sq = zeros(N_steps), zeros(N_steps)
    for n in 1:N_steps
        a, a_sq = [], []
        for r in 1:num
            for p in r+1:num
                hh = ham_dist(step_msa[n][:, r], step_msa_B[n][:, p])
                push!(a, hh)
                push!(a_sq, hh^2)
            end
        end
        res[n] = mean(a)
        res_sq[n] = mean(a_sq)
    end
    return res, res_sq
end

"""
    pairwise_ham_dist(msa::Array{Int8,2}; n_seq = 100, all = false)

Calculate the pairwise Hamming distances for a given MSA.

# Arguments
- `msa::Array{Int8,2}`: MSA where each column represents a sequence.
- `n_seq::Int`: Number of sequences to include in the calculation (default is 100).
- `all::Bool`: If `true`, returns all pairwise distances; otherwise, returns the mean (default is `false`).

# Returns
- `mean_dist::Float64` or `distances::Array{Float64}`: Mean pairwise Hamming distance if `all=false`, otherwise an array of pairwise distances.
"""
function pairwise_ham_dist(msa::Array{Int8,2}; n_seq = 100, all = false)
    L,M = size(msa)
    new_msa = msa[:,sample(1:M,n_seq, replace = false)]
    res = []
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res, ham_dist(new_msa[:, i], new_msa[:, j]))
        end
    end
    return all ? res : mean(res)
end


function pairwise_ham_dist2(msa::Array{Int8,2}; n_seq = 100, all = false)
    L, M = size(msa)
    new_msa = msa[:, sample(1:M, n_seq, replace=false)]

    # Create all pairs (i, j) where i < j
    pairs = [(i, j) for i in 1:n_seq for j in i+1:n_seq]

    # Compute distances in parallel using @tasks
    res = @tasks for (i, j) in pairs
        ham_dist(new_msa[:, i], new_msa[:, j])
    end

    return all ? res : mean(res)
end


function pairwise_ham_dist3(msa::Array{Int8,2}; n_seq = 100, all = false)
    L, M = size(msa)
    new_msa = msa[:, sample(1:M, n_seq, replace=false)]
    println("new")

    # Generate all unique pairs (i, j)
    pairs = [(i, j) for i in 1:n_seq for j in i+1:n_seq]

    # Chunk size: number of pairs per thread
    n_threads = Threads.nthreads()
    chunk_size = cld(length(pairs), n_threads)

    # Split into chunks
    chunks = [pairs[i:min(i+chunk_size-1, end)] for i in 1:chunk_size:length(pairs)]

    # Compute Hamming distances in parallel chunks
    partial_results = @tasks for chunk in chunks
        map(p -> ham_dist(new_msa[:, p[1]], new_msa[:, p[2]]), chunk)
    end

    # Flatten results
    res = reduce(vcat, partial_results)

    return all ? res : mean(res)
end





function pairwise_ham_dist(msa::Array{Int,2}; n_seq = 100, all = false)
    L,M = size(msa)
    new_msa = msa[:,sample(1:M,n_seq, replace = false)]

    res = []
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res, ham_dist(new_msa[:, i], new_msa[:, j]))
        end
    end
    return all ? res : mean(res)
end

function pairwise_ham_dist(msa::Array{String,2}; n_seq = 100, all = false)
    L,M = size(msa)
    new_msa = msa[:,sample(1:M,n_seq, replace = false)]

    res = []
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res, ham_dist(new_msa[:, i], new_msa[:, j]))
        end
    end
    return all ? res : mean(res)
end

function pairwise_ham_dist_deg(msa::Array{String,2}; n_seq = 100, all = false)
    L,M = size(msa)
    new_msa = msa[:,sample(1:M,n_seq, replace = false)]

    res = []
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res, ham_dist_deg(new_msa[:, i], new_msa[:, j]))
        end
    end
    return all ? res : mean(res)
end





function ham_dist_inter_msa(msa1::Array{T,2}, msa2::Array{T,2}; shuffles = 10) where {T}
    res = []
    for _ in 1:shuffles
        m1 = msa1[:,shuffle(1:size(msa1, 2))]
        m2 = msa2[:,shuffle(1:size(msa2, 2))]
        dd = ham_dist(m1,m2)
        append!(res, dd)
    end

    return res
end

