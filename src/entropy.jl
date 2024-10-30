# entropy_calculations.jl
# Functions for calculating entropy and related properties.

"""
    get_entropy(f::Array{T,2}; q::Int = 21) where {T}

Compute the entropy at each position in a sequence alignment frequency matrix.

# Arguments
- `f::Array{T,2}`: Frequency matrix where `f[a, i]` is the frequency of state `a` at position `i`.
- `q::Int`: Number of states for each position (default is 21).

# Returns
- `entr::Array{Float64}`: Array of entropies for each position in the sequence.
"""
function get_entropy(f::Array{T,2}; q::Int = 21) where {T}
    N = length(f[1, :])
    entr = zeros(Float64, N)
    for i in 1:N
        for a in 1:q
            if f[a, i] > 0
                entr[i] -= f[a, i] * log(f[a, i])
            end
        end
    end
    return entr / log(2)
end

"""
    single_site_prob_cond(k::Int, mutated_seq::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}, L::Int; q::Int = 21) where {T}

Calculate the conditional probability for each state at a single position, given a sequence.

# Arguments
- `k::Int`: Position in the sequence.
- `mutated_seq::Array{Int8,1}`: Sequence with a mutation applied.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `L::Int`: Length of the sequence.
- `q::Int`: Number of states for each position (default is 21).

# Returns
- `prob::Array{T}`: Array of conditional probabilities for each state at position `k`.
"""
function single_site_prob_cond(k::Int, mutated_seq::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}, L::Int; q::Int = 21) where {T}
    prob = T.(zeros(q))
    for i in 1:q
        q_k = i
        log_proba = h[q_k, k]
        for j in 1:L
            log_proba += J[mutated_seq[j], j, q_k, k]
        end
        prob[i] = exp(log_proba)
    end
    return normalize(prob, 1)
end

"""
    cont_dep_entr(background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21) where {T}

Compute context-dependent entropy for each position in a sequence.

# Arguments
- `background::Array{Int8,1}`: Background sequence.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `q::Int`: Number of states for each position (default is 21).

# Returns
- `entropy::Array{T}`: Context-dependent entropy for each position.
"""
function cont_dep_entr(background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21) where {T}
    L = size(background, 1)
    prob = hcat([ProbabilityWeights(single_site_prob_cond(site, background, h, J, L, q = q)) for site in 1:L]...)
    return get_entropy(prob, q = q)[:]
end

"""
    single_site_prob_cond_with_T(k::Int, mutated_seq::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}, L::Int; q::Int = 21, temp = 1.0) where {T}

Calculate the conditional probability for each state at a position, with a temperature factor.

# Arguments
- `k::Int`: Position in the sequence.
- `mutated_seq::Array{Int8,1}`: Sequence with a mutation applied.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `L::Int`: Length of the sequence.
- `q::Int`: Number of states for each position (default is 21).
- `temp::Float64`: Temperature parameter (default is 1.0).

# Returns
- `prob::Array{T}`: Array of conditional probabilities for each state at position `k`.
"""
function single_site_prob_cond_with_T(k::Int, mutated_seq::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}, L::Int; q::Int = 21, temp = 1.0) where {T}
    prob = T.(zeros(q))
    for i in 1:q
        q_k = i
        log_proba = h[q_k, k]
        for j in 1:L
            log_proba += J[mutated_seq[j], j, q_k, k]
        end
        prob[i] = exp(log_proba / temp)
    end
    return normalize(prob, 1)
end

"""
    cont_dep_entr_with_T(background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21, temp = 1.0) where {T}

Compute context-dependent entropy for each position with a temperature factor.

# Arguments
- `background::Array{Int8,1}`: Background sequence.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `q::Int`: Number of states for each position (default is 21).
- `temp::Float64`: Temperature parameter (default is 1.0).

# Returns
- `entropy::Array{T}`: Context-dependent entropy for each position.
"""
function cont_dep_entr_with_T(background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21, temp = 1.0) where {T}
    L = size(background, 1)
    prob = hcat([ProbabilityWeights(single_site_prob_cond_with_T(site, background, h, J, L, q = q, temp = temp)) for site in 1:L]...)
    return get_entropy(prob, q = q)[:]
end

"""
    cde_1site(site::Int, background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21) where {T}

Compute context-dependent entropy at a specific position.

# Arguments
- `site::Int`: Position in the sequence.
- `background::Array{Int8,1}`: Background sequence.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `q::Int`: Number of states for each position (default is 21).

# Returns
- `entropy::T`: Context-dependent entropy at the specified position.
"""
function cde_1site(site::Int, background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q = 21) where {T}
    L = length(background)
    prob = ProbabilityWeights(single_site_prob_cond(site, background, h, J, L, q = q))
    return get_entropy(prob, q = q)
end

"""
    CIE(msa::Array{Int8,2}; q = 21)

Compute the context-independent entropy for a multiple sequence alignment.

# Arguments
- `msa::Array{Int8,2}`: Multiple sequence alignment, where each column is a sequence.
- `q::Int`: Number of states for each position (default is 21).

# Returns
- `entropy::Array{Float64}`: Array of context-independent entropies for each position.
"""
function CIE(msa::Array{Int8,2}; q = 21)
    L = size(msa, 1)
    f = reshape(DCAUtils.compute_weighted_frequencies(msa, q+1, 0.2)[1], (q, L))
    return get_entropy(f, q = q)
end
