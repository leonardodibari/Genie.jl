# energy.jl
# Functions for calculating energy-related properties, including average energy and energy variance.

"""
    energy(seq::Array{Int8}, h::Array{T,2}, J::Array{T,4}, L::Int) where {T}

Calculate the total energy for a given sequence `seq`, using the local interaction matrix `h` and the pairwise interaction matrix `J`.

# Arguments
- `seq::Array{Int8}`: Sequence of integers representing states or indices.
- `h::Array{T,2}`: Local field interactions for each position in the sequence.
- `J::Array{T,4}`: Pairwise interactions between sequence positions.
- `L::Int`: Length of the sequence.

# Returns
- `sum::T`: The total computed energy for the sequence.
"""
function energy(seq::Array{Int8}, h::Array{T,2}, J::Array{T,4}, L::Int) where {T}
    sum = zero(T)
    @inbounds for i in 1:L
        sum -= h[seq[i], i]
        @inbounds for j in i+1:L
            sum -= J[seq[i], i, seq[j], j]
        end
    end
    return sum
end

"""
    energy(seq::Array{Int8}, h::Array{T,2}, J::Array{T,4}) where {T}

Compute the total energy for a sequence `seq`, automatically determining its length.

# Arguments
- `seq::Array{Int8}`: Sequence of integers representing states or indices.
- `h::Array{T,2}`: Local field interactions for each position in the sequence.
- `J::Array{T,4}`: Pairwise interactions between sequence positions.

# Returns
- `sum::T`: The total computed energy for the sequence.
"""
function energy(seq::Array{Int8}, h::Array{T,2}, J::Array{T,4}) where {T}
    L = size(seq, 1)
    return energy(seq, h, J, L)
end

"""
    energy(msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4}) where {T}

Calculate the energy for each sequence in a multiple sequence alignment (`msa`) matrix.

# Arguments
- `msa::Array{Int8,2}`: Multiple sequence alignment matrix where each column represents a sequence.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions between sequence positions.

# Returns
- `energies::Array{T}`: Array of energies for each sequence in the alignment.
"""
function energy(msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4}) where {T}
    L, M = size(msa)
    return [energy(msa[:, m], h, J, L) for m in 1:M]
end

"""
    single_mut_dE(seq::Array{Int8, 1}, h::Array{T,2}, J::Array{T,4}, new_aa, mut_pos::Int, L::Int) where {T}

Calculate the change in energy (`delta_E`) due to a single mutation at a specified position.

# Arguments
- `seq::Array{Int8, 1}`: Original sequence.
- `h::Array{T,2}`: Local field interactions.
- `J::Array{T,4}`: Pairwise interactions.
- `new_aa`: New amino acid (or state) introduced by the mutation.
- `mut_pos::Int`: Position of the mutation in the sequence.
- `L::Int`: Length of the sequence.

# Returns
- `delta_E::T`: Change in energy due to the mutation.
"""
function single_mut_dE(seq::Array{Int8, 1}, h::Array{T,2}, J::Array{T,4}, new_aa, mut_pos::Int, L::Int) where {T}
    delta_E = h[seq[mut_pos], mut_pos] - h[new_aa, mut_pos]
    @inbounds for j in 1:L
        delta_E += J[seq[mut_pos], mut_pos, seq[j], j] - J[new_aa, mut_pos, seq[j], j]
    end
    return delta_E
end
