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

function energy(seq::Array{Int8}, h::Array{T,2}, J::Array{T,4}) where {T}
    L = size(seq,1)
    return energy(seq, h, J, L)
end

function energy(msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4}) where {T}
    L,M = size(msa)
    return [energy(msa[:,m], h, J, L) for m in 1:M]
end

function single_mut_dE(seq::Array{Int8, 1}, h::Array{T,2}, J::Array{T,4}, new_aa, mut_pos::Int, L::Int) where {T}
    delta_E = h[seq[mut_pos], mut_pos] - h[new_aa, mut_pos]
    @inbounds for j in 1:L
        delta_E += J[seq[mut_pos], mut_pos, seq[j], j] - J[new_aa, mut_pos, seq[j], j]
    end
    return delta_E
end
    
    
function single_mut_dE(seq::Array{Int, 1}, h::Array{T,2}, J::Array{T,4}, new_aa, mut_pos::Int, L::Int) where {T}
    delta_E = h[seq[mut_pos], mut_pos] - h[new_aa, mut_pos]
    @inbounds for j in 1:L
        delta_E += J[seq[mut_pos], mut_pos, seq[j], j] - J[new_aa, mut_pos, seq[j], j]
    end
    return delta_E
end
   

