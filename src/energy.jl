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
   

#### old functions



function Delta_energy(h::Array{Float64,2}, J::Array{Float64,4}, 
        S::Array{<:Integer,1}, ref::Array{<:Integer, 1})
    q, N = size(h)
    E = 0.0
    
    index_v = collect(1:N)
    common = (S .== ref)
    idx = findfirst(isequal(false), common)
    common = index_v[common]
    E -= (h[S[idx],idx] - h[ref[idx],idx])
    @fastmath for j = 1:N
        if j > idx
            @inbounds  E -= (J[S[j],S[idx],j,idx] - J[ref[j],ref[idx],j,idx] )
        end
    end
    
    @fastmath for i in common
        if idx > i
            @inbounds  E -= (J[S[idx],S[i],idx,i] - J[ref[idx],ref[i],idx,i] )
        end
    end
    
    return E
end






###old sequences

function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)
    N = size(h)[2]
    q = size(h)[1]
    E = 0.0
    for i = 1:N
        E -= h[S[i],i]
        for j = (i+1):N
			E -= J[S[i],S[j],i,j]
		end
	end
return E
end

