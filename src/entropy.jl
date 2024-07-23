function get_entropy(f::Array{T,2}; q::Int = 21) where {T}
    N=length(f[1,:])
    entr = zeros(Float64, N)
    for i in 1:N
        for a in 1:q
            if(f[a,i]>0)
                entr[i]-=f[a,i]*log(f[a,i])
            end
        end
    end
    
    return entr / log(2)
end


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
	return normalize(prob,1)
    
end


function cont_dep_entr(background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q =21) where {T}
    
    L = size(background,1)
    
    prob = hcat([ProbabilityWeights(single_site_prob_cond(site, background, h, J, L, q=q)) for site in 1:L]...)

    return get_entropy(prob, q = q)[:]
end

function cde_1site(site::Int, background::Array{Int8,1}, h::Array{T,2}, J::Array{T,4}; q =21) where {T}
  
    L = length(background)
    prob = ProbabilityWeights(single_site_prob_cond(site,  background, h, J, L, q = q))
    return get_entropy(prob, q = q)
end


function CIE(msa::Array{Int8,2}; q = 21)
    L = size(msa,1)
    f = reshape(DCAUtils.compute_weighted_frequencies(msa, q+1, 0.2)[1], (q, L))
    return get_entropy(f, q = q)
end
