function log_prob!(log_prob::Array{T,1}, 
        chain, 
        h::Array{T,2}, 
        J::Array{T,4},
        codon_usage::Dict{String, T},
        codon_list::Array{String,1},
        temp::T,
        seq_site::Int,
        L::Int) where {T}
    
    @inbounds for x in 1:length(log_prob)
        a = cod2amino[codon_list[x]] 
        log_prob[x] = h[a,seq_site] +temp*log(codon_usage[codon_list[x]])
        @inbounds for j in 1:L
            log_prob[x] += J[chain.seq[j], j, a, seq_site]
        end
    end
    
end
      

function get_prob!(chain, 
        amino_list::Array{Int8,1}, 
        codon_list::Array{String,1},
        log_prob::Array{T,1},
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        nucleo_site::Int,
        codon_usage::Dict{String, T},
        temp::T,
        L::Int) where {T}
       
    if all_equal(amino_list) == true
        chain.seq[seq_site] == amino_list[1]
        chain.DNA[nucleo_site] == rand(codon_list)
    else
        log_prob!(log_prob, chain, h, J, codon_usage, codon_list, temp, seq_site, L)
        log_prob ./= temp
        loc_softmax!(log_prob)
        loc_sample!(chain.generator, log_prob, codon_list, chain.seq, chain.DNA, seq_site)
    end
end

function prob_cond!(chain, 
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        nucleo_site::Int,
        codon_net::Dict{String, Dict{Int64, Vector{String}}}, 
        codon_usage::Dict{String, T},
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        temp::T,
        L::Int) where {T}
    
    L_moves = length_of_moves[chain.DNA[seq_site], nucleo_site]
    
    if L_moves == 2
        
        chain.codon_list2 .= accessible_codons(chain.DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list2[idx] = cod2amino[chain.codon_list2[idx]]
        end
        get_prob!(chain, chain.amino_list2, chain.codon_list2, chain.log_prob2, h, J, seq_site, 
            nucleo_site, codon_usage, temp, L)
       
    elseif L_moves == 3
        
        chain.codon_list3 .= accessible_codons(chain.DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list3[idx] = cod2amino[chain.codon_list3[idx]]
        end
        get_prob!(chain, chain.amino_list3, chain.codon_list3, chain.log_prob3, h, J, seq_site, 
            nucleo_site, codon_usage, temp, L)
        
    elseif L_moves == 4
        
        chain.codon_list4 .= accessible_codons(chain.DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list4[idx] = cod2amino[chain.codon_list4[idx]]
        end
        get_prob!(chain, chain.amino_list4, chain.codon_list4, chain.log_prob4, h, J, seq_site, 
            nucleo_site, codon_usage, temp, L)
    end

end

function return_pos_rng(rng, seq::Array{Int8,1},L::Int)
    idx = rand(rng, 1:L)
    if seq[idx] == 21
        return_pos_rng(rng, seq, L)
    else
        return idx
    end
end


function run_gibbs_sampling!(chains, 
        h::Array{T,2}, 
        J::Array{T,4}, 
        codon_net::Dict{String, Dict{Int64, Vector{String}}}, 
        codon_usage::Dict{String, T},
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        N_chains::Int,
        temp::T,
        L::Int) where {T}
    
    @tasks for n in 1:N_chains
        seq_site = return_pos_rng(chains[n].generator, chains[n].seq, L)   
        #seq_site = rand(chains[n].generator, findall(x -> x != 21, chains[n].seq))
        nucleo_site = rand(chains[n].generator,1:3)
        prob_cond!(chains[n], h, J, seq_site, nucleo_site, codon_net, codon_usage, length_of_moves, temp, L)
    end
end





