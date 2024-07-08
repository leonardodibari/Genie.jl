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
       
    #println("New codon and amino list")
    #println(codon_list, amino_list)
    
    if all_equal(amino_list) == true
        chain.seq[seq_site] == amino_list[1]
        chain.DNA[nucleo_site] == rand(codon_list)
    else
        log_prob!(log_prob, chain, h, J, codon_usage, codon_list, temp, seq_site, L)
        log_prob ./= temp
        loc_softmax!(log_prob)
        #println("New probability")
        #println(log_prob)
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
    
    
       
    #amino_list2, codon_list2 =  get_accessible_nucleo_muts_DNA_det_bal(chain.DNA[seq_site], nucleo_site)
    
    #println("Old codon and amino list")
    #println(codon_list2, amino_list2)
    
    #print("Old probability")
    #println(nucleo_cond_proba_DNA_gibbs_new(seq_site, amino_list2, codon_list2, chain.seq, h, permutedims(J, [1,3,2,4]), L,  temp = temp))
    
    #=if length(unique(amino_list)) == 1
        chain.seq[seq_site] == cod2amino[codon_list[1]]
        chain.DNA[nucleo_site] == rand(codon_list)
    else
        log_prob = T.(zeros(length(codon_list)))
        log_prob!(log_prob, chain, h, J, codon_usage, codon_list, temp, seq_site, L)
        log_prob ./= temp
        loc_softmax!(log_prob)
        #println("New probability")
        #println(log_prob)
        loc_sample!(chain.generator, log_prob, codon_list, chain.seq, chain.DNA, seq_site)
end=#
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
        seq_site = rand(chains[n].generator, findall(x -> x != 21, chains[n].seq))
        nucleo_site = rand(chains[n].generator,1:3)
        prob_cond!(chains[n], h, J, seq_site, nucleo_site, codon_net, codon_usage, length_of_moves, temp, L)
    end
end





### old functions

function get_accessible_nucleo_muts_DNA_det_bal(old_codon, nucleo_pos::Integer)
    old_codon = [string(old_codon[i]) for i in 1:3 ]
	codon_list = Vector{AbstractString}(undef, 4)
	new_codon = deepcopy(old_codon)
	for (j, nucl) in enumerate(["A", "C", "G", "T"]) 
		new_codon[nucleo_pos] = nucl
		codon_list[j] = join(new_codon)
	end
    
    codon_list = filter!(aa->aa != "TAA", codon_list)
    codon_list = filter!(aa->aa != "TAG", codon_list)
    codon_list = filter!(aa->aa != "TGA", codon_list)
    
	amino_list = get.(Ref(cod2amino), codon_list, 0)
	#amino_list = filter!(aa->aa != 21, amino_list)

	return amino_list, codon_list
end


function sum_couplings_inbounds(mutated_seq, q_k, k, J::Array{Float64,4}, N)
    result = 0
    @inbounds for i in 1:N
        result += J[mutated_seq[i], q_k ,i, k]
    end
    return result
end 


function log_p(mutated_seq, q_k, k, h, J::Array{Float64,4}, N::Integer; temp=1)
    
    log_proba = h[q_k, k] -temp*log(length(amino2cod[q_k]))
    log_proba += sum_couplings_inbounds(mutated_seq, q_k, k, J, N)
    return log_proba  
end

function nucleo_cond_proba_DNA_gibbs_new(k, q, codon_list, mutated_seq, h::Array{Float64,2}, J::Array{Float64,4}, N::Integer;  temp = 1)
	prob = zeros(length(codon_list))
        
    d = Dict{Int64, Float64}()
    
	@inbounds for (index, q_k) in enumerate(q)
            
        if haskey(d, q_k)
            log_proba = d[q_k]
            #println("not computed")
        else
            log_proba = log_p(mutated_seq, q_k, k, h, J, N, temp = temp)
            d[q_k] = log_proba
            #println("computed")
        end
		prob[index] = exp(log_proba/temp)
	end
        
	return normalize(prob,1)
end

