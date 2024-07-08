function proposed_codon(chain,
    all_codons::Array{String},
    seq_site::Int)
    
    beta = 1/64
    #println(chain.DNA[seq_site])
    
    if chain.DNA[seq_site] == "---"
        #if i start from a gap, i go to one of the 64 codons randomly
        #println("POSSIBLE DELETION")
        return rand(all_codons)
    else
        if rand(chain.generator) > 1-beta
        #if i start from a codon, i can go to a gap with probability beta
            return "---"
        end
    end
end

function proposed_codon(codon,
    all_codons::Array{String})
    
    beta = 1/64
    #println(chain.DNA[seq_site])
    
    if codon == "---"
        #if i start from a gap, i propose one of the 64 codons randomly
        #println("POSSIBLE DELETION")
        return rand(all_codons)
    else
        if rand() > 1-beta
        #if i start from a codon, i propose a gap with probability beta
            return "---"
        end
    end
end

function metropolis_indels!(chain,
        new_codon::Union{String, Nothing}, 
        seq_site::Int,
        h::Array{T,2}, 
        J::Array{T,4},
        codon_usage::Dict{String, T},
        temp::T,
        L::Int) where {T}
     
    
    if new_codon == "TAA"  ||  new_codon == "TAG"   || new_codon == "TGA"  || new_codon == chain.seq[seq_site] || new_codon == nothing
        ## in this cases you just refuse the move
    else
        
        new_amino = cod2amino[new_codon]
        dE = single_mut_dE(chain.seq, h, J, new_amino, seq_site, L)
        
        #=new_seq = deepcopy(chain.seq)
        new_seq[seq_site] = new_amino
        J_old = permutedims(J, [1,3,2,4]);
        deltaE = compute_energy_single_sequence(h, J_old, new_seq) - compute_energy_single_sequence(h, J_old, chain.seq)
        println(energy(new_seq,h,J,L) - energy(chain.seq, h, J, L) )
        println("Real deltaE : $(deltaE)")
        println("New method: dE $(dE)")
        
        deltaE = Delta_energy(h, J_old, new_seq, chain.seq)
        println("Old method 2: dE $(deltaE)")
        accept_proba = (length(amino2cod[chain.seq[seq_site]])/length(amino2cod[new_amino]))*exp(-(1/temp)*deltaE)
        println("old method: acceptance probability $(accept_proba)")
        =#
        
        accept_proba = (codon_usage[new_codon]/codon_usage[chain.DNA[seq_site]])*exp(-dE/temp)
        #println("New method: acceptance probability $(accept_proba)")
        if rand(chain.generator) < accept_proba
            chain.seq[seq_site] = new_amino
            chain.DNA[seq_site] = new_codon  
        end
    end  
    
end



function run_metropolis_indels!(chains, 
        h::Array{T,2}, 
        J::Array{T,4},
        all_codons::Array{String},
        codon_usage::Dict{String, T},
        N_chains::Int,
        temp::T,
        L::Int) where {T}
    
    @tasks for n in 1:N_chains
        seq_site = rand(chains[n].generator, 1:L)
        new_codon = proposed_codon(chains[n], all_codons, seq_site)
        #mutated_seq = SeqToEvolve(chains[n].seq, chains[n].DNA)
        #metro_del_ins_step(chains[n], seq_site, new_codon, h, J, L, temp)
        metropolis_indels!(chains[n], new_codon, seq_site, h, J, codon_usage, temp, L)
    end 
    
end




###old functions

function get_accessible_nucleo_for_del_ins(old_cod)

    if old_cod .== "---"
        amino_list = [i for i in 0:20]
        amino_list[1] = 21
        codon_list = reduce(vcat,[amino2cod[a] for a in amino_list])
        #push!(codon_list, old_cod )
        push!(codon_list, "TAG")
        push!(codon_list, "TAA")
        push!(codon_list, "TGA")
        return amino_list, codon_list
    else
        codon_list = ["---"]
        push!(codon_list, old_cod)
        amino_list = get.(Ref(cod2amino), codon_list, 0)
        #println(amino_list)
        return amino_list, codon_list
    end
end
   

function del_ins_codon_sampling(arr)
    n = length(arr)
    beta = 1/64
    alpha = 1-64*beta
    gamma = 1-beta
    if n==2
        probabilities = [beta, gamma]
        index = sample(1:n, Weights(probabilities))
        return arr[index]
    else 
        probabilities = fill(beta, n)
        probabilities[1] = alpha
        index = sample(1:n, Weights(probabilities))
        return arr[index]
    end
end


function metro_del_ins_step(chain, pos_mut, new_codon, h, J, N, T)
   
        
    #pos_mut = rand(1:length(mutated_seq.Amino))

	old_codon = chain.DNA[pos_mut]
    old_amino = chain.seq[pos_mut]
    #println("Old method: mutating codon $(old_codon) at site $(pos_mut) ")
    amino_list, codon_list = get_accessible_nucleo_for_del_ins(old_codon) 
    new_codon = del_ins_codon_sampling(codon_list)
    
    #println("Old method: Proposed codon  $(new_codon)")
    
    if new_codon == "TAA"  ||  new_codon == "TAG"   || new_codon == "TGA"  || new_codon == old_codon || new_codon == nothing
        
    else
        new_amino = cod2amino[new_codon]  
        new_sequence = deepcopy(chain.seq)
        new_sequence[pos_mut] = new_amino
        deltaE = Delta_energy(h, J, new_sequence, chain.seq)
        #println("Old method: dE $(deltaE)")
        accept_proba = (length(amino2cod[old_amino])/length(amino2cod[new_amino]))*exp(-(1/T)*deltaE)
        #println("Old method: acceptance probability $(accept_proba)")
        
        if rand() < accept_proba
            aa = cod2amino[new_codon]
            chain.DNA[pos_mut] = new_codon	
            chain.seq[pos_mut] = aa
        end
    end
end