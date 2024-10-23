Base.@kwdef mutable struct Seq <: TreeNodeData # Create a custom data type
    seq::Array{Int,1} = [0,0]
    DNA::Array{String,1} = ["A","b"]
end

Base.copy(s::Seq) = Seq(copy(s.seq), copy(s.DNA))

struct Storage{Ti, T}
    tree::Tree{Seq}
    generator::Xoshiro
    codon_list2::Vector{String}
    codon_list3::Vector{String}
    codon_list4::Vector{String}
    amino_list2::Vector{Ti}
    amino_list3::Vector{Ti}
    amino_list4::Vector{Ti}
    log_prob2::Vector{T}
    log_prob3::Vector{T}
    log_prob4::Vector{T}
end


function Storage(start_seq::Array{Int,1}, tree_file::String, generator::Xoshiro; T::DataType = Float64, Ti::DataType=Int)
    
    tree = read_tree(tree_file, node_data_type = Seq)
    dna_seq = String.(amino_seq2dna_seq(start_seq))
    for l in keys(tree.lnodes)
        data!(tree[l], Seq(seq = copy(start_seq), DNA = copy(dna_seq)))
    end
    codon_list2=Vector{String}(undef, 2);codon_list3=Vector{String}(undef, 3);codon_list4 =Vector{String}(undef, 4);
    amino_list2=Ti.([1,1]); amino_list3=Ti.([1,1,1]); amino_list4=Ti.([1,1,1,1]);
    log_prob2 = T.([0.,0.]); log_prob3 = T.([0.,0.,0.]); log_prob4 = T.([0.,0.,0.,0.]);
    Storage{Ti, T}(tree, generator, codon_list2, codon_list3, codon_list4, 
        amino_list2, amino_list3, amino_list4, 
        log_prob2, log_prob3, log_prob4)
end



function log_prob_tree!(log_prob::Array{T,1}, 
        seq::Array{Int,1}, 
        DNA::Array{String,1}, 
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
            log_prob[x] += J[seq[j], j, a, seq_site]
        end
    end
    
end
      


function get_prob_tree!(seq::Array{Int,1}, 
        DNA::Array{String,1}, 
        amino_list::Array{Int,1}, 
        codon_list::Array{String,1},
        log_prob::Array{T,1},
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        nucleo_site::Int,
        codon_usage::Dict{String, T},
        rng::Xoshiro,
        temp::T,
        L::Int) where {T} 
    
    if all_equal(amino_list) == true
        seq[seq_site] == amino_list[1]
        DNA[nucleo_site] == rand(codon_list)
    else
        log_prob_tree!(log_prob, seq, DNA, h, J, codon_usage, codon_list, temp, seq_site, L)
        log_prob ./= temp
        loc_softmax!(log_prob)
        loc_sample!(rng, log_prob, codon_list, seq, DNA, seq_site)
    end
end


function prob_cond_tree!(chain::Storage{Int64, Float64}, 
        seq::Array{Int,1}, 
        DNA::Array{String,1},
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        nucleo_site::Int,
        codon_net::Dict{String, Dict{Int64, Vector{String}}}, 
        codon_usage::Dict{String, T},
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        temp::T,
        L::Int) where {T}
    
    L_moves = length_of_moves[DNA[seq_site], nucleo_site]
    
    if L_moves == 2
        
        chain.codon_list2 .= accessible_codons(DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list2[idx] = cod2amino[chain.codon_list2[idx]]
        end
        get_prob_tree!(seq, DNA, chain.amino_list2, chain.codon_list2, chain.log_prob2, h, J, seq_site, 
            nucleo_site, codon_usage, chain.generator, temp, L)
       
    elseif L_moves == 3
        
        chain.codon_list3 .= accessible_codons(DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list3[idx] = cod2amino[chain.codon_list3[idx]]
        end
        get_prob_tree!(seq, DNA, chain.amino_list3, chain.codon_list3, chain.log_prob3, h, J, seq_site, 
            nucleo_site, codon_usage, chain.generator, temp, L)
        
    elseif L_moves == 4
        
        chain.codon_list4 .= accessible_codons(DNA[seq_site], codon_net, nucleo_site)
        for idx in 1:L_moves
            chain.amino_list4[idx] = cod2amino[chain.codon_list4[idx]]
        end
        get_prob_tree!(seq, DNA, chain.amino_list4, chain.codon_list4, chain.log_prob4, h, J, seq_site, 
            nucleo_site, codon_usage, chain.generator, temp, L)
    end
    
end


function return_pos(seq::Array{Int,1},L::Int)
    idx = rand(1:L)
    if seq[idx] == 21
        return_pos(seq,L)
    else
        return idx
    end
end

function run_gibbs_sampling_tree!(chain::Storage{Int64, Float64}, 
        seq::Array{Int,1}, 
        DNA::Array{String,1}, 
        h::Array{T,2}, 
        J::Array{T,4}, 
        codon_net::Dict{String, Dict{Int64, Vector{String}}}, 
        codon_usage::Dict{String, T},
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        temp::T,
        L::Int) where {T}

    seq_site = return_pos(seq, L)    
    #println("Site $(seq_site) codon $(DNA[seq_site]) amino $(seq[seq_site])")
    nucleo_site = rand(1:3);
    prob_cond_tree!(chain, seq, DNA, h, J, seq_site, nucleo_site, codon_net, codon_usage, length_of_moves, temp, L)
end



function proposed_codon_tree(DNA::Array{String,1}, 
    chain::Storage{Int64, Float64},
    all_codons::Array{String},
    seq_site::Int)
    
    beta = 1/64
    #println(chain.DNA[seq_site])
    
    if DNA[seq_site] == "---"
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


function metropolis_indels_tree!(chain::Storage{Int64, Float64}, 
        seq::Array{Int,1},
        DNA::Array{String,1},
        new_codon::Union{String, Nothing}, 
        seq_site::Int,
        h::Array{T,2}, 
        J::Array{T,4},
        codon_usage::Dict{String, T},
        temp::T,
        L::Int) where {T}
     
    
    if new_codon == "TAA"  ||  new_codon == "TAG"   || new_codon == "TGA"  || new_codon == seq[seq_site] || new_codon == nothing
        ## in this cases you just refuse the move
    else
        
        new_amino = cod2amino[new_codon]
        dE = single_mut_dE(seq, h, J, new_amino, seq_site, L)
        accept_proba = (codon_usage[new_codon]/codon_usage[DNA[seq_site]])*exp(-dE/temp)
        if rand(chain.generator) < accept_proba
            seq[seq_site] = new_amino
            DNA[seq_site] = new_codon  
        end
    end  
    
end



function run_metropolis_indels_tree!(chain::Storage{Int64, Float64}, 
        seq::Array{Int,1},
        DNA::Array{String,1},
        h::Array{T,2}, 
        J::Array{T,4},
        all_codons::Array{String},
        codon_usage::Dict{String, T},
        temp::T,
        L::Int) where {T}
    
    seq_site = rand(chain.generator, 1:L)
    new_codon = proposed_codon_tree(DNA, chain, all_codons, seq_site)
    metropolis_indels_tree!(chain, seq, DNA, new_codon, seq_site, h, J, codon_usage, temp, L)
end


function assign_sequences!(node::TreeNode{Seq}, 
        chain::Storage{Int64, Float64},
        h::Array{T,2}, 
        J::Array{T,4}, 
        codon_net::Dict{String, Dict{Int64, Vector{String}}}, 
        all_codons::Array{String},
        codon_usage::Dict{String, T},
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        temp::T,
        mu::Float64,
        p::Float64,
        L::Int) where {T}
    
    if isempty(node.child)
        return 0
    end
    

    for a in node.child
        
        for i in 1:L
            data(a).seq[i] = data(a.anc).seq[i] 
            data(a).DNA[i] = data(a.anc).DNA[i]
        end
                
        estimate_steps = L*mu*branch_length(a)
        steps = floor(Int, estimate_steps)
        if rand() < (estimate_steps - steps)
            steps += 1
        end
        
        for _ in 1:steps
            if rand() < 1 - p
                #sampling gibbs with probability p
                run_gibbs_sampling_tree!(chain, data(a).seq, data(a).DNA, h, J, 
                    codon_net, codon_usage, length_of_moves, temp, L)
            else
                #sampling metropolis with probability 1-p
                run_metropolis_indels_tree!(chain, data(a).seq, data(a).DNA, h, J, 
                    all_codons, codon_usage, temp, L)
            end 
            #println("Node $(a.label) $(sum(data(a).seq .== data(a.anc).seq))")
        end
        
        assign_sequences!(a, chain, h, J, 
            codon_net, all_codons, codon_usage, length_of_moves, temp, mu, p, L)
    end
end




function run_evolution_ontree(start_seq::Union{Array{Int,1}, Array{Int,2}}, tree_file::String, h::Array{T,2},   J::Array{T,4};
        temp::Float64 = 1.0, 
        mu::Float64 = 1.0,
        p::Float64 = 0.5, 
        q = 21, 
        codon_bias::Union{Nothing, Dict{String, Float64}} = nothing, 
        verbose = false) where {T}
    
    
    L = size(start_seq, 1)
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    end
    
    
    codon_net = create_nested_codon_dict()
    length_of_moves = create_length_dict(codon_net)
    temp = T(temp)
    all_codons = vcat([amino2cod[i] for i in 1:20]...)
    push!(all_codons, "TAG")
    push!(all_codons, "TAA")
    push!(all_codons, "TGA")
    
    if codon_bias == nothing
        no_cod_bias = Dict(x => T(1/length(amino2cod[cod2amino[x]])) for x in keys(cod2amino))
        codon_usage = no_cod_bias
    end
    
    N_trees = size(start_seq,2)    
    
    if N_trees == 1
        rng = random_gens(N_trees+1)  
        @time chains = [Storage(start_seq, tree_file, rng[n]) for n in 1:N_trees+1]
        @time assign_sequences!(chains[1].tree.root, chains[1], h, J, 
            codon_net, all_codons, codon_usage, length_of_moves, temp, mu, p, L)
        return chains[1].tree
        
    elseif N_trees > 1
        
        trees = []
        rng = random_gens(N_trees)  
        chains = [Storage(start_seq[:,n], tree_file, rng[n]) for n in 1:N_trees]
        @tasks for n in 1:N_trees
            assign_sequences!(chains[n].tree.root, chains[n], h, J, 
                codon_net, all_codons, codon_usage, length_of_moves, temp, mu, p, L)
        end 

        return [chains[i].tree for i in 1:N_trees]
        
    elseif (N_trees == 0) || (N_trees < 0)
        error("N_trees must be >= 1")
    end
end
    



function msa_from_leafs(tree)
    msa = []; for a in keys(tree.lleaves)
    push!(msa, data(tree[a]).seq) end; msa = hcat(msa...);
    return msa
end

function msa_from_nodes(tree)
    msa = []; for a in keys(tree.lnodes)
    push!(msa, data(tree[a]).seq) end; msa = hcat(msa...);
    return msa
end


function find_optimal_mu(tree_file::String, mus::Array{Float64,1}, start_seq::Array{Int,1}, nat_msa::Array{Int8,2}, h::Array{Float64,2}, J::Array{Float64,4}; n_seq = 1000)
    
    f_nat = pair_dist_freq(nat_msa, n_seq = n_seq); n_max = argmax(f_nat); res = [];
    for mu in mus
        tree = run_evolution_ontree(start_seq, tree_file, h, J, mu = mu, p = 0.5); 
        msa = msa_from_leafs(tree); 
        f_sim = pair_dist_freq(msa, n_seq = n_seq); 
        push!(res,kldivergence(f_sim[1:n_max]./sum(f_sim[1:n_max]), f_nat[1:n_max]./sum(f_nat[1:n_max])))
    end
    println("Optimal mu: $(mus[argmin(res)])")
    return (mus = mus, score = res)
end
