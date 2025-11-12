struct SeqWithParent
    seq::Vector{String}
    parent::Union{Vector{String}, Nothing}
end

Base.:(==)(a::SeqWithParent, b::SeqWithParent) = a.seq == b.seq && a.parent == b.parent
Base.hash(x::SeqWithParent, h::UInt) = hash((x.seq, x.parent), h)



function random_mutations_nucleo_mcmc_inspired!(
    msa::Matrix{Int8}, msa_dna::Matrix{String}, 
    codon_net::Dict{String, Dict{Int, Vector{String}}}, 
    length_of_moves::Dict{Tuple{String, Int64}, Int64},
    codon_list::Vector{Vector{String}},
    abund_dict::Dict{SeqWithParent, Int64},   # << change here!
    mu::T; 
    q::Int = 4,
    verbose = false
) where {T}

    L, M = size(msa)
    distrib = Binomial(L, prob_mut_one_site(mu, q = q))

    all_nucleos = ['A','C', 'G', 'T']
    empty!(abund_dict)  # clear dict before filling

    for m in 1:M
        parent_seq = [msa_dna[seq_site, m] for seq_site in 1:L]  # collect original sequence as Vector{String}
        K = rand(distrib)  # number of mutations 
        pos = sample(1:3*L, K, replace=false)  # mutation positions

        child_seq = copy(parent_seq)  # start from parent sequence

        for k in 1:K
            seq_site = div(pos[k]-1, 3) + 1
            nucleo_site = mod(pos[k]-1, 3) + 1

            if child_seq[seq_site] == "21"
                # no mutation if 21 (stop codon?) - adjust logic as needed
            else
                L_moves = length_of_moves[(child_seq[seq_site], nucleo_site)]
                if L_moves == 1
                    child_seq[seq_site] = accessible_codons(child_seq[seq_site], codon_net, nucleo_site)[1]
                else
                    codon_list[L_moves] .= accessible_codons(child_seq[seq_site], codon_net, nucleo_site)
                    child_seq[seq_site] = sample(codon_list[L_moves])
                end
            end
        end

        # Create struct storing child and parent sequences
        seq_struct = SeqWithParent(child_seq, parent_seq)

        # Update dictionary counts
        abund_dict[seq_struct] = get(abund_dict, seq_struct, 0) + 1
    end

    clean_dict_nucleo!(abund_dict)  # optional, if you want to clean zeros

    if verbose
        println("$(100 * length(abund_dict) / M) % unique seqs after mutagenesis")
    end
end

function random_mutations_nucleo_mcmc!(
    msa::Matrix{Int8}, msa_dna::Matrix{String}, 
    codon_net::Dict{String, Dict{Int, Vector{String}}}, 
    length_of_moves::Dict{Tuple{String, Int64}, Int64},
    codon_list::Vector{Vector{String}},
    abund_dict::Dict{SeqWithParent, Int64}, 
    mu::T; 
    q::Int = 4,
    temp::T = 1.,
    verbose = false
) where {T}

    L, M = size(msa)
    distrib = Binomial(L, prob_mut_one_site(mu, q = q))

    all_nucleos = ['A','C', 'G', 'T']
    empty!(abund_dict)  # clear dict before filling

    for m in 1:M
        parent_seq = [msa_dna[seq_site, m] for seq_site in 1:L]  # collect original sequence as Vector{String}
        K = rand(distrib)  # number of mutations 
        pos = sample(1:3*L, K, replace=false)  # mutation positions

        child_seq = copy(parent_seq)  # start from parent sequence

        for k in 1:K
            seq_site = div(pos[k]-1, 3) + 1
            nucleo_site = mod(pos[k]-1, 3) + 1

            if child_seq[seq_site] == "21"
                # no mutation if 21 (stop codon?) - adjust logic as needed
            else
                L_moves = length_of_moves[(child_seq[seq_site], nucleo_site)]
                if L_moves == 1
                    child_seq[seq_site] = accessible_codons(child_seq[seq_site], codon_net, nucleo_site)[1]
                    en_child = energy_dna2(child_seq, h, J, L)
                    en_parent = energy_dna2(parent_seq, h, J, L)
                    if rand()<exp(-(en_child-en_parent)/temp)
                        # Create struct storing child and parent sequences
                        seq_struct = SeqWithParent(child_seq, parent_seq)
                        # Update dictionary counts
                        abund_dict[seq_struct] = get(abund_dict, seq_struct, 0) + 1
                    else
                        # Create struct storing parent sequence
                        seq_struct = SeqWithParent(parent_seq, parent_seq)
                        # Update dictionary counts
                        abund_dict[seq_struct] = get(abund_dict, seq_struct, 0) + 1
                    end
                        
                else
                    codon_list[L_moves] .= accessible_codons(child_seq[seq_site], codon_net, nucleo_site)
                    child_seq[seq_site] = sample(codon_list[L_moves])
                    en_child = energy_dna2(child_seq, h, J, L)
                    en_parent = energy_dna2(parent_seq, h, J, L)
                    if rand()<exp(-(en_child-en_parent)/temp)
                        # Create struct storing child and parent sequences
                        seq_struct = SeqWithParent(child_seq, parent_seq)
                        # Update dictionary counts
                        abund_dict[seq_struct] = get(abund_dict, seq_struct, 0) + 1
                    else
                        # Create struct storing parent sequence
                        seq_struct = SeqWithParent(parent_seq, parent_seq)
                        # Update dictionary counts
                        abund_dict[seq_struct] = get(abund_dict, seq_struct, 0) + 1
                    end
                end
            end
        end
    end

    clean_dict_nucleo!(abund_dict)  # optional, if you want to clean zeros

    if verbose
        println("$(100 * length(abund_dict) / M) % unique seqs after mutagenesis")
    end
end


function selection_nucleo_mcmc_inspired!(
    h::Array{T,2}, J::Array{T,4}, abund_dict::Dict{SeqWithParent, Int64}, 
    codon_usage::Dict{String, Float64}, L::Int, M::Int; 
    temp::T = 1., verbose = false
) where {T}

    seq_structs = collect(keys(abund_dict))

    aa = zeros(length(seq_structs))
    surv_prob = Vector{Float64}(undef, length(seq_structs))
    ens = zeros(length(seq_structs))
    mu_bind = zeros(length(seq_structs))

    @tasks for i in 1:length(seq_structs)
        child_seq = seq_structs[i].seq
        parent_seq = seq_structs[i].parent

        ens[i] = energy_dna2(child_seq, h, J, L)

        # Compute parent's energy (if parent exists)
        if parent_seq === nothing
            mu_bind[i] = 0.0  # or some default value if no parent
        else
            mu_bind[i] = energy_dna2(parent_seq, h, J, L)
        end
    end

    @tasks for i in 1:length(seq_structs)
        seq = seq_structs[i].seq
        if any(codon -> codon == "TAG" || codon == "TAA" || codon == "TGA", seq)
            surv_prob[i] = 0.0
        else
            aa[i] = exp(-((ens[i] - mu_bind[i]) / temp))
            p = 1.0  # or your inv_degeneracy(seq, codon_usage)
            surv_prob[i] = (p * aa[i]) / (1 + aa[i])
        end

        distrib = Binomial(abund_dict[seq_structs[i]], surv_prob[i])
        abund_dict[seq_structs[i]] = rand(distrib)
    end

    clean_dict_nucleo!(abund_dict)

    if verbose
        println("$(100 * length(abund_dict) / M) % unique seqs after selection")
    end
end


function amplification_nucleo_mcmc_inspired!(
    final_msa::Array{Int8,2}, final_msa_dna::Array{String,2}, 
    abund_dict::Dict{SeqWithParent, Int64}
)
    L, Mf = size(final_msa)
    
    seq_structs = collect(keys(abund_dict))
    
    # Extract child sequences from SeqWithParent keys
    seqs_dna = [s.seq for s in seq_structs]
    
    w = [abund_dict[s] for s in seq_structs]
    
    M = length(seqs_dna)
    
    idxs = sample(1:M, Weights(w), Mf)
    
    @tasks for m in 1:Mf
        for i in 1:L
            final_msa_dna[i, m] = seqs_dna[idxs[m]][i]
            final_msa[i, m] = cod2amino[final_msa_dna[i, m]]
        end
    end
end


function run_dir_evol_nucleo_mcmc_inspired(
        start_msa::Array{Int8,2}, start_msa_dna::Array{String, 2}, h::Array{T,2}, J::Array{T,4};
        rounds::Int = 4, 
        each_step::Union{Int, Nothing} = nothing,
        seq_steps::Union{Int, Nothing} = nothing,
        seq_reads::Int = 100,
        temp::Float64 = 1.0,  
        mu::T = 10^-2,
        q::Int = 21,
        codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
        mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
        verbose = false
    ) where {T}
    
    L,M = size(start_msa);
    
    if M < seq_reads
        error("Increase number of sequences, seq_reads must be smaller than amount of sequences")
    end
    
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    elseif (size(J,1) !== size(h,1)) || (size(J,3) !== size(h,1))
        error("Issues with color size")
    elseif seq_steps !== nothing && each_step !== nothing
        error("You cannot specify both seq_steps and each_step: 
            \n-if you specify seq_steps, you get a msa on the times you suggested along the trajectory
            \n-if you specify each_step, you get a msa every each_step along the trajectory
            \n-if you don't specify anything then you get only the msa at the final round \n")
    end
      
    if each_step !== nothing
        if each_step > rounds
            error("each_step must be smaller than rounds")
        end
        steps = [i for i in 1:each_step:rounds] 
        step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
        step_msa_dna = [Matrix{String}(undef, L, seq_reads) for i in 1:length(steps)]
        pre_sel_step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
        pre_sel_step_msa_dna = [Matrix{String}(undef, L, seq_reads) for i in 1:length(steps)]
        post_sel_step_msa = Vector{Matrix{Int8}}()
        post_sel_step_msa_dna = Vector{Matrix{String}}()
    end
    
    if seq_steps !== nothing
        if maximum(seq_steps) > rounds
            error("Sequencing steps must be smaller than rounds")
        elseif is_positive_increasing(seq_steps) == false
            error("Sequencing steps must be strictly positive and increasing")
        end
        steps = deepcopy(seq_steps)
        step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
        step_msa_dna = [Matrix{String}(undef, L, seq_reads) for i in 1:length(steps)]
        pre_sel_step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
        pre_sel_step_msa_dna = [Matrix{String}(undef, L, seq_reads) for i in 1:length(steps)]
        post_sel_step_msa = Vector{Matrix{Int8}}()
        post_sel_step_msa_dna = Vector{Matrix{String}}()
    end
    
    codon_net = create_nested_codon_dict_no_same()
    length_of_moves = create_length_dict(codon_net)
    count = 0
    count2 = 0
    count3 = 0
    all_codons = vcat([amino2cod[i] for i in 1:20]...)
    push!(all_codons, "TAG")
    push!(all_codons, "TAA")
    push!(all_codons, "TGA")
    final_msa = deepcopy(start_msa);
    final_msa_dna = deepcopy(start_msa_dna);
    
    # Initialize abund_dict_nucleo with SeqWithParent (parent = nothing for initial sequences)
    abund_dict_nucleo = Dict{SeqWithParent, Int}()
    for m in 1:M
        seq_copy = copy(start_msa_dna[:, m])
        abund_dict_nucleo[SeqWithParent(seq_copy, nothing)] = 1
    end
    
    codon_list = [Vector{String}(undef, i) for i in 1:4]
    prob_codon_list = [Vector{Float64}(undef, i) for i in 1:4]

    if codon_bias == nothing
        no_cod_bias = Dict(x => T(1/length(amino2cod[cod2amino[x]])) for x in keys(cod2amino))
        codon_usage = no_cod_bias
    else
        codon_usage = codon_bias
    end
    
    if mut_bias !== nothing
        trans_dict = create_codon_transition_dict(mut_bias)
    end
        
    for r in 1:rounds
        if verbose == true
            println()
            println("Round $(r)")
        end
        
        Genie.random_mutations_nucleo_mcmc_inspired!(final_msa, final_msa_dna, codon_net, 
                length_of_moves, codon_list, abund_dict_nucleo, mu, verbose = verbose);
    
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            count2+=1
            idxs = rand(1:M,seq_reads)
            pre_sel_step_msa[count2] .= final_msa[:,idxs]
            pre_sel_step_msa_dna[count2] .= final_msa_dna[:,idxs]
        end 
        
        Genie.selection_nucleo_mcmc_inspired!(h, J, abund_dict_nucleo, codon_usage, L, M, temp = temp, verbose = verbose);
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            flattened_dna = Vector{Vector{String}}()

            # Flatten using .seq field from SeqWithParent
            for (swp, count) in abund_dict_nucleo
                append!(flattened_dna, fill(swp.seq, count))
            end
            
            mat_dna = hcat(flattened_dna...)
            println(size(mat_dna))
            
            push!(post_sel_step_msa_dna, mat_dna)
            
            mat = zeros(Int8, size(mat_dna))
            
            for i in 1:size(mat_dna,1)
                for m in 1:size(mat_dna,2)
                    mat[i,m] = cod2amino[mat_dna[i,m]]
                end
            end
           
            push!(post_sel_step_msa, mat)
        end
        
        clean_dict_nucleo!(abund_dict_nucleo)
        Genie.amplification_nucleo_mcmc_inspired!(final_msa, final_msa_dna, abund_dict_nucleo)  
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            count+=1
            idxs = rand(1:M,seq_reads)
            step_msa[count] .= final_msa[:,idxs]
            step_msa_dna[count] .= final_msa_dna[:,idxs]
        end            
    end
    
    if (seq_steps !== nothing) || (each_step !== nothing)
        return (start_msa = start_msa, steps = steps, step_msa = step_msa, pre_sel_step_msa = pre_sel_step_msa, pre_sel_step_msa_dna = pre_sel_step_msa_dna, post_sel_step_msa = post_sel_step_msa, post_sel_step_msa_dna = post_sel_step_msa_dna, step_msa_dna = step_msa_dna, mu = mu, temp = temp)
    else
        return (start_msa = start_msa, rounds = rounds, final_msa = final_msa, final_msa_dna = final_msa_dna, mu = mu, temp = temp)
    end
end


function run_dir_evol_nucleo_mcmc_inspired(start_seq::Array{Int8,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
                   q::Int = 21,
                   verbose = false) where {T}
    
    
    start_msa = hcat([start_seq for i in 1:N_start]...);
    start_dna = amino_seq2dna_seq(start_seq);
    start_msa_dna = hcat([start_dna for i in 1:N_start]...);
    return run_dir_evol_nucleo(start_msa, start_msa_dna, h, J,
                   rounds = rounds,
                   each_step = each_step,
                   seq_reads = seq_reads,
                   temp = temp,  
                   mu = mu,
                   codon_bias = codon_bias,
                   mut_bias = mut_bias,
                   q = q,
                   verbose = verbose) 
end

function run_dir_evol_nucleo_mcmc_inspired(start_seq::Array{String,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
                   q::Int = 21,
                   verbose = false) where {T}
    
    
    
    start_amino_seq = [cod2amino[x] for x in start_seq];
    start_msa = hcat([start_amino_seq for i in 1:N_start]...);
    start_msa_dna = hcat([start_seq for i in 1:N_start]...);
    return run_dir_evol_nucleo(start_msa, start_msa_dna, h, J,
                   rounds = rounds,
                   each_step = each_step,
                   seq_reads = seq_reads,
                   temp = temp,  
                   mu = mu,
                   codon_bias = codon_bias,
                   mut_bias = mut_bias,
                   q = q,
                   verbose = verbose) 
end


