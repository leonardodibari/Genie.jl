function selection_nucleo_boltz!(h::Array{T,2}, J::Array{T,4}, abund_dict::Dict{Vector{String}, Int64}, codon_usage::Dict{String, Float64}, L::Int, M::Int; temp::T = 1., mu_bind::T = 18.6, verbose = false) where {T}
        

    seq_dna = collect(keys(abund_dict))  # Turn to indexable vector

    aa = zeros(length(seq_dna))
    surv_prob = Vector{Float64}(undef, length(seq_dna))

    @tasks for i in 1:length(seq_dna)
        seq = seq_dna[i]
        if "TAG" in seq || "TAA" in seq || "TGA" in seq
            surv_prob[i] = 0.
        else
            surv_prob[i] = exp( - (energy_dna2(seq, h, J, L) / temp))
        end
        distrib = Binomial(abund_dict[seq], surv_prob[i])
        abund_dict[seq_dna[i]] = rand(distrib)
    end
    
    surv_prob ./= maximum(surv_prob)
    
    @tasks for i in 1:length(seq_dna)
        seq = seq_dna[i]
        distrib = Binomial(abund_dict[seq], surv_prob[i])
        abund_dict[seq_dna[i]] = rand(distrib)
    end
  

    println("Extrema of surv_prob is $(extrema(surv_prob))")
    
    clean_dict_nucleo!(abund_dict)
    
    if verbose == true
        println("$(100*length(abund_dict)/M) % unique seqs after selection")
    end
end

function run_dir_evol_nucleo_boltz(start_msa::Array{Int8,2}, start_msa_dna::Array{String, 2}, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mu_bind::Float64 = 18.6,
                   q::Int = 21,
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
                   mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
                   verbose = false) where {T}
    
    L,M = size(start_msa);
    
    
    if M < seq_reads
        error("Increase number of sequences, seq_reads must be smaller than amount of sequeces")
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

    all_codons = vcat([amino2cod[i] for i in 1:20]...)
    push!(all_codons, "TAG")
    push!(all_codons, "TAA")
    push!(all_codons, "TGA")
    final_msa = deepcopy(start_msa);
    final_msa_dna = deepcopy(start_msa_dna);
    abund_dict_nucleo = Dict{Vector{String},Int}()
    codon_list = [Vector{String}(undef, i) for i in 1:4]
    prob_codon_list = [Vector{Float64}(undef, i) for i in 1:4]

    if codon_bias == nothing
        no_cod_bias = Dict(x => T(1/length(amino2cod[cod2amino[x]])) for x in keys(cod2amino))
        codon_usage = no_cod_bias
    end
    
    if mut_bias !== nothing
        trans_dict = create_codon_transition_dict(mut_bias)
    end
        
    for r in 1:rounds
        if verbose == true
            println()
            println("Round $(r)")
        end
        
        if mut_bias == nothing 
            Genie.random_mutations_nucleo!(final_msa, final_msa_dna, codon_net, 
                length_of_moves, codon_list, abund_dict_nucleo, mu, verbose = verbose);
        else
            Genie.random_mutations_nucleo_biased!(final_msa, final_msa_dna, codon_net, 
                length_of_moves, codon_list, prob_codon_list, abund_dict_nucleo, mu, 
                trans_dict, verbose = verbose);
        end

        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            count2+=1
            idxs = rand(1:M,seq_reads)
            pre_sel_step_msa[count2] .= final_msa[:,idxs]
            pre_sel_step_msa_dna[count2] .= final_msa_dna[:,idxs]
        end 
        
        Genie.selection_nucleo_boltz!(h, J, abund_dict_nucleo, codon_usage, L, M, temp = temp, mu_bind = mu_bind, verbose = verbose);
        
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            flattened_dna = Vector{Vector{String}}()

            for (vec, count) in abund_dict_nucleo
                append!(flattened_dna, fill(vec, count))
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
        Genie.amplification_nucleo!(final_msa, final_msa_dna, abund_dict_nucleo)  
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            count+=1
            idxs = rand(1:M,seq_reads)
            step_msa[count] .= final_msa[:,idxs]
            step_msa_dna[count] .= final_msa_dna[:,idxs]
        end            
    end
    
    if (seq_steps !== nothing) || (each_step !== nothing)
        return (start_msa = start_msa, steps = steps, step_msa = step_msa, pre_sel_step_msa = pre_sel_step_msa, pre_sel_step_msa_dna = pre_sel_step_msa_dna, post_sel_step_msa = post_sel_step_msa, post_sel_step_msa_dna = post_sel_step_msa_dna, step_msa_dna = step_msa_dna, mu = mu, mu_bind = mu_bind, temp = temp)
    else
        return (start_msa = start_msa, rounds = rounds, final_msa = final_msa, mu = mu, mu_bind = mu_bind, temp = temp)
    end
end


function run_dir_evol_nucleo_boltz(start_seq::Array{Int8,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mu_bind::Float64 = 18.6,
                   mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
                   q::Int = 21,
                   verbose = false) where {T}
    
    
    start_msa = hcat([start_seq for i in 1:N_start]...);
    start_dna = amino_seq2dna_seq(start_seq);
    start_msa_dna = hcat([start_dna for i in 1:N_start]...);
    return run_dir_evol_nucleo_boltz(start_msa, start_msa_dna, h, J,
                   rounds = rounds,
                   each_step = each_step,
                   seq_reads = seq_reads,
                   temp = temp,  
                   mu = mu,
                   mu_bind = mu_bind,
                   codon_bias = codon_bias,
                   mut_bias = mut_bias,
                   q = q,
                   verbose = verbose) 
end

function run_dir_evol_nucleo_boltz(start_seq::Array{String,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mut_bias::Union{Nothing, Dict{Tuple{Char, Char}, Float64}} = nothing,
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing,
                   mu_bind::Float64 = 18.6,
                   q::Int = 21,
                   verbose = false) where {T}
    
    
    
    start_amino_seq = [cod2amino[x] for x in start_seq];
    start_msa = hcat([start_amino_seq for i in 1:N_start]...);
    start_msa_dna = hcat([start_seq for i in 1:N_start]...);
    return run_dir_evol_nucleo_boltz(start_msa, start_msa_dna, h, J,
                   rounds = rounds,
                   each_step = each_step,
                   seq_reads = seq_reads,
                   temp = temp,  
                   mu = mu,
                   mu_bind = mu_bind,
                   codon_bias = codon_bias,
                   mut_bias = mut_bias,
                   q = q,
                   verbose = verbose) 
end


