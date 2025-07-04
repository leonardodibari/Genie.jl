
function replace_char(s::String, i::Int, new_char::Char)
    return s[1:i-1] * string(new_char) * s[i+1:end]
end

function random_mutations_nucleo!(msa::Matrix{Int8}, msa_dna::Matrix{String}, 
        codon_net::Dict{String, Dict{Int, Vector{String}}}, 
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        codon_list::Vector{Vector{String}},
        abund_dict::Dict{Vector{String}, Int64}, 
        mu::T; 
        q::Int = 4,
        verbose = false) where {T}
    
    L,M = size(msa)
    distrib = Binomial(L, prob_mut_one_site(mu, q = q))
       
    all_nucleos = ['A','C', 'G', 'T'];
    @tasks for m in 1:M
        K = rand(distrib) ## number of mutation 
        pos = sample(1:3*L,K, replace = false); ## positions of mutations
        for k in 1:K
            seq_site = div(pos[k]-1,3) + 1
            nucleo_site = mod(pos[k]-1,3) + 1
            
            
            #=
            
            ##New version
            new_nucleo = sample(filter(x -> x != msa_dna[seq_site,m][nucleo_site], all_nucleos), 1)
            
            new_codon = replace_char(msa_dna[seq_site, m], nucleo_site, new_nucleo[1])
            if new_codon != "TAG" && new_codon != "TAA" && new_codon != "TGA"
                msa_dna[seq_site, m] = new_codon
            end
            msa[seq_site,m] = cod2amino[msa_dna[seq_site,m]] =#
            
            
            if msa[seq_site,m] == 21
                
            else
                
            ##Old Version
                L_moves = length_of_moves[msa_dna[seq_site,m], nucleo_site]
                if L_moves == 1
                    msa_dna[seq_site,m] = accessible_codons(msa_dna[seq_site,m], codon_net, nucleo_site)[1]
                    msa[seq_site,m] = cod2amino[msa_dna[seq_site,m]]
                else
                    codon_list[L_moves] .= accessible_codons(msa_dna[seq_site,m], codon_net, nucleo_site)         
                    msa_dna[seq_site,m] = sample(codon_list[L_moves])
                    msa[seq_site,m] = cod2amino[msa_dna[seq_site,m]]
                end
            end
            
        end
    end
       
    empty!(abund_dict)
    seq_count_dict_nucleo!(abund_dict, msa_dna)
    clean_dict_nucleo!(abund_dict)

    if verbose == true
        println("$(100*length(abund_dict)/M) % unique seqs after mutagenesis")
    end
    
end



function random_mutations_nucleo_biased!(msa::Matrix{Int8}, msa_dna::Matrix{String}, 
        codon_net::Dict{String, Dict{Int, Vector{String}}}, 
        length_of_moves::Dict{Tuple{String, Int64}, Int64},
        codon_list::Vector{Vector{String}},
        prob_codon_list::Vector{Vector{Float64}},
        abund_dict::Dict{Vector{String}, Int64}, 
        mu::T,
        trans_dict::Dict{Tuple{String, String}, Float64}; 
        q::Int = 4,
        verbose = false) where {T}
    
    L,M = size(msa)
    distrib = Binomial(L, prob_mut_one_site(mu, q = q))
    
    @tasks for m in 1:M
        K = rand(distrib) 
        pos = sample(1:3*L,K, replace = false);
        for k in 1:K
            seq_site = div(pos[k]-1,3) + 1
            nucleo_site = mod(pos[k]-1,3) + 1
            L_moves = length_of_moves[msa_dna[seq_site,m], nucleo_site]
            if L_moves == 1
                msa_dna[seq_site,m] = accessible_codons(msa_dna[seq_site,m], codon_net, nucleo_site)[1]
                msa[seq_site,m] = cod2amino[msa_dna[seq_site,m]]
            else
                codon_list[L_moves] .= accessible_codons(msa_dna[seq_site,m], codon_net, nucleo_site)
                for i in 1:length(codon_list[L_moves])
                    prob_codon_list[L_moves][i] = trans_dict[msa_dna[seq_site,m], codon_list[L_moves][i]]
                end                    
                msa_dna[seq_site,m] = sample(codon_list[L_moves], Weights(prob_codon_list[L_moves]))
                msa[seq_site,m] = cod2amino[msa_dna[seq_site,m]]
            end
        end
    end
    
    empty!(abund_dict)
    seq_count_dict_nucleo!(abund_dict, msa_dna)
    clean_dict_nucleo!(abund_dict)

    if verbose == true
        println("$(100*length(abund_dict)/M) % unique seqs after mutagenesis")
    end
end


function inv_degeneracy(seq::Array{String,1}, codon_usage::Dict{String, Float64})
    p = 1.
    for x in seq
        p *= codon_usage[x]
    end
    return p
end 
     
#### METTI 1/DEGENERAZIONE SEQUENZA DAVANTI A (a/1+a), forse idea non buona, circa 10^-100 per ogni sequenza

function selection_nucleo!(h::Array{T,2}, J::Array{T,4}, abund_dict::Dict{Vector{String}, Int64}, codon_usage::Dict{String, Float64}, L::Int, M::Int; temp::T = 1., mu_bind::T = 18.6, verbose = false) where {T}
        

    seq_dna = collect(keys(abund_dict))  # Turn to indexable vector

    aa = zeros(length(seq_dna))
    surv_prob = Vector{Float64}(undef, length(seq_dna))

    @tasks for i in 1:length(seq_dna)
        seq = seq_dna[i]
        if "TAG" in seq || "TAA" in seq || "TGA" in seq
            surv_prob[i] = 0.
        else
            aa[i] = exp( - ((energy_dna2(seq, h, J, L) - mu_bind ) / temp))
            p = 1. #inv_degeneracy(seq_dna[i], codon_usage)
            surv_prob[i] = (p*aa[i])/(1+aa[i])
        end
        distrib = Binomial(abund_dict[seq], surv_prob[i])
        abund_dict[seq_dna[i]] = rand(distrib)
    end
  
    
    #=
    surv_prob = [ ]
    for seq_dna in keys(abund_dict)
        if "TAG" in seq_dna || "TAA" in seq_dna || "TGA" in seq_dna
            push!(surv_prob,0.)
        else
            a = exp( - ((energy_dna2(seq_dna, h, J, L) - mu_bind ) / temp))
            p = 1. #inv_degeneracy(seq_dna, codon_usage)
            push!(surv_prob, (p*a)/(1+a))
        end
    end =#
        
    
    #println("Extrema of surv_prob is $(extrema(surv_prob))")
    
    #=idx = 1; 
    for seq in keys(abund_dict)
        distrib = Binomial(abund_dict[seq], surv_prob[idx])
        abund_dict[seq] = rand(distrib)
        idx+=1    
    end =#
    
    clean_dict_nucleo!(abund_dict)
    
    if verbose == true
        println("$(100*length(abund_dict)/M) % unique seqs after selection")
    end
end


function amplification_nucleo!(final_msa::Array{Int8,2}, final_msa_dna::Array{String,2}, abund_dict::Dict{Vector{String}, Int64})
    
    L, Mf = size(final_msa);
        
    seqs_dna = []
    for kk in keys(abund_dict)
        push!(seqs_dna, kk)
    end
    
    w = [abund_dict[seq] for seq in keys(abund_dict)];
    #w2 = [abund_dict[seq] for seq in seqs];
    
    M = length(seqs_dna);
    
    #println(sum(abs2, w .- w2))
    
    idxs = sample(1:M, Weights(w), Mf);
    
    @tasks for m in 1:Mf
        for i in 1:L
            final_msa_dna[i,m] = seqs_dna[idxs[m]][i]
            final_msa[i,m] = cod2amino[final_msa_dna[i,m]]
        end
    end
end


function seq_count_dict_nucleo!(dict::Dict{Vector{String}, Int64}, msa_dna::Array{String,2})

    M = size(msa_dna,2)
    
    for m=1:M
        haskey(dict,msa_dna[:,m]) ?  dict[msa_dna[:,m]]+=1 : dict[msa_dna[:,m]]=1
    end

    return dict

end



function clean_dict_nucleo!(d::Dict{Vector{String}, Int64})
    for k in keys(d)  # Iterate over keys
        if d[k] == 0
            delete!(d, k)  # Remove key with value 0
        end
    end
end
    


function run_dir_evol_nucleo(start_msa::Array{Int8,2}, start_msa_dna::Array{String, 2}, h::Array{T,2}, J::Array{T,4};
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
    count3 = 0
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
        
        Genie.selection_nucleo!(h, J, abund_dict_nucleo, codon_usage, L, M, temp = temp, mu_bind = mu_bind, verbose = verbose);
        
        
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
        return (start_msa = start_msa, rounds = rounds, final_msa = final_msa, final_msa_dna = final_msa_dna, mu = mu, mu_bind = mu_bind, temp = temp)
    end
end


function run_dir_evol_nucleo(start_seq::Array{Int8,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
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
    return run_dir_evol_nucleo(start_msa, start_msa_dna, h, J,
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

function run_dir_evol_nucleo(start_seq::Array{String,1}, N_start::Int, h::Array{T,2}, J::Array{T,4};
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
    return run_dir_evol_nucleo(start_msa, start_msa_dna, h, J,
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


