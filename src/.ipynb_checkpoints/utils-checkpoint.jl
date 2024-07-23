function loc_sample!(rng::AbstractRNG, wv, codon_list::Array{String,1}, dest_seq::Array{Ti, 1}, dest_DNA::Array{String,1}, seq_site::Int) where {Ti<:Integer}
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = one(Ti)
    cw = wv[1]
    while cw < t && i < n
        i += one(Ti)
        @inbounds cw += wv[i]
    end
    
    dest_DNA[seq_site] = codon_list[i]
    dest_seq[seq_site] = cod2amino[codon_list[i]]
end


function loc_softmax!(out::AbstractArray{T}, x::AbstractArray{T}) where {T}
    max_ = T(maximum(x))
    if isfinite(max_)
        @fastmath out .= exp.(x .- max_)
    else
        _zero, _one, _inf = T(0), T(1), T(Inf)
        @fastmath @. out = ifelse(isequal(max_,_inf), ifelse(isequal(x,_inf), _one, _zero), exp(x - max_))
    end
    #tmp = dims isa Colon ? sum(out) : sum!(max_, out)
    out ./= sum(out)#tmp
end

loc_softmax!(x::AbstractArray) = loc_softmax!(x, x)


compute_freq(Z::Matrix) = compute_weighted_frequencies(Matrix{Int8}(Z), fill(1/size(Z,2), size(Z,2)), 22)

function compute_freq!(f1, f2, Z::Matrix)
    f1 .= compute_weighted_frequencies(Matrix{Int8}(Z), fill(1/size(Z,2), size(Z,2)), 22)[1]
    f2 .= compute_weighted_frequencies(Matrix{Int8}(Z), fill(1/size(Z,2), size(Z,2)), 22)[2]
end



function random_gens(num_generators::Int) 
    rng_array = []
    for seed in 1:num_generators
        push!(rng_array, Random.Xoshiro(seed))
    end
    return rng_array
end


function quickread(fastafile; moreinfo=false)  
    Weights, Z, N, M, _ = ReadFasta(fastafile, 0.9, :auto, true, verbose = false);
    moreinfo && return Weights, Z, N, M
    return Matrix{Int8}(Z), Weights
end


function quickread(fastafile, n_seq::Int; moreinfo=false)  
    Weights, Z, N, M, _ = ReadFasta(fastafile, 0.9, :auto, true, n_seq, verbose = false);
    moreinfo && return Weights, Z, N, M
    return Matrix{Int8}(Z), Weights
end

function ReadFasta(filename::AbstractString,max_gap_fraction::Real, theta::Any, remove_dups::Bool;verbose=true)
    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z,verbose=verbose)
    end
    N, M = size(Z)
    q = round(Int,maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W , Meff = compute_weights(Z,q,theta,verbose=verbose)
    println("Meff = $(Meff)")
    rmul!(W, 1.0/Meff)
    Zint=round.(Int,Z)
    return W,Zint,N,M,q
end

function ReadFasta(filename::AbstractString,max_gap_fraction::Real, theta::Any, remove_dups::Bool, n_seq::Int;verbose=true)
    Z = read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z,verbose=verbose)
    end
    ZZ = Z[:, sample(1:size(Z,2), n_seq, replace=false, ordered=true)]
    N, M = size(Z)
    q = round(Int,maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W , Meff = compute_weights(ZZ,q,theta,verbose=verbose)
    println("Meff = $(Meff)")
    rmul!(W, 1.0/Meff)
    Zint=round.(Int,ZZ)
    return W,Zint,N,M,q
end


function amino_seq2dna_seq(seq)
    seq_dna=[]
    for a in seq
        push!(seq_dna, sample(amino2cod[a]))
    end
    return seq_dna
end

ecoli_usage = Dict(
    "TTT" => 0.58, "TTC" => 0.42,
    "TTA" => 0.14, "TTG" => 0.13,
    "TCT" => 0.17, "TCC" => 0.15, "TCA" => 0.14, "TCG" => 0.14,
    "TAT" => 0.59, "TAC" => 0.41,
    "TGT" => 0.46, "TGC" => 0.54, "TGG" => 1.00,
    "CTT" => 0.12, "CTC" => 0.10, "CTA" => 0.04, "CTG" => 0.47,
    "CCT" => 0.18, "CCC" => 0.13, "CCA" => 0.20, "CCG" => 0.49,
    "CAT" => 0.57, "CAC" => 0.43,
    "CAA" => 0.34, "CAG" => 0.66,
    "CGT" => 0.36, "CGC" => 0.36, "CGA" => 0.07, "CGG" => 0.1,
    "ATT" => 0.50, "ATC" => 0.39, "ATA" => 0.11,
    "ACT" => 0.19, "ACC" => 0.40, "ACA" => 0.16, "ACG" => 0.25,
    "AAT" => 0.49, "AAC" => 0.51,
    "AAA" => 0.74, "AAG" => 0.26,
    "AGT" => 0.16, "AGC" => 0.24, "AGA" => 0.07, "AGG" => 0.04,
    "GTT" => 0.28, "GTC" => 0.20, "GTA" => 0.17, "GTG" => 0.35,
    "GCT" => 0.18, "GCC" => 0.26, "GCA" => 0.23, "GCG" => 0.33,
    "GAT" => 0.63, "GAC" => 0.37,
    "GAA" => 0.68, "GAG" => 0.32, "ATG" => 1.00,
    "GGT" => 0.35, "GGC" => 0.37, "GGA" => 0.13, "GGG" => 0.15, "---" => 1.
)

cod2amino::Dict{String, Int8} = Dict(
    "ATA" => Int8(8), "ATC" => Int8(8), "ATT" => Int8(8), "ATG" => Int8(11),
    "ACA" => Int8(17), "ACC" => Int8(17), "ACG" => Int8(17), "ACT" => Int8(17),
    "AAC" => Int8(12), "AAT" => Int8(12), "AAA" => Int8(9), "AAG" => Int8(9),
    "AGC" => Int8(16), "AGT" => Int8(16), "AGA" => Int8(15), "AGG" => Int8(15),
    "CTA" => Int8(10), "CTC" => Int8(10), "CTG" => Int8(10), "CTT" => Int8(10),
    "CCA" => Int8(13), "CCC" => Int8(13), "CCG" => Int8(13), "CCT" => Int8(13),
    "CAC" => Int8(7), "CAT" => Int8(7), "CAA" => Int8(14), "CAG" => Int8(14),
    "CGA" => Int8(15), "CGC" => Int8(15), "CGG" => Int8(15), "CGT" => Int8(15),
    "GTA" => Int8(18), "GTC" => Int8(18), "GTG" => Int8(18), "GTT" => Int8(18),
    "GCA" => Int8(1), "GCC" => Int8(1), "GCG" => Int8(1), "GCT" => Int8(1),
    "GAC" => Int8(3), "GAT" => Int8(3), "GAA" => Int8(4), "GAG" => Int8(4),
    "GGA" => Int8(6), "GGC" => Int8(6), "GGG" => Int8(6), "GGT" => Int8(6),
    "TCA" => Int8(16), "TCC" => Int8(16), "TCG" => Int8(16), "TCT" => Int8(16),
    "TTC" => Int8(5), "TTT" => Int8(5), "TTA" => Int8(10), "TTG" => Int8(10),
    "TAC" => Int8(20), "TAT" => Int8(20), "TGC" => Int8(2), "TGT" => Int8(2), "TGG" => Int8(19),
    "---" => Int8(21)
)
#=
cod2amino = Dict( "ATA" => Int8(8), "ATC" => Int8(8), "ATT"=> Int8(8), "ATG"=> Int8(11), 
        "ACA"=>Int8(17), "ACC"=>Int8(17), "ACG"=>Int8(17), "ACT"=> Int8(17), 
        "AAC"=>Int8(12), "AAT"=>Int8(12), "AAA"=>Int8(9), "AAG"=>Int8(9), 
        "AGC"=>Int8(16), "AGT"=> Int8(16), "AGA"=> Int8(15), "AGG"=> Int8(15),                  
        "CTA"=>Int8(10), "CTC"=>Int8(10), "CTG"=>Int8(10), "CTT"=>Int8(10), 
        "CCA"=>Int8(13), "CCC"=>Int8(13), "CCG"=>Int8(13), "CCT"=>Int8(13), 
        "CAC"=>Int8(7), "CAT"=>Int8(7), "CAA"=>Int8(14), "CAG"=>Int8(14), 
        "CGA"=>Int8(15), "CGC"=>Int8(15), "CGG"=>Int8(15), "CGT"=>Int8(15), 
        "GTA"=>Int8(18), "GTC"=>Int8(18), "GTG"=>Int8(18), "GTT"=>Int8(18), 
        "GCA"=>Int8(1), "GCC"=>Int8(1), "GCG"=>Int8(1), "GCT"=>Int8(1), 
        "GAC"=>Int8(3), "GAT"=>Int8(3), "GAA"=>Int8(4), "GAG"=>Int8(4), 
        "GGA"=>Int8(6), "GGC"=>Int8(6), "GGG"=>Int8(6), "GGT"=>Int8(6), 
        "TCA"=>Int8(16), "TCC"=>Int8(16), "TCG"=>Int8(16), "TCT"=>Int8(16), 
        "TTC"=>Int8(5), "TTT"=>Int8(5), "TTA"=>Int8(10), "TTG"=>Int8(10), 
        "TAC"=>Int8(20), "TAT"=>Int8(20), "TGC"=> Int8(2), "TGT"=>Int8(2) , "TGG"=> Int8(19), 
        "---" => Int8(21) )
=#

amino2cod = Dict()
for amino in 1:21
    codons = Array{String, 1}()
    for (key, val) in cod2amino
       val == amino && push!(codons, key)
    end
    amino2cod[amino] = codons
end



function codon_neighbors()
    codons = ["A", "C", "G", "T"]
    codon_dict = Dict{String, Vector{String}}()

    for a in codons, b in codons, c in codons
        old_codon = string(a, b, c)
        if old_codon != "TAA" && old_codon != "TAG" && old_codon != "TGA" 
            codon_list = String[]
        
            for i in 1:3
                for nucl in codons
                    if nucl != old_codon[i]
                        new_codon = string(old_codon[1:i-1], nucl, old_codon[i+1:end])
                        if new_codon != "TAA" && new_codon != "TAG" && new_codon != "TGA" 
                            push!(codon_list, new_codon)
                        end
                    end
                end
            end
            codon_dict[old_codon] = unique(codon_list)
        end
    end

    return codon_dict
end


# Function to construct the nested codon dictionary
function create_nested_codon_dict()
    nucleotides = ["A", "C", "G", "T"]
    codon_dict = Dict{String, Dict{Int, Vector{String}}}()

    for a in nucleotides, b in nucleotides, c in nucleotides
        old_codon = string(a, b, c)
        if old_codon != "TAA" && old_codon != "TAG" && old_codon != "TGA"
            codon_dict[old_codon] = Dict{Int, Vector{String}}()
            for i in 1:3
                codon_list = String[]
                for nucl in nucleotides
                    if nucl != old_codon[i]
                        new_codon = string(old_codon[1:i-1], nucl, old_codon[i+1:end])
                        if new_codon != "TAA" && new_codon != "TAG" && new_codon != "TGA"
                            push!(codon_list, new_codon)
                        end
                    end
                end
                codon_dict[old_codon][i] = codon_list
            end
        end
    end

    return codon_dict
end


function create_length_dict(codon_net::Dict{String, Dict{Int64, Vector{String}}})
    length_dict = Dict{Tuple{String, Int64}, Int}()

    for (key1, sub_dict) in codon_net
        for (key2, vec) in sub_dict
            length_dict[(key1, key2)] = length(vec)
        end
    end

    return length_dict
end


function all_equal(arr::AbstractVector)
    # Check if the array is empty
    if isempty(arr)
        return true
    end
    # Get the first element to compare against
    first_element = arr[1]
    # Iterate through the array and compare each element
    for i in 2:length(arr)
        if arr[i] != first_element
            return false
        end
    end
    return true
end




function accessible_codons(old_codon::String, codon_net::Dict{String, Dict{Int64, Vector{String}}}, nucleo_pos::Int)
    codon_changes = codon_net[old_codon][nucleo_pos]
    return codon_changes
end


# Function to get accessible nucleotide mutations using the nested codon dictionary
function accessible_codons_old(old_codon::String, codon_net::Dict{String, Dict{Int64, Vector{String}}}, nucleo_pos::Int)
    
    return codon_net[old_codon][nucleo_pos], Int8.(get.(Ref(cod2amino), codon_net[old_codon][nucleo_pos], 0))
end



function read_dist_from_file(filename::String, zipped::Bool)

    if zipped
        filename = filename * ".gz" 
        file = GZip.open(filename,"r")
        dist_vec = Int.(readdlm(file))
        close(file)
    else
        file = open(filename,"r")
        dist_vec = Int.(readdlm(file))
        close(file)
    end

    return dist_vec
    
end



function exchange_parameters(J::Array{Float64, 4}, t::Float64, n::Int)
    q, L, q2, L2 = size(J)
    
    if q != q2 || L != L2
        error("The dimensions of J should be q x L x q x L")
    end
    
    J_copy = copy(J)
    
    # Step 1: Identify elements above the threshold
    indices_above_threshold = []
    for a in 1:q, i in 1:L, b in 1:q, j in 1:L
        if abs.(J[a, i, b, j]) > t
            push!(indices_above_threshold, (a, i, b, j))
        end
    end

    # Step 2: Select `n` elements randomly
    if length(indices_above_threshold) < n
        error("Not enough elements above the threshold to exchange.")
    end
    selected_indices = rand(indices_above_threshold, n)

    # Step 3: Exchange elements while maintaining symmetry
    for k in 1:2:n  # Step in pairs
        if k < length(selected_indices)
            # Get the pair of indices to swap
            idx1 = selected_indices[k]
            idx2 = selected_indices[k+1]

            # Perform the swap
            J_copy[idx1[1], idx1[2], idx1[3], idx1[4]], J_copy[idx2[1], idx2[2], idx2[3], idx2[4]] = 
            J_copy[idx2[1], idx2[2], idx2[3], idx2[4]], J_copy[idx1[1], idx1[2], idx1[3], idx1[4]]

            # Maintain symmetry
            J_copy[idx1[3], idx1[4], idx1[1], idx1[2]], J_copy[idx2[3], idx2[4], idx2[1], idx2[2]] = 
            J_copy[idx2[3], idx2[4], idx2[1], idx2[2]], J_copy[idx1[3], idx1[4], idx1[1], idx1[2]]
        end
    end

    return J_copy
end