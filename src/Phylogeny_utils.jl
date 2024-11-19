aminoalphabet = Dict{String,Int64}()
nucleoalphabet = Dict{String,Int64}()
aminoalphabet = Dict('A' => 1, 'B' => 21, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 'G' => 6, 'H' => 7, 'I' => 8, 'J' => 21, 'K' => 9, 'L' => 10, 'M' => 11, 'N' => 12, 'O' => 21, 'P' => 13, 'Q' => 14, 'R' => 15, 'S' => 16, 'T' => 17, 'U' => 21, 'V' => 18, 'W' => 19, 'X' => 21, 'Y' => 20, 'Z' => 21, '-' => 21)
nucleoalphabet = Dict('A' => 1, 'C' => 2, 'G' => 3, 'U' => 4, '-' => 5)

aminonum = Dict{Int64, String}()
nucleonum = Dict{Int64, String}()
aminonum = Dict(1 => 'A', 2 => 'C', 3 => 'D', 4 => 'E', 5 => 'F', 6 => 'G', 7 => 'H', 8 => 'I', 9 => 'K', 10 => 'L', 11 => 'M', 
                12 => 'N', 13 => 'P', 14 => 'Q', 15 => 'R', 16 => 'S', 17 => 'T', 18 => 'V', 19 => 'W', 20 => 'Y', 21 => '-')
nucleonum = Dict(1 => 'A', 2 => 'C', 3 => 'G', 4 => 'U', 5 => '-')


codontoamino = Dict("AAA" => 'K', "AAC" => 'N', "AAG" => 'K', "AAT" => 'N',
                    "ACA" => 'T', "ACC" => 'T', "ACG" => 'T', "ACT" => 'T',
                    "AGA" => 'R', "AGC" => 'S', "AGG" => 'R', "AGT" => 'S',
                    "ATA" => 'I', "ATC" => 'I', "ATG" => 'M', "ATT" => 'I',
                    "CAA" => 'Q', "CAC" => 'H', "CAG" => 'Q', "CAT" => 'H',
                    "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCT" => 'P',
                    "CGA" => 'R', "CGC" => 'R', "CGG" => 'R', "CGT" => 'R',
                    "CTA" => 'L', "CTC" => 'L', "CTG" => 'L', "CTT" => 'L',
                    "GAA" => 'E', "GAC" => 'D', "GAG" => 'E', "GAT" => 'D',
                    "GCA" => 'A', "GCC" => 'A', "GCG" => 'A', "GCT" => 'A',
                    "GGA" => 'G', "GGC" => 'G', "GGG" => 'G', "GGT" => 'G',
                    "GTA" => 'V', "GTC" => 'V', "GTG" => 'V', "GTT" => 'V',
                    "TAC" => 'Y', "TAT" => 'Y',
                    "TCA" => 'S', "TCC" => 'S', "TCG" => 'S', "TCT" => 'S',
                    "TGC" => 'C', "TGG" => 'W', "TGT" => 'C',
                    "TTA" => 'L', "TTC" => 'F', "TTG" => 'L', "TTT" => 'F',
                    "---" => '-')

function translate_codontoamino(nucleosequence::String)   
    L = length(nucleosequence)
    if L % 3 == 0
        AminoLength = L ÷ 3
    else
        println("Length is not a multiple of 3")
        return 0
    end
    
    codons = [nucleosequence[n*3+1:n*3+3] for n in 0:AminoLength-1]
    
    aminostring = ""
    for i in 1:AminoLength
        aminostring *= codontoamino[codons[i]]
    end

    return aminostring

end


function H_distance(seq1::Vector, seq2::Vector)
    ### Function computing the Hamming distance between two vectors

    L = length(seq1)
    @assert L == length(seq2) "Error: the two vectors do not have the same length"
    d = 0
    for i in 1:L
        if seq1[i] != seq2[i]
            d += 1
        end
    end
    return d
end


function H_distance_nogaps(seq1::Vector, seq2::Vector)
    ### Function computing the Hamming distance between two vectors without considering gaps

    L = length(seq1)
    @assert L == length(seq2) "Error: the two vectors do not have the same length"
    d = 0
    for i in 1:L
        if seq1[i] != 21 && seq2[i] != 21
            if seq1[i] != seq2[i]
                d += 1
            end
        end
    end
    return d
end


function read_fasta(filename::AbstractString, max_gap_fraction::Real, theta::Any, remove_dups::Bool)
    """
    filename: name of the fasta file 
    max_gap_fraction : maximum fraction of gaps per sequence
    theta : parameter of the reweigthing# function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)
    remove_dups : remove duplicate sequences
    """
    Z = read_fasta_alignment(filename, max_gap_fraction) #reads file and traslate it immediately to numbers (but the gap is at 21 (?)) - function from the package DCAutils
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z) #removes duplicate sequences - function from the package DCAutils
    end
    N, M = size(Z)
    q = round(Int, maximum(Z)) #number of symbols, I think
    W, Meff = compute_weights(Z, theta) #computes the weights of each sequence based on similarity with the others, and gives the effective number of sequences Meff
    return W, Z, N, M, q
end


function read_fasta_dict(FileSeq::String)
    """
    Reads the sequences in the input file and saves them in a dictionary with the name of the sequences as keys and the sequences (in numeric form) as values
    """
    seq_text = readlines(FileSeq)
    Z_text = Dict()
    seq = ""
    name = ""
    for i in 1:length(seq_text)
        if seq_text[i][1] == '>'
            name = split(seq_text[i],">")[2]
            seq = ""
        else
            seq *= seq_text[i]
        end
        Z_text[name] = seq
    end

    Z = Dict{}()
    for key in keys(Z_text)
        num_seq = [aminoalphabet[Z_text[key][i]] for i in 1:length(Z_text[key])]
        # if !(num_seq in values(Z))
        #     Z[key] = num_seq
        # end
        Z[key] = num_seq
    end
    return Z
end


function compute_empirical_freqs(Z::Array{Int64,2}, q::Int; eps::Float64 = 0.0) 
    """
    Function to compute the frequency of occurrence of each aminoacid in the MSA.

    parameters:
    - Z: MSA in format Length_of_sequences x Number_of_sequences
    - q: length of the alphabet
    - eps: pseudocount value

    output:
    - f: matrix Length_of_sequences x q with the frequence of each aminoacid at each site

    """
    L, M = size(Z)
    f = zeros(L, q)
    for i in 1:L
        for s in 1:M
            f[i, Z[i, s]] += 1
        end
    end
    f ./= M

    for s in 1:L
        for a in 1:q
            f[s,a] += eps
        end
        f[s,:] ./= sum(f[s,:])
    end
    return f
    
end

#################################################################
############### FELSENSTEIN ALGORITHM ###########################
#################################################################

### Define the type of data that I want to assign to the tree nodes. In my case it's just a vector of probabilities. Can it be done better?
Base.@kwdef mutable struct ProbabilityOnTree <: TreeNodeData # Create a custom data type
    prob :: Array{Float64,1} = zeros(Float64,21)
end


function FelsensteinSingle!(T, n, m, p, q)
"""
Function to run the Felsenstein algorithm on just one site. 
Input:
    - T: Tree we want to run over (using Pierre's TreeTools package)
    - n: Node we want to infer (can also be an inside node)
    - m: Value of μ (timescale of evolution)
    - p: Stationary probability on the leaves
    - q: Alphabet length (21 in proteins) 
Output: 
    - T: Modified tree with inferred probability for every amino acid at each internal node
    - NodeProb: Probability of each aminoacid on the requested node
"""
    
    LogNodeProb = zeros(Float64,q)
    for child in children(T[n])
        if isleaf(child)
            ChildProb = data(child).prob
        else
            ChildProb = FelsensteinSingle!(T, label(child), m, p, q) # Iterate over the child
        end
        t = branch_length(child)
        if t == 0
            t += 10^-5 # ADD CORRECTION FOR 0 BRANCH LENGTH CASE. TO BE MODIFIED CORRECTLY
        end 
        W = exp(- m * t) 
        TransitionProb = repeat(p', q) * (1-W) + W * I(q) # This is equivalent to P(b|a,t) = δ_{a,b} * exp(- μt) + (1-exp(-μt)) w(b)
        LogNodeProb += log.(TransitionProb * ChildProb) # This is equivalent to log(\sum_b P(b|a,t) ρ(b))
    end
    
    loc_softmax!(LogNodeProb)
    #NodeProb = exp.(LogNodeProb) ./ sum(exp.(LogNodeProb))
    data!(T[n], ProbabilityOnTree(prob = LogNodeProb))
    return LogNodeProb
end


function Felsenstein(Z::AbstractDict, S::Array{Float64,2}, T, n, m::Float64; q::Int = 21, verbose::Bool = true, eps::Float64 = 10^-5)
"""
Function to run the Felsenstein algorithm on a given tree. 
Input:
    - Z: Dictionary containing the sequences that occupy the leaves of the tree
    - S: Stationary distribution of aminoacid on each site (NOT obtained with Z, but from all the natural sequences, that may be different from Z)
    - T: Tree we want to run over (using Pierre TreeTools package)
    - n: Node we want to infer (can also be an inside node)
    - m: Value of μ (timescale of evolution)
    - q (optional): Alphabet length (21 in proteins)
    - verbose (optional): Prints each time a site is inferred (true)
    - eps (optional): Pseudocount value (10^-5) 
Output: 
    - Ancestor: Amino acid sequence that maximizes the probability of each site
    - AncestorProb: Probability of each amino acid on the requested node
"""
    
    @assert in(n,T) "The selected node is not present inside the tree"
    @assert abs(sum(S[1,:]) - 1) < 1e-8 "Stationary probability is not normalized, $(sum(S[1,:]))"
    @assert length(leaves(T)) == length(Z) "Number of leaves in the tree T does not match the number of sequences in Z"
    
    M = length(Z)
    L = length(collect(values(Z))[1])
    
    AncestorProb = zeros(Float64,(L,q))
    Ancestor = zeros(Int,L)
    @tasks for s in 1:L
        _T = copy(T)
        for leaf in leaves(_T)
            ρ = zeros(Float64,q)
            ρ[Z[label(leaf)][s]] = 1
            data!(leaf, ProbabilityOnTree(prob = ρ))
        end
        AncestorProb[s,:] = FelsensteinSingle!(_T, n, m, view(S,s,:), q)
        Ancestor[s] = argmax(view(AncestorProb,s,:))
        if verbose == true 
            println("Reconstruction of site $s completed, a[$s] = $(Ancestor[s])")
        end
    end

    return Ancestor, AncestorProb
end


function Felsenstein(FileNat::String, FileTree::String, FileSequences::String, m::Float64; verbose::Bool = true, eps::Float64 = 10^-5)
"""
Function to run the Felsenstein algorithm given the files containing the sequences and tree. 
Input:
    - FileNat: File containing the natural sequences (from which we estimate the stationary probability)
    - FileTree: File containing the tree
    - FileSequences: File containing the sequences that sit at the leaves of the tree
    - m: Value of μ (timescale of evolution)
    - verbose (optional): Prints each time a site is inferred (true)
    - eps (optional): Pseudocount value (10^-5) 
Output: 
    - AncestorSequence: Amino acid sequence that maximizes the probability of each site
    - AncestorProbability: Probability of each amino acid on the requested node
"""  
    _, NatMSA, _, _, q = read_fasta(FileNat,1.0,0.2,false)
    StationaryProb = compute_empirical_freqs(Int.(NatMSA), q; eps = eps)
    
    LeavesSequences = read_fasta_dict(FileSequences)
    MyTree = read_tree(FileTree, node_data_type = ProbabilityOnTree)

    AncestorSequence, AncestorProbability = Felsenstein(LeavesSequences, StationaryProb, MyTree, label(root(MyTree)), m; q = q, verbose = verbose, eps=eps)

    return AncestorSequence, AncestorProbability
end