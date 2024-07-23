struct Chain{Ti, T}
    seq::Array{Ti,1}
    DNA::Array{String,1}
    L::Int
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


function Chain(seq::Array{<:Integer,1}, DNA::Array{String,1}, q::Int, generator::Xoshiro; T::DataType = Float64, Ti::DataType=Int8)
    L = size(seq,1)
    seq = Ti.(seq)
    codon_list2=Vector{String}(undef, 2);codon_list3=Vector{String}(undef, 3);codon_list4 =Vector{String}(undef, 4);
    amino_list2=Ti.([1,1]); amino_list3=Ti.([1,1,1]); amino_list4=Ti.([1,1,1,1]);
    log_prob2 = T.([0.,0.]); log_prob3 = T.([0.,0.,0.]); log_prob4 = T.([0.,0.,0.,0.]);
    Chain{Ti, T}(seq, DNA, L, generator, codon_list2, codon_list3, codon_list4, 
        amino_list2, amino_list3, amino_list4, 
        log_prob2, log_prob3, log_prob4)
end

function Chain(seq::Array{<:Integer,1}, q::Int, generator::Xoshiro; T::DataType = Float64, Ti::DataType=Int8)
    L = size(seq,1)
    DNA = amino_seq2dna_seq(seq)
    seq = Ti.(seq)
    codon_list2=Vector{String}(undef, 2);codon_list3=Vector{String}(undef, 3);codon_list4 =Vector{String}(undef, 4);
    amino_list2=Ti.([1,1]); amino_list3=Ti.([1,1,1]); amino_list4=Ti.([1,1,1,1]);
    log_prob2 = T.([0.,0.]); log_prob3 = T.([0.,0.,0.]); log_prob4 = T.([0.,0.,0.,0.]);
    Chain{Ti, T}(seq, DNA, L, generator, codon_list2, codon_list3, codon_list4, 
        amino_list2, amino_list3, amino_list4, 
        log_prob2, log_prob3, log_prob4)
end

function Chain(DNA::Array{String,1}, q::Int, generator::Xoshiro; T::DataType = Float64, Ti::DataType=Int8)
    L = size(DNA,1)
    seq= Ti.([cod2amino[x] for x in DNA])
    codon_list2=Vector{String}(undef, 2);codon_list3=Vector{String}(undef, 3);codon_list4 =Vector{String}(undef, 4);
    amino_list2=Ti.([1,1]); amino_list3=Ti.([1,1,1]); amino_list4=Ti.([1,1,1,1]);
    log_prob2 = T.([0.,0.]); log_prob3 = T.([0.,0.,0.]); log_prob4 = T.([0.,0.,0.,0.]);
    Chain{Ti, T}(seq, DNA, L, generator, codon_list2, codon_list3, codon_list4, 
        amino_list2, amino_list3, amino_list4, 
        log_prob2, log_prob3, log_prob4)
end


mutable struct SeqToEvolve
    Amino   
    DNA 
end


struct AminoChain{Ti}
    seq::Array{Ti,1}
    L::Int
    generator::Xoshiro
end


function AminoChain(seq::Array{<:Integer,1}, q::Int, generator::Xoshiro; Ti::DataType=Int8)
    L = size(seq,1)
    seq = Ti.(seq)
    AminoChain{Ti}(seq, L, generator)
end

