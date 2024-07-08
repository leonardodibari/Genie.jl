module Genie


using Revise
using OhMyThreads
using LinearAlgebra
using PyPlot
using DelimitedFiles
using Random
using JLD2
using StatsBase
using DCAUtils
using MultivariateStats

import KitMSA: fasta2matrix, letter2num, num2letter, extract_params, read_par_BM, set_max_field_to_0


include("gibbs.jl")
include("metropolis.jl")
include("utils.jl")
include("run.jl")
include("types.jl")
include("hamming.jl")
include("energy.jl")
include("equilibration.jl")


# Write your package code here.
export cod2amino, amino2cod, codon_usage, Chain, run_evolution, ham_dist, prob_cond!, log_prob!, log_prob2!
export random_gens, prob_cond2!, create_length_dict
export accessible_codons, energy, check_equilibration



end
