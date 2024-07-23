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
using DensityPlot
using GZip
using Distributions

import KitMSA: fasta2matrix, letter2num, num2letter, extract_params, read_par_BM, set_max_field_to_0


include("gibbs.jl")
include("metropolis.jl")
include("utils.jl")
include("run.jl")
include("types.jl")
include("hamming.jl")
include("energy.jl")
include("equilibration.jl")
include("PCA_APM.jl")
include("timescales.jl")
include("entropy.jl")
include("g4analysis.jl")

include("random_potts.jl")


# Write your package code here.
export cod2amino, amino2cod, codon_usage, Chain, run_evolution, ham_dist, pairwise_ham_dist, prob_cond!, log_prob!
export random_gens, create_length_dict
export accessible_codons, energy, check_equilibration, perform_pca, G4, G4_nat, cont_dep_entr, CIE
export read_dist_from_file, run_potts, Potts_fixed_mean_conn, exchange_parameters



end
