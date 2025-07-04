module Genie


using Revise
using OhMyThreads
using LinearAlgebra
using PyCall
using PyPlot
using DelimitedFiles
using Random
using JLD2
using StatsBase
using DCAUtils
using MultivariateStats
using KernelDensity
using GZip
using Distributions
using TreeTools
using FastaIO
using LsqFit

import Flux: relu

include("gibbs.jl")
include("metropolis.jl")
include("utils.jl")
include("run.jl")
include("run_ontree.jl")
include("types.jl")
include("hamming.jl")
include("energy.jl")
include("equilibration.jl")
include("PCA_APM.jl")
include("timescales.jl")
include("entropy.jl")
include("g4analysis.jl")
include("read_write.jl")
include("random_potts.jl")
include("lattice_glass.jl")
include("directed_evolution.jl")
include("metropolis_directed_evolution.jl")
include("analyze_msa.jl")
include("directed_evolution_nucleo.jl")
include("directed_evolution_nucleo_mcmc_inspired.jl")
include("boltz_directed_evolution_nucleo.jl")
include("directed_evolution_nucleo_variant.jl")
include("mut_algo_for_nobu.jl")
include("mut_algo_for_nobu_beyond_revers.jl")
include("global_mut_algo.jl")

# Write your package code here.
export cod2amino, amino2cod, codon_usage, Chain, run_evolution, ham_dist, pairwise_ham_dist, prob_cond!, log_prob!
export random_gens, create_length_dict
export accessible_codons, energy, check_equilibration, perform_pca, G4, G4_nat, cont_dep_entr, CIE
export read_dist_from_file, run_potts, Potts_fixed_mean_conn, exchange_parameters

export neighbours, node_energy, lattice_energy, lattice_delta_energy, swap_occupation!, initialize_lattice
export swap_options, run_lattice_dynamics, ham_dist_AB_rp, ham_dist_AB_rp3, G4_AB, Seq, run_evolution_ontree, check_pairwise, msa_from_leafs, msa_from_nodes, pair_dist_freq, find_optimal_mu, check_energy, check_pca

export fasta2matrix, letter2num, num2letter, extract_params, read_par_BM, vec2string

export analyze_msa, analyze_target_msa, read_par_BM_0gapsave, read_par_BM_lettersave, symmetrize_Jsave, set_max_field_to_0save

end
