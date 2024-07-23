# Genie


# Evolution Simulation Package

This package provides tools for simulating the evolution of multiple sequence alignments (MSAs) under different conditions. The main function in this package is `run_evolution`, which simulates the evolutionary process based on given parameters and options.

## Parameters

- **start_msa::Array{T,2}**: Initial multiple sequence alignment (MSA) as a 2D array, T can be both Int8 if it's composed of aminoacids or a String if it's composed of nucleotide codons. If you give aminoacids as inputs, the relative codons will be sampled randomly among the ones coding for the repective aminoacids. The format size must be (L,M), where L is the length of the sequence, so columns represent different sequences.
- **h::Array{T,2}**: 2D array of size (q, L) representing field parameters.
- **J::Array{T,4}**: 4D array of size (q, L, q, L) representing coupling parameters.

## Optional parameters
- **N_steps::Int**: Number of steps for the simulation (default is 100).
- **temp**: Temperature parameter for the simulation (default is 1.0).
- **p**: Probability for choosing Metropolis Sampling (default is 0.5, if p = 0 it only does single nucleotide mutations with Gibbs Sampling, if p = 1 it does only indels with Metropolis Sampling).
- **N_points::Union{Int, Nothing}**: Number of points to save MSA in logarithmic scale along the trajectory (default is nothing, you can put N_points or each_step, but not both).
- **each_step::Union{Int, Nothing}**: Interval to save MSA every each_step steps along the trajectory (default is nothing, you put N_points or each_step, but not both).
- **rand_init::Bool**: Whether to initialize sequences randomly (default is false).
- **q::Int**: Number of unique aminoacids in the sequences (default is 21).
- **codon_bias::Union{Nothing, Dict{String, Float64}}**: Codon usage bias dictionary (default is nothing, as default we assume no codon bias).
- **verbose::Bool**: Whether to print progress information (default is false).

## Returns

A named tuple containing the results of the simulation. If N_points or each_step remain to default nothing you get:

- **msa::Array{Int8, 2}**: Final multiple sequence alignment.
- **msa_dna::Array{String, 2}**: Final multiple sequence alignment in DNA format.
- **codon_usage::Dict{String, Float64}**: Codon usage used in the simulation.
- **p::Float64**: Probability of choosing Metropolis used for the simulation.
- **temp::Float64**: Temperature used for the simulation.

Instead, if either N_points or each_step get specified you get:

- **step_msa::Array{Array{Int8, 2},1}**: List of multiple sequence alignment at different time points.
- **msa_dna::Array{Array{String, 2},1}**: List multiple sequence alignment in DNA format at different time points.
- **codon_usage::Dict{String, Float64}**: Codon usage used in the simulation.
- **p::Float64**: Probability of choosing Metropolis used for the simulation.
- **temp::Float64**: Temperature used for the simulation.
- **steps::Array{Int, 1}** (if N_points or each_step is specified): Steps at which MSAs were saved.

## Description

The `run_evolution` function simulates the evolution of a given multiple sequence alignment (MSA) over a specified number of steps. It uses a combination of Gibbs sampling and Metropolis sampling to evolve the sequences. The function supports options for random initialization, codon usage bias, and saving intermediate MSAs at specified intervals.


[![Build Status](https://github.com/leonardodibari/Genie.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/leonardodibari/Genie.jl/actions/workflows/CI.yml?query=branch%3Amain)
