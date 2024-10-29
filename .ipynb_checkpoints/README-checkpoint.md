
# Genie.jl Evolution Simulation Package

**Genie.jl** is a Julia package designed to simulate the evolution of multiple sequence alignments (MSAs) based on specified parameters and conditions. Its primary function, `run_evolution`, enables users to simulate evolutionary trajectories for amino acid or nucleotide codon sequences under complex interaction parameters.

## Key Function: `run_evolution`

### Overview
The `run_evolution` function simulates the evolution of a given multiple sequence alignment (MSA) over a specified number of steps. It uses a combination of Gibbs sampling and Metropolis sampling to evolve the sequences, supporting options for random initialization, codon usage bias, and saving intermediate MSAs at specified intervals.

### Parameters
- **`start_msa::Array{T,2}`**: Initial MSA as a 2D array where `T` can be either `Int8` for amino acids or `String` for nucleotide codons. If amino acids are provided, the corresponding codons will be randomly sampled among those coding for the amino acids. The array dimensions must be `(L, M)`, where `L` is the sequence length, and each column represents a different sequence.

- **`h::Array{T,2}`**: A 2D array of size `(q, L)` representing the field parameters.

- **`J::Array{T,4}`**: A 4D array of size `(q, L, q, L)` representing the coupling parameters.

### Optional Parameters
- **`N_steps::Int`**: Number of steps for the simulation (default is 100).
- **`temp`**: Temperature parameter for the simulation (default is 1.0).
- **`p`**: Probability for choosing Metropolis Sampling (default is 0.5). A value of `0` will use only Gibbs Sampling with single nucleotide mutations, while `1` will use only Metropolis Sampling with indels.
- **`N_points::Union{Int, Nothing}`**: Number of points to save the MSA in logarithmic scale along the trajectory (default is `nothing`). Specify either `N_points` or `each_step`, but not both.
- **`each_step::Union{Int, Nothing}`**: Interval to save the MSA every `each_step` steps along the trajectory (default is `nothing`). Specify either `N_points` or `each_step`, but not both.
- **`rand_init::Bool`**: Whether to initialize sequences randomly (default is `false`).
- **`q::Int`**: Number of unique amino acids in the sequences (default is 21).
- **`codon_bias::Union{Nothing, Dict{String, Float64}}`**: Codon usage bias dictionary (default is `nothing`; assumes no codon bias).
- **`verbose::Bool`**: Whether to print progress information (default is `false`).

### Returns
The function returns a named tuple containing the results of the simulation. The structure of the output depends on whether `N_points` or `each_step` is specified.

1. **If `N_points` or `each_step` are not specified**:
   - **`msa::Array{Int8, 2}`**: Final MSA in amino acid form.
   - **`msa_dna::Array{String, 2}`**: Final MSA in DNA format.
   - **`codon_usage::Dict{String, Float64}`**: Codon usage dictionary used in the simulation.
   - **`p::Float64`**: Probability of choosing Metropolis Sampling.
   - **`temp::Float64`**: Temperature used in the simulation.

2. **If either `N_points` or `each_step` are specified**:
   - **`step_msa::Array{Array{Int8, 2}, 1}`**: List of MSAs at different time points in amino acid format.
   - **`msa_dna::Array{Array{String, 2}, 1}`**: List of MSAs in DNA format at different time points.
   - **`codon_usage::Dict{String, Float64}`**: Codon usage dictionary used in the simulation.
   - **`p::Float64`**: Probability of choosing Metropolis Sampling.
   - **`temp::Float64`**: Temperature used in the simulation.
   - **`steps::Array{Int, 1}`**: Steps at which MSAs were saved.

### Description
The `run_evolution` function simulates the evolution of an initial MSA over a specified number of steps. It employs both Gibbs sampling and Metropolis sampling to generate evolved sequences, with options for random initialization, codon usage bias, and saving intermediate MSAs at specified intervals. The function is highly customizable and suitable for simulating complex evolutionary dynamics under various conditions.

