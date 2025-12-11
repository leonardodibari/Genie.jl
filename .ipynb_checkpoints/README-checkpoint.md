
# Genie.jl Evolution Simulation Package

**Genie.jl** is a Julia package designed to simulate the evolution of multiple sequence alignments (MSAs) based on specified parameters and conditions. Its primary function, `run_evolution`, enables users to simulate evolutionary trajectories for amino acid or nucleotide codon sequences under complex interaction parameters.

The complete description of the algorithm is available at
> Emergent time scales of epistasis in protein evolution   
> Leonardo Di Bari, Matteo Bisardi, Sabrina Cotogno, Martin Weigt, Francesco Zamponi;
> doi: https://www.pnas.org/doi/10.1073/pnas.2406807121

Please cite this article if you use Genie.jl! 

# Genie.jl Evolution Simulation Package

**Genie.jl** is a Julia package designed to simulate the evolution of multiple sequence alignments (MSAs) based on specified parameters and conditions. Its primary function, `run_evolution`, enables users to simulate evolutionary trajectories for amino acid or nucleotide codon sequences under complex interaction parameters.

# Installing the Genie Package via Git Clone

To install the Genie package directly from its Git repository, follow these steps:

Open a terminal and run the following command to clone the Genie repository:

```bash
git clone https://github.com/leonardodibari/Genie.jl.git
```

then enter in the Genie.jl folder and call Julia with

```bash
julia --project
```

once inside the Julia REPL install the dependencies


```julia
using Pkg
Pkg.add("Revise")
Pkg.add(url = "https://github.com/matteobisardi/KitMSA.git")
Pkg.activate(".")
Pkg.instantiate()
```


## Running MCMC from the Terminal

Genie.jl can be run directly from the terminal using the `launch_mcmc.jl` script. This is useful if you want to run simulations with independent markov chain monte carlowithout opening Julia REPL or Jupyter notebooks.

### Quick Start (Beginner-Friendly)

```bash
julia --project=. --threads n  scripts/launch_mcmc.jl \
    -wt_name TEM \
    -Nseqs 1000 \
    -Nsteps 10000 \
    -Npairs 100 \
    -temp 1.0 \
    -file_wt "../TEM1_pfam_DNA.fasta" \
    -file_model "../model_natural_PSE.dat" \
    -output_folder "../"
```


### What it does

- Reads the command-line arguments.  
- Loads the wild-type sequences and model parameters.  
- Performs the evolution simulation using MCMC.  
- Saves results (MSAs, energys, hamming from wildtype and pairwise hamming) in the specified output folder.

### Main Options

| Option | Description |
|--------|-------------|
| `-wt_name` | Name of the wild-type sequence. |
| `-Nseqs` | Number of sequences to simulate. |
| `-Nsteps` | Number of evolutionary MCMC steps. |
| `-Npairs` | Number of sequence pairs to consider for pairwise hamming computations. |
| `-temp` | Temperature parameter for the simulation (inverse selective pressure). |
| `-file_wt` | Path to the wild-type FASTA file (DNA format). |
| `-file_model` | Path to the precomputed DCA model parameters. |
| `-output_folder` | Directory where results will be saved. |

### Parallel Execution Note (`--threads n`)

The Julia command includes the flag `--threads n`, where `n` is the number of parallel threads. This flag is crucial for performance in MCMC simulations because it enables the simulation of the specified `Nseqs` (number of sequences) as independent Markov chains running in parallel. Utilizing multiple CPU cores via this flag significantly reduces the total wall time required for the simulation. For optimal speed, it is recommended to set `n` as equal or lower than to the number of available physical or logical cores on your machine.


## Running PopGen Simulations from the Terminal

Genie.jl can be run directly from the terminal using the `popgen` workflow. This allows simulating sequence evolution under population genetics parameters without opening Julia REPL or Jupyter notebooks.

### Quick Start (Beginner-Friendly)

```bash
julia --project=. --threads n scripts/launch_popgen.jl \
    -wt_name TEM \
    -Nseqs 100 \
    -rounds 12 \
    -mu 0.025 \
    -mu_bind 100.0 \
    -Npairs 100 \
    -temp 1.0 \
    -file_wt "../../data/TEM1_pfam_DNA.fasta" \
    -file_model "../../data/model_natural_PSE.dat" \
    -output_folder "../TEM/popgen/"
```


### What it does

- Reads the command-line arguments.  
- Loads the wild-type sequences and model parameters.  
- Performs population-genetics evolution using `run_dir_evol_nucleo`.  
- Saves results including:
  - Final MSAs  
  - Sequence energies  
  - Hamming distances from the wild-type  
  - Pairwise Hamming distances  

### Main Options

| Option | Description |
|--------|-------------|
| `-wt_name` | Name of the wild-type sequence. |
| `-Nseqs` | Number of sequences to simulate. |
| `-rounds` | Number of population-genetics rounds to perform. |
| `-mu` | Jukes cantor mutation rate. |
| `-mu_bind` | survival DCA score threshold. |
| `-Npairs` | Number of sequence pairs for pairwise Hamming computation. |
| `-temp` | Temperature parameter (inverse selective pressure). |
| `-file_wt` | Path to the wild-type FASTA file (DNA format). |
| `-file_model` | Path to the precomputed DCA model parameters. |
| `-output_folder` | Directory where results will be saved. |

## Running MCMC on a Phylogenetic Tree from the Terminal

Genie.jl can run MCMC evolution simulations directly along a phylogenetic tree using the `mcmctree` workflow. This allows simulating sequence evolution under tree constraints without opening Julia REPL or notebooks.

### Quick Start (Beginner-Friendly)

```bash
julia --project=. --threads n scripts/launch_mcmctree.jl \
    -wt_name TEM \
    -Nsteps 100 \
    -Npairs 100 \
    -temp 1.0 \
    -file_wt "../../data/TEM1_pfam_DNA.fasta" \
    -file_model "../../data/model_natural_PSE.dat" \
    -file_tree "../../data/TreeSingleThread_TEM1_round12.nwk" \
    -output_folder "../TEM/mcmctree/"
```

### What it does

- Reads the command-line arguments.  
- Loads the wild-type sequences, model parameters, and the phylogenetic tree.  
- Performs MCMC evolution along the tree using `run_evolution_ontree`.  
- Saves results including:
  - Final MSAs  
  - Sequence energies  
  - Hamming distances from the wild-type  
  - Pairwise Hamming distances  

### Main Options

| Option | Description |
|--------|-------------|
| `-wt_name` | Name of the wild-type sequence. |
| `-Nsteps` | Number of evolutionary steps along the tree. |
| `-Npairs` | Number of sequence pairs for pairwise Hamming computation. |
| `-temp` | Temperature parameter (inverse selective pressure). |
| `-file_wt` | Path to the wild-type FASTA file (DNA format). |
| `-file_model` | Path to the precomputed DCA model parameters. |
| `-file_tree` | Path to the phylogenetic tree file (Newick format). |
| `-output_folder` | Directory where results will be saved. |




# Running the Genie Package: Using the Example Notebook or Julia REPL

Once you have installed the Genie package, you can either use the example notebook provided in the `examples` folder or work directly from the Julia REPL with parallel processing.

## Option 1: Using the Example Notebook
- **Navigate to the `examples` Folder**:
   Locate the `examples` folder inside the Genie package directory. This folder contains an example Jupyter notebook designed to help you explore Genie’s features.

## Option 2: working from the Julia REPL
- **Navigate to the Genie package folder**:
   Open Julia in the local environment with **n threads** over which the MCMC samoling can be parallelized by doing

```bash   
../julia-1.10.0/bin/julia --project=. --thread n
```
- then you can directly copy each cell of the 'examples' folder in your terminal to explore Genie's features (remember to adjust the paths of the data).


## Additional details on key Functions: `run_evolution`, 

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

