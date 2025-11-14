using Genie, DelimitedFiles, DCAUtils, JLD2, PyPlot, Statistics, LinearAlgebra
import KitMSA: fasta2matrix, matrix2fasta, cod2amino, amino2cod
using TreeTools # Required for MCMCTree simulation
import KitMSA # Ensures functions like cod2amino are available

# Define a struct to hold all parameters read from the CSV
struct SimulationParams
    family_name::String
    file_model::String
    file_wt::String
    Nseqs::Int
    Npairs::Int
    
    # MCMC
    N_steps_mcmc::Int
    T_mcmc::Float64

    # MCMCTree
    N_steps_mcmct::Int
    T_mcmct::Float64
    file_tree_mcmct::String

    # PopGen
    rounds::Int
    temp_popgen::Float64
    mu_popgen::Float64
    mu_bind_popgen::Float64

    # Neutral
    mu_neutral::Float64
end

# Function to read parameters from the CSV file
function read_parameters(family_name::String, csv_file::String)::SimulationParams
    println("Reading parameters for $family_name from $csv_file...")
    
    # Read all lines from the CSV file
    lines = readlines(csv_file)
    header = split(lines[1], ',')
    
    # Find the correct row
    data_row = nothing
    for i in 2:length(lines)
        row = split(lines[i], ',')
        if row[1] == family_name
            data_row = Dict(zip(header, row))
            break
        end
    end

    if isnothing(data_row)
        error("Protein family '$family_name' not found in the parameter file.")
    end

    # Parse and construct the SimulationParams struct
    try
        return SimulationParams(
            data_row["Family_Name"],
            data_row["FILE_MODEL"],
            data_row["FILE_WT"],
            parse(Int, data_row["NSEQS"]),
            parse(Int, data_row["NPAIRS"]),
            parse(Int, data_row["N_STEPS_MCMC"]),
            parse(Float64, data_row["T_MCMC"]),
            parse(Int, data_row["N_STEPS_MCMCT"]),
            parse(Float64, data_row["T_MCMCT"]),
            data_row["FILE_TREE_MCMCT"],
            parse(Int, data_row["ROUNDS"]),
            parse(Float64, data_row["TEMP_POPGEN"]),
            parse(Float64, data_row["MU_POPGEN"]),
            parse(Float64, data_row["MU_BIND_POPGEN"]),
            parse(Float64, data_row["MU_NEUTRAL"])
        )
    catch e
        error("Error parsing parameters from CSV: $e")
    end
end

# Function to load common data (DCA model and WT sequence)
function setup_data(params::SimulationParams)


    println("Reading model parameters")
    if params.family_name == "AAC" || params.family_name == "DHFR"
            h_tmp, J_tmp = read_par_BM_0gapsave(file_model) 
    else
            h_tmp, J_tmp = read_par_BM_lettersave(file_model)
    end
    h = set_max_field_to_0save(h_tmp);
    J_tmp2 = symmetrize_Jsave(J_tmp); 
    J = permutedims(J_tmp2, [1,3,2,4]);

    # Reading wild-type (WT) sequence
    println("Reading wild-type sequence from $(params.file_wt)")
    wt = join(readdlm(params.file_wt, skipstart = 1));
    L_big = round(Int, length(wt) / 3)

    # Replicate the original DNA slicing logic
    if params.family_name == "AAC" || params.family_name == "DHFR"
        wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 1:L_big];
    else
        wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 3:L_big-1]; 
    end
    
    return h, J, wt_DNA_seq
end

# --- Simulation 1: MCMC ---
function run_mcmc(params::SimulationParams, h, J, wt_DNA_seq, output_folder::String)
    println("\n--- Running MCMC Simulation ---")
    
    Nsteps = params.N_steps_mcmc
    temp = params.T_mcmc
    
    
    # Performing evolution
    println("Performing evolution (MCMC)")
    start_msa = hcat([wt_DNA_seq for i in 1:params.Nseqs]...); 
    
    res = run_evolution(start_msa, 
        h, 
        J,
        p = 0., 
        temp = temp, 
        N_steps = Nsteps,  
        verbose = false);

    # Writing the result
    println("Writing the MCMC result in files")
    
    file_msa = output_folder*params.family_name*"_mcmc_samples"*string(params.Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).fa"
    file_en = output_folder*"en_"*params.family_name*"_mcmc_samples"*string(params.Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).txt"
    file_ham = output_folder*"ham_"*params.family_name*"_mcmc_samples"*string(params.Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).txt"
    
    KitMSA.matrix2fasta(file_msa, res.sequences)
    writedlm(file_en, res.energies)
    writedlm(file_ham, res.ham_dist)
    println("MCMC Simulation Complete.")
end

# --- Simulation 2: MCMCTree ---
function run_mcmctree(params::SimulationParams, h, J, wt_DNA_seq, output_folder::String)
    println("\n--- Running MCMCTree Simulation ---")
    
    Nsteps = params.N_steps_mcmct
    temp = params.T_mcmct
    file_tree = params.file_tree_mcmct
    
    # Calculate mu, which depends on the tree and Nsteps
    println("Reading the tree from $file_tree")
    
    # Using the first three functions from PhyloTools/TreeTools
    tree = read_tree(file_tree); 
    branch_d = mean([distance(tree.root, a) for a in leaves(tree)]);
    
    # L calculation (sequence length in amino acids)
    L_wt = [KitMSA.cod2amino[x] for x in wt_DNA_seq];
    L = length(L_wt);
    
    mu = (Nsteps / L) / branch_d;
    
    # Performing evolution
    println("Performing evolution (MCMCTree)")
    
    # Assuming run_evolution_ontree is available
    res = run_evolution_ontree(wt_DNA_seq, file_tree, 
        h, J, temp = temp, mu = mu, p = 0.); 
    msa = Int.(KitMSA.msa_from_leafs(res)); 
    
    # Writing the result
    println("Writing the MCMCTree result in files")
    
    file_msa = output_folder*params.family_name*"_mcmctree_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).fa"
    file_en = output_folder*"en_"*params.family_name*"_mcmctree_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).txt"
    file_ham = output_folder*"ham_"*params.family_name*"_mcmctree_T"*string(temp)* "_steps$(Nsteps)_pairs$(params.Npairs).txt"
    
    KitMSA.matrix2fasta(file_msa, msa)
    writedlm(file_en, res.energies)
    writedlm(file_ham, res.ham_dist)
    println("MCMCTree Simulation Complete.")
end

# --- Simulation 3: PopGen ---
function run_popgen(params::SimulationParams, h, J, wt_DNA_seq, output_folder::String)
    println("\n--- Running PopGen Simulation ---")
    
    rounds = params.rounds
    temp = params.temp_popgen
    mu = params.mu_popgen
    mu_bind = params.mu_bind_popgen
    
    # Performing evolution
    println("Performing evolution (PopGen)")
    
    # Assuming Genie.run_dir_evol_nucleo is available
    res = Genie.run_dir_evol_nucleo(wt_DNA_seq, params.Nseqs, h, J;
                       rounds = rounds,
                       seq_reads = params.Nseqs, 
                       temp = temp,  
                       mu = mu, 
                       mu_bind = mu_bind);

    # Writing the result
    println("Writing the PopGen result in files")
    
    file_msa = output_folder*params.family_name*"_popgen_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).fa"
    file_en = output_folder*"en_"*params.family_name*"_popgen_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).txt"
    file_ham = output_folder*"ham_"*params.family_name*"_popgen_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).txt"
    
    KitMSA.matrix2fasta(file_msa, res.sequences)
    writedlm(file_en, res.energies)
    writedlm(file_ham, res.ham_dist)
    println("PopGen Simulation Complete.")
end

# --- Simulation 4: Neutral ---
function run_neutral(params::SimulationParams, h, J, wt_DNA_seq, output_folder::String)
    println("\n--- Running Neutral Simulation ---")
    
    rounds = params.rounds
    # Following the original neutral script setup: 
    # Temp is set to 0.0 (no selection bias on energy), and mu_bind to 0.0.
    temp = 0.0 
    mu = params.mu_neutral
    mu_bind = 0.0 
    
    # Performing evolution
    println("Performing evolution (Neutral)")
    
    # Assuming Genie.run_dir_evol_nucleo is available
    res = Genie.run_dir_evol_nucleo(wt_DNA_seq, params.Nseqs, h, J;
                       rounds = rounds,
                       seq_reads = params.Nseqs, 
                       temp = temp,  
                       mu = mu, 
                       mu_bind = mu_bind,
                       neutral = true); # Key argument for neutral simulation

    # Writing the result
    println("Writing the Neutral result in files")
    
    file_msa = output_folder*params.family_name*"_neutral_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).fa"
    file_en = output_folder*"en_"*params.family_name*"_neutral_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).txt"
    file_ham = output_folder*"ham_"*params.family_name*"_neutral_samples"*string(params.Nseqs)*"_T"*string(temp)* "_mu$(mu)_mubind$(mu_bind)_pairs$(params.Npairs).txt"
    
    KitMSA.matrix2fasta(file_msa, res.sequences)
    writedlm(file_en, res.energies)
    writedlm(file_ham, res.ham_dist)
    println("Neutral Simulation Complete.")
end


# Main function to execute everything
function main()
    # Check for command line arguments (expecting the protein family name)
    if length(ARGS) < 1
        error("Usage: <ProteinFamilyName>. Please provide the family name (e.g., TEM).")
    end

    if length(ARGS) < 2
        error("Usage: Please provide the output_folder.")
    end

    family_name = ARGS[1]

    output_folder = ARGS[2]
    
    # The parameter file name from the user upload
    csv_file = "../best_pars.csv"

    # 1. Read parameters
    params = read_parameters(family_name, csv_file)

    # 2. Setup common data (DCA model h, J and WT sequence)
    h, J, wt_DNA_seq = setup_data(params)
    
    # 3. Run all simulations
    run_mcmc(params, h, J, wt_DNA_seq, output_folder)
    run_mcmctree(params, h, J, wt_DNA_seq, output_folder)
    run_popgen(params, h, J, wt_DNA_seq, output_folder)
    run_neutral(params, h, J, wt_DNA_seq, output_folder)
    
    println("\nAll four simulations completed successfully for family: $family_name.")
end

# Execute the main function
main()