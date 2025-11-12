
using Genie, DelimitedFiles, DCAUtils, JLD2, PyPlot, Statistics, LinearAlgebra
import KitMSA: fasta2matrix, matrix2fasta

wt_name = "TEM"
Nseqs = 100
Nsteps = 100
Npairs = 100
temp = 1.0
file_wt = "../data/TEM1_pfam_DNA.fasta"
file_model = "../data/model_natural_PSE.dat"

println("Reading arguments")
L_args = length(ARGS)
@assert L_args > 0 "Please specify the info for the simulation"
for i in 1:L_args
    if ARGS[i] == "-wt_name"
        global wt_name = parse(String,ARGS[i+1])
    elseif ARGS[i] == "-Nseqs"
        global Nseqs = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-Nsteps"
        global Nsteps = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-Npairs"
        global Npairs = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-temp"
        global temp = parse(Float64,ARGS[i+1])     
    elseif ARGS[i] == "-file_wt"
        global file_wt = parse(String,ARGS[i+1])
    elseif ARGS[i] == "-file_model"
        global file_model = parse(Int,ARGS[i+1])       
    end
end


# Read model parameters
println("Reading model parameters")
h_tmp, J_tmp = read_par_BM_lettersave(file_model);
h = set_max_field_to_0save(h_tmp);
J_tmp2 = symmetrize_Jsave(J_tmp); 
J = permutedims(J_tmp2, [1,3,2,4]);


# Reading sequencies in the seed
println("Reading sequences in the seed")
wt = join(readdlm(file_wt, skipstart = 1));
L_big = round(Int,length(wt)/3)
wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 3:L_big-1];
wt_seq = [cod2amino[x] for x in wt_DNA_seq];


# Performing evolution and writing files
println("Performing evolution and writing files")

start_msa = hcat([wt_DNA_seq for i in 1:Nseqs]...); 
res = run_evolution(start_msa, 
    h, 
    J,
    p = 0., 
    temp = temp, 
    N_steps = Nsteps,  
    verbose = false);


# Writing the result in a file
println("Writing the result in a file")

dir_name = string(wt_name)*"_MCMC_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)"
if isdir(dir_name)
    # Rimuove tutti i contenuti della cartella
    for file in readdir(dir_name; join=true)
        isdir(file) ? rm(file; recursive=true) : rm(file)
    end
else
    # Se non esiste, la crea
    mkdir(dir_name)
end

file_msa = "../simulation_data/" *string(wt_name)*"_MCMC_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps).fa"
file_en = "../simulation_data/en_" *string(wt_name)*"_MCMC_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps).txt"
file_ham = "../simulation_data/ham_" *string(wt_name)*"_MCMC_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps).txt"
file_pair_ham = "../simulation_data/pair_ham_" *string(wt_name)*"_MCMC_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps).txt"

matrix2fasta(file_msa, Int8.(res.msa'))
writedlm(file_en, Genie.energy(res.msa, h, J))
writedlm(file_ham, ham_dist(wt_seq, res.msa))
writedlm(file_pair_ham, pairwise_ham_dist(res.msa, n_seq = Npairs, all = true))



