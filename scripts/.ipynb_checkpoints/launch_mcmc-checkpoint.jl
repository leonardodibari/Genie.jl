using Genie, DelimitedFiles, DCAUtils, JLD2, PyPlot, Statistics, LinearAlgebra
import KitMSA: fasta2matrix, matrix2fasta

wt_name = "TEM"
Nseqs = 100
Nsteps = 100
Npairs = 100
temp = 1.0
file_wt = "../../data/TEM1_pfam_DNA.fasta"
file_model = "../../data/model_natural_PSE.dat"
output_folder = "../../"

println("Reading arguments")
L_args = length(ARGS)
@assert L_args > 0 "Please specify the info for the simulation"
for i in 1:L_args
    if ARGS[i] == "-wt_name"
        global wt_name = ARGS[i+1]
    elseif ARGS[i] == "-Nseqs"
        global Nseqs = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-Nsteps"
        global Nsteps = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-Npairs"
        global Npairs = parse(Int,ARGS[i+1])
    elseif ARGS[i] == "-temp"
        global temp = parse(Float64,ARGS[i+1])     
    elseif ARGS[i] == "-file_wt"
        global file_wt = ARGS[i+1]
    elseif ARGS[i] == "-file_model"
        global file_model = ARGS[i+1]
    elseif ARGS[i] == "-output_folder"
	global output_folder = ARGS[i+1]       
    end
end

println(output_folder*string(wt_name)*"_mcmc_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(Npairs)")


# Read model parameters
println("Reading model parameters")
if wt_name == "AAC" || wt_name == "DHFR"
    	h_tmp, J_tmp = read_par_BM_0gapsave(file_model)
else
	h_tmp, J_tmp = read_par_BM_lettersave(file_model)
end
h = set_max_field_to_0save(h_tmp);
J_tmp2 = symmetrize_Jsave(J_tmp); 
J = permutedims(J_tmp2, [1,3,2,4]);



# Reading sequences in the seed
println("Reading sequences in the seed")
wt = join(readdlm(file_wt, skipstart = 1));
L_big = round(Int,length(wt)/3)
if wt_name == "AAC" || wt_name == "DHFR"
	wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 1:L_big];
else
	wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 3:L_big-1];
end
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
println("Writing the result in files")

file_msa = output_folder*string(wt_name)*"_mcmc_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(Npairs).fa"
file_en = output_folder*"en_"*string(wt_name)*"_mcmc_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(Npairs).txt"
file_ham = output_folder*"ham_"*string(wt_name)*"_mcmc_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(Npairs).txt"
file_pair_ham = output_folder*"pair_ham_" *string(wt_name)*"_mcmc_samples"*string(Nseqs)*"_T"*string(temp)* "_steps$(Nsteps)_pairs$(Npairs).txt"

matrix2fasta(file_msa, Int8.(res.msa'))
writedlm(file_en, Genie.energy(res.msa, h, J))
writedlm(file_ham, ham_dist(wt_seq, res.msa))
writedlm(file_pair_ham, pairwise_ham_dist(res.msa, n_seq = Npairs, all = true))



