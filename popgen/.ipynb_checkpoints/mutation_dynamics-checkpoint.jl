ENV["PYTHON"] = "" 
using Pkg
Pkg.build("PyCall")  # Rebuild PyCall to use the internal Python
using Genie, DelimitedFiles, DCAUtils, JLD2, PyPlot, Statistics, LinearAlgebra
import KitMSA: fasta2matrix, matrix2fasta

ENV["PYTHON"] = "" 
using Pkg
Pkg.build("PyCall")  # Rebuild PyCall to use the internal Python
using Genie, DelimitedFiles, DCAUtils, JLD2, PyPlot, Statistics, LinearAlgebra
import KitMSA: fasta2matrix, matrix2fasta



using PyPlot


wt_name = "PSE"
Nseqs = 100000
#file_wt = "../../data/PSE1_pfam_DNA.fasta"
#file_model = "../../data/model_natural_PSE.dat"

file_wt = "../data_directed_evolution/PSE1_pfam_DNA.fasta"
file_model = "../model_natural_PSE.dat"

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


Nsteps = 5*120
temp = 1.3
start_msa = hcat([wt_DNA_seq for i in 1:Nseqs]...); 
mcmc = run_evolution(start_msa, 
    h, 
    J,
    p = 0., 
    temp = temp, 
    N_steps = Nsteps, 
    each_step = 6,
    verbose = false);

f1_mcmc = [];
for x in mcmc.step_msa
    a,b = compute_weighted_frequencies(Int8.(x), 22, 0.)
    push!(f1_mcmc, a[:])
end

f1_mcmc = reduce(vcat, (permutedims(v) for v in f1_mcmc))  # matrix (steps × freqs)

N1, N2 = size(f1_mcmc);



function plot_evolving_frequencies(f1_matrix::AbstractMatrix, title_str::String; threshold=0.1)
    n_steps, n_freqs = size(f1_matrix)

    figure(figsize=(8,6))

    for j in 1:n_freqs
        traj = f1_matrix[:, j]   # frequency trajectory for freq j
        if maximum(traj) ≥ threshold
            plot(1:n_steps, traj, alpha=0.6)
        end
    end

    xlabel("Step", fontsize = 18)
    ylabel("Frequency",fontsize = 18)
    title(title_str, fontsize = 18)
    #grid(true)
    tight_layout()

    # Save figure
    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_trajectories.png", dpi=300)

    println("✅ Figure saved as $(filename)_trajectories.png in $(pwd())")

    show()
end



plot_evolving_frequencies(f1_mcmc, "MCMC")
plot_evolving_frequencies(f1_mcmc[1:round(Int,N1/5),:], "MCMC")


popgen = Genie.run_dir_evol_nucleo(wt_DNA_seq, Nseqs, h, J;
                   rounds = 5*20,
                   each_step = 1,
                   seq_reads = Nseqs,
                   temp = 12.0,
                   mu = 0.03,
                   mu_bind = -50.0);



f1_popgen = [];
for x in popgen.step_msa
    a,b = compute_weighted_frequencies(Int8.(x), 22, 0.)
    push!(f1_popgen, a[:])
end

f1_popgen = reduce(vcat, (permutedims(v) for v in f1_popgen))  # matrix (steps × freqs)

N1, N2 = size(f1_popgen);

plot_evolving_frequencies(f1_popgen, "Pop Gen")
plot_evolving_frequencies(f1_popgen[1:round(Int,N1/5),:], "Pop Gen")




