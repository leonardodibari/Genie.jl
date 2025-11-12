using DCAUtils, PyPlot, Statistics, DelimitedFiles
using KitMSA

# ------------------------------------------------------------
# Compare replicates (scatter plots f1 and f2)
# ------------------------------------------------------------
function compare_replicates(file1::String, file2::String, title_str::String)
    println("\n🔹 [$title_str] Loading alignments...")
    aln1 = read_fasta_alignment(file1, 0.9)
    aln2 = read_fasta_alignment(file2, 0.9)

    println("🔹 [$title_str] Computing frequencies...")
    f1_a, f2_a = compute_weighted_frequencies(aln1, 22, 0.0)
    f1_b, f2_b = compute_weighted_frequencies(aln2, 22, 0.0)

    println("🔹 [$title_str] Computing correlations...")
    corr_f1 = cor(f1_a[:], f1_b[:])
    corr_f2 = cor(f2_a[:], f2_b[:])

    println("🔹 [$title_str] Plotting scatter plots...")
    figure(figsize=(10,5))

    subplot(1,2,1)
    scatter(f1_a[:], f1_b[:], s=10, alpha=0.6)
    xscale("log"); yscale("log")
    title("$title_str - f1")
    xlabel("Replicate 1"); ylabel("Replicate 2")
    text(0.05, 0.95, "corr = $(round(corr_f1, digits=3))", transform=gca().transAxes, va="top", ha="left")

    subplot(1,2,2)
    scatter(f2_a[:], f2_b[:], s=10, alpha=0.6)
    xscale("log"); yscale("log")
    title("$title_str - f2")
    xlabel("Replicate 1"); ylabel("Replicate 2")
    text(0.05, 0.95, "corr = $(round(corr_f2, digits=3))", transform=gca().transAxes, va="top", ha="left")

    tight_layout()
    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_scatter_log.png", dpi=300)
    println("✅ Figure saved as $(filename)_scatter_log.png")
    show()

    return corr_f1, corr_f2
end

# ------------------------------------------------------------
# Plot ranked frequencies
# ------------------------------------------------------------
function plot_frequencies_vs_rank(files::Vector{String}, labels::Vector{String}, colors::Vector{String}, title_str::String)
    println("\n🔹 [$title_str] Plotting frequencies vs rank...")
    figure(figsize=(12,5))

    # --- f1 subplot ---
    subplot(1,2,1)
    for (file, lbl, col) in zip(files, labels, colors)
        aln = read_fasta_alignment(file, 0.9)
        f1, f2 = compute_weighted_frequencies(aln, 22, 0.0)
        freqs1 = sort(f1[:], rev=true)
        scatter(1:length(freqs1), freqs1, s=8, alpha=0.6, label=lbl, color=col)
    end
    yscale("log"); xlabel("Rank"); ylabel("Frequency")
    title("$title_str - 1-point"); grid(true, which="both"); legend()

    # --- f2 subplot ---
    subplot(1,2,2)
    for (file, lbl, col) in zip(files, labels, colors)
        aln = read_fasta_alignment(file, 0.9)
        f1, f2 = compute_weighted_frequencies(aln, 22, 0.0)
        freqs2 = sort(f2[:], rev=true)
        scatter(1:length(freqs2), freqs2, s=8, alpha=0.6, label=lbl, color=col)
    end
    yscale("log"); xlabel("Rank"); ylabel("Frequency")
    title("$title_str - 2-point"); grid(true, which="both"); legend()

    tight_layout()
    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_f1_f2_vs_rank.png", dpi=300)
    println("✅ Figure saved as $(filename)_f1_f2_vs_rank.png")
    show()
end



# ------------------------------------------------------------
# Plot site frequency spectrum
# ------------------------------------------------------------

using StatsBase  # for Histogram

function plot_site_frequency_spectrum_1point(files::Vector{String}, labels::Vector{String}, colors::Vector{String}, title_str::String)

    
    println("\n🔹 [$title_str] Plotting 1-point site frequency spectrum (SFS)...")
    figure(figsize=(6,5))

    for (file, lbl, col) in zip(files, labels, colors)
        aln = read_fasta_alignment(file, 0.9)
        f1, f2 = compute_weighted_frequencies(aln, 22, 0.0)

        # flatten and remove zeros
        freqs1 = vec(f1)
        freqs1 = freqs1[freqs1 .> 0]

        # histogram with 100 bins
        h = fit(Histogram, freqs1; nbins=100)
        mids = (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2
        counts = h.weights

        #scatter(mids, counts, s=8, alpha=0.6, label=lbl, color=col)
        plot(mids, counts, alpha=0.6, label=lbl, color=col)
    end

    xscale("log"); yscale("log")
    xlabel("Frequency"); ylabel("Count")
    title("$title_str - 1-point SFS")
    grid(true, which="both"); legend()

    tight_layout()
    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_sfs_1point.png", dpi=300)
    println("✅ Figure saved as $(filename)_sfs_1point.png")
    show()
end


# NOTE: The functions 'cod2amino', 'read_fasta_alignment', and 'compute_weighted_frequencies' 
# are assumed to be defined elsewhere.
const q_states = 21 # The true number of states (q) for indexing the flattened frequency vector

function plot_site_frequency_spectrum_1point_no_wt(wt_name::String, ref_file::String, files::Vector{String}, labels::Vector{String}, colors::Vector{String}, title_str::String)

    println("Reading sequences in the seed")
    wt = join(readdlm(ref_file, skipstart = 1)); # Read the reference sequence (DNA)
    L_big = round(Int,length(wt)/3)

    # 1. Determine the relevant DNA codons
    if wt_name == "AAC" || wt_name == "DHFR"
        wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 1:L_big];
    else
        wt_DNA_seq = [wt[((i-1)*3 +1):(i*3)] for i in 3:L_big-1];
    end
    
    # 2. Convert DNA codons to the integer-encoded wild-type sequence (wt_seq)
    L_sites = length(wt_DNA_seq)
    
    # ***PLACEHOLDER: Replace this line with your actual integer-encoding logic***
    # wt_seq elements must be integers from 1 to q_states (1 to 21).
    wt_seq = [cod2amino[x] for x in wt_DNA_seq];

    println("\n🔹 [$title_str] Plotting 1-point site frequency spectrum (SFS) (excluding WT frequencies)...")
    figure(figsize=(6,5))

    # The function is called with 22, so we keep that in the call, but use 21 internally.
    q_param_call = q_states + 1 

    for (file, lbl, col) in zip(files, labels, colors)
        aln = read_fasta_alignment(file, 0.9) 
        # f1 is a 1D vector of length L*21
        f1, f2 = compute_weighted_frequencies(aln, q_param_call, 0.0) 
        
        # --- MODIFICATION START: Filter out WT frequencies from the 1D vector ---
        
        # Check size consistency
        expected_length = L_sites * q_states
        if length(f1) != expected_length
            error("Frequency vector f1 has length $(length(f1)), but expected $expected_length for $L_sites sites and $q_states states. Check your 'compute_weighted_frequencies' output.")
        end

        ff1 = reshape(f1, q_states, L_sites)

        for (i, aa) in enumerate(wt_seq)
            ff1[aa,i] = NaN
        end

        freqs1 = filter(isfinite, vec(ff1))

                
        # The remaining part is the standard plotting logic
        
        # histogram with 100 bins
        h = fit(Histogram, freqs1; nbins=100)
        mids = (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2
        counts = h.weights

        plot(mids, counts, linewidth = 3.0, alpha=0.6, label=lbl, color=col)
    end

    xscale("log"); yscale("log")
    xlabel("Frequency", fontsize = 15); ylabel("Count", fontsize = 15)
    #title("SFS (Mutations Only)")
    legend()

    tight_layout()
    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_sfs_1point_no_wt.png", dpi=300)
    println("✅ Figure saved as $(filename)_sfs_1point_no_wt.png")
    show()
end

# ------------------------------------------------------------
# Compute site entropy
# ------------------------------------------------------------
function site_entropy(f1_vec::Vector{Float64}, q::Int; pseudocount=1e-6)
    L = div(length(f1_vec), q)
    f1 = reshape(copy(f1_vec), q, L)
    λ = pseudocount
    for l in 1:L, a in 1:q
        f1[a,l] = (f1[a,l] + λ/q) / (1 + λ)
    end
    S = zeros(L)
    for l in 1:L, a in 1:q
        p = f1[a,l]
        if p > 0
            S[l] -= p * log2(p)
        end
    end
    return S
end

# ------------------------------------------------------------
# Compare site entropy between two alignments
# ------------------------------------------------------------
function compare_entropy(file1::String, file2::String, title_str::String; pseudocount=1e-6)
    println("\n🔹 [$title_str] Loading alignments...")
    aln1 = read_fasta_alignment(file1, 0.9)
    aln2 = read_fasta_alignment(file2, 0.9)

    println("🔹 [$title_str] Computing frequencies...")
    f1_a, _ = compute_weighted_frequencies(aln1, 22, 0.0)
    f1_b, _ = compute_weighted_frequencies(aln2, 22, 0.0)

    println("🔹 [$title_str] Computing site entropies...")
    S_a = site_entropy(f1_a, 22; pseudocount=pseudocount)
    S_b = site_entropy(f1_b, 22; pseudocount=pseudocount)

    corr_S = cor(S_a, S_b)

    println("🔹 [$title_str] Plotting entropy scatter...")
    figure(figsize=(5,5))
    scatter(S_a, S_b, s=15, alpha=0.6)
    title("$title_str - site entropy")
    xlabel("Replicate 1 entropy (bits)")
    ylabel("Replicate 2 entropy (bits)")
    text(0.05, 0.95, "corr = $(round(corr_S, digits=3))", transform=gca().transAxes, va="top", ha="left")
    tight_layout()

    filename = lowercase(replace(title_str, " " => "_"))
    savefig("$(filename)_entropy.png", dpi=300)
    println("✅ Entropy figure saved as $(filename)_entropy.png")
    show()

    return corr_S
end

# ------------------------------------------------------------
# Example usage (adapt paths)
# ------------------------------------------------------------

#=sequence_overlap(
    "../simulation_data/PSE/popgen/PSE_popgen_samples1000000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",
    "../simulation_data/PSE/popgen/PSE_popgen_samples1000000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
    "Pop Gen (10^6)"
)

sequence_overlap(
    "../simulation_data/PSE/mcmc/PSE_mcmc_samples100000_T1.3_steps120_pairs10000.fa",
    "../simulation_data/PSE/mcmc/PSE_mcmc_samples100000_T1.31_steps120_pairs10000.fa",
    "MCMC"
)

sequence_overlap(
    "../simulation_data/PSE/mcmctree/PSE_mcmctree_T1.4_steps110_pairs10000.fa",
    "../simulation_data/PSE/mcmctree/PSE_mcmctree_T1.41_steps110_pairs10000.fa",
    "MCMC on tree"
)

sequence_overlap(
    "../simulation_data/PSE/popgen/PSE_popgen_samples100000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",
    "../simulation_data/PSE/popgen/PSE_popgen_samples100000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
    "Pop Gen"
)

    # --- Compare replicates ---
  #=  compare_replicates(
        "../PSE/popgen/PSE_popgen_samples1000000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",
        "../PSE/popgen/PSE_popgen_samples1000000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
        "Pop Gen (10^6)"
    )

    compare_replicates(
        "../PSE/mcmc/PSE_mcmc_samples100000_T1.3_steps120_pairs10000.fa",
        "../PSE/mcmc/PSE_mcmc_samples100000_T1.31_steps120_pairs10000.fa",
        "MCMC"
    )


compare_replicates(
           "../PSE/mcmctree/PSE_mcmctree_T1.4_steps110_pairs10000.fa",
               "../PSE/mcmctree/PSE_mcmctree_T1.41_steps110_pairs10000.fa",
                   "MCMC on tree"
       )

compare_replicates(
     "../PSE/popgen/PSE_popgen_samples100000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",      "../PSE/popgen/PSE_popgen_samples100000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
     "Pop Gen"
	) =#
   =# 
# --- Frequencies vs rank ---
    ref_file = "../data/PSE1_pfam_DNA.fasta"
    wt_name = "PSE"
    files = [
        "../data/pse_round20_noflankgaps_nodupl.fa",
        "../simulation_data/PSE/mcmc/PSE_mcmc_samples100000_T1.3_steps120_pairs10000.fa",
        "../simulation_data/PSE/popgen/PSE_popgen_samples100000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",
        "../simulation_data/PSE/mcmctree/PSE_mcmctree_T1.4_steps110_pairs10000.fa",
        "../simulation_data/PSE/popgen/PSE_popgen_samples10000_T0.0_mu0.03_mubind-50.0_pairs10000.fa"
    ]
    labels = ["Exp", "MCMC", "Pop Gen", "MCMC on tree", "Neutral evolution"]
    colors = ["black", "red", "green", "blue", "grey"]


    plot_site_frequency_spectrum_1point_no_wt(wt_name, ref_file, files, labels, colors, "PSE")

    #plot_frequencies_vs_rank(files, labels, colors, "PSE")



    # --- Compare site entropies ---
  #=  compare_entropy(
        "../PSE/popgen/PSE_popgen_samples1000000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",
        "../PSE/popgen/PSE_popgen_samples1000000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
        "Pop Gen Entropy"
    )

compare_entropy(
     "../PSE/popgen/PSE_popgen_samples1000000_T12.0_mu0.03_mubind-50.0_pairs10000.fa",     
     "../PSE/popgen/PSE_popgen_samples1000000_T12.01_mu0.03_mubind-50.0_pairs10000.fa",
     "Pop Gen (10^6)"
)

compare_entropy(
    "../PSE/mcmc/PSE_mcmc_samples100000_T1.3_steps120_pairs10000.fa",
    "../PSE/mcmc/PSE_mcmc_samples100000_T1.31_steps120_pairs10000.fa",
    "MCMC"
)

compare_entropy(
   "../PSE/mcmctree/PSE_mcmctree_T1.4_steps110_pairs10000.fa",
   "../PSE/mcmctree/PSE_mcmctree_T1.41_steps110_pairs10000.fa",
   "MCMC on tree"
) =#


