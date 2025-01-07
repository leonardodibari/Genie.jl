
# timescales.jl
# Functions for calculating timescales related to correlation and relaxation.


"""
    check_timescales(step_msa::Array{Array{Int8,2},1}, nat_msa::Array{Int8,2}, 
                     h::Array{T,2}, J::Array{T,4}, steps::Array{Int,1}, 
                     filepath::String) where {T}

Generates and saves plots to analyze timescales based on provided multiple sequence alignments (MSA).

# Arguments
- `step_msa`: An array of multiple sequence alignments at different evolutionary steps. *(non-optional)*
- `nat_msa`: The natural MSA used for comparison. *(non-optional)*
- `h`: A 2D array of external field parameters. *(non-optional)*
- `J`: A 4D array of interaction parameters, should be of shape `(q, L, q, L)`. *(non-optional)*
- `steps`: An array containing the Monte Carlo steps for the analysis. *(non-optional)*
- `filepath`: The path where the generated plots will be saved. *(non-optional)*

# Returns
This function does not return a value; it saves the generated plots to the specified file path.

# Errors
- May throw errors related to plotting or file-saving issues, depending on the provided `filepath`.
"""

function check_timescales(step_msa::Array{Array{Int8,2},1}, nat_msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4}, steps::Array{Int,1}, filepath::String) where {T}
    
    L, M = size(step_msa[1])
    
    wt = step_msa[1][:,1]
    
    cie = CIE(nat_msa)
    
    cde_wt = cont_dep_entr(wt, h, J)
    
    freqs = [reshape(compute_weighted_frequencies(x,22,0.)[1],(21, L)) 
        for x in step_msa]
    entr = [get_entropy(f) for f in freqs];
    entr = hcat(entr...)';
    
  
    epis= [x in sortperm(cie .- cde_wt, rev=true)[1:10] for x in 1:length(cie)]
    d_from_neg_bisec = (cde_wt .+ cie .- 2) ./ sqrt(2)
    varr = [x in sortperm(d_from_neg_bisec, rev=true)[1:10] for x in 1:length(cie)]
    cons = [x in sortperm(d_from_neg_bisec, rev=false)[1:10] for x in 1:length(cie)]
    X = steps;

    pointsize = 80
    transp = 0.8
    transp2 = 1
    lab = [0,1,2,3,4]
    lab2 = [10^i for i in 2:7]
    axis_width = 2.

    ticks_font = 22
    axis_font = 26
    
    close("all")

    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(20)
 
    ax1 = plt.subplot2grid(shape=(1, 3), loc=(0, 0), colspan = 1)
    ax2 = plt.subplot2grid(shape=(1, 3), loc=(0, 1), colspan = 2)

    fig.subplots_adjust(wspace=0.2)

    ax1.plot([0, 4.0], [0, 4.0], linestyle="--", alpha = 0.4, color = "darkgrey")
    ax1.scatter(cde_wt, cie, color = "grey", s = pointsize, alpha = transp)
    ax1.scatter(cde_wt[varr], cie[varr], color = "red", s = pointsize, alpha = transp2)
    ax1.scatter(cde_wt[epis], cie[epis], color = "blue",  s = pointsize, alpha = transp2)
    ax1.scatter(cde_wt[cons], cie[cons], color = "green",  s = pointsize, alpha = transp2)
    ax1.set_xlabel("\$CDE_{WT_{background}}\$", fontsize=axis_font)
    ax1.set_ylabel("CIE", fontsize=axis_font)
    ax1.grid(color = "grey", linestyle = "--", linewidth = 0.5, alpha = 0.4)
    ax1.set_xticks([0,1,2,3,4])
    ax1.set_yticks([0,1,2,3,4])
    ax1.set_xticklabels(lab, fontsize = ticks_font)
    ax1.set_yticklabels(lab, fontsize = ticks_font)
    ax1.spines["top"].set_visible(false)
    ax1.spines["right"].set_visible(false)
    ax1.spines["left"].set_linewidth(axis_width)
    ax1.spines["bottom"].set_linewidth(axis_width)
    ax1.set_aspect("equal")
    ax1.set_title("wt background", fontsize = axis_font)

    ax2.plot(X, entr, color = "grey", alpha = 0.35)
    ax2.plot(X, entr[:, varr], color = "red", alpha = 0.35)
    ax2.plot(X, entr[:, epis], color = "blue", alpha = 0.35)
    ax2.plot(X, entr[:, cons], color = "green", alpha = 0.35)
    ax2.plot(X, mean(entr[:, varr], dims = 2)[:], color = "red", linewidth=4.5, label = "Variable" )
    ax2.plot(X, mean(entr[:, epis], dims = 2)[:], color = "blue", linewidth = 4.5, label = "Epistatic")
    ax2.plot(X, mean(entr[:, cons], dims = 2)[:], color = "green", linewidth = 4.5, label = "Conserved")
    ax2.plot(X, mean(entr[:, varr], dims = 2)[:], color = "darkred", linewidth=4.5 )
    ax2.plot(X, mean(entr[:, epis], dims = 2)[:], color = "midnightblue", linewidth = 4.5)
    ax2.plot(X, mean(entr[:, cons], dims = 2)[:], color = "darkgreen", linewidth = 4.5)
    ax2.set_xlabel("Monte Carlo steps", fontsize=axis_font)
    ax2.set_ylabel("Site Entropy", fontsize=axis_font)
    ax2.grid(color = "grey", linestyle = "--", linewidth = 0.5, alpha = 0.4)
    ax2.set_yticks([0,1,2,3,4])
    ax2.set_xticks(lab2)
    ax2.set_xticklabels(lab2, fontsize = ticks_font)
    ax2.set_yticklabels(lab, fontsize = ticks_font)
    ax2.set_xlim(minimum(X),maximum(X))
    ax2.spines["top"].set_visible(false)
    ax2.spines["right"].set_visible(false)
    ax2.spines["left"].set_linewidth(axis_width)
    ax2.spines["bottom"].set_linewidth(axis_width)
    ax2.set_xscale("log")

    tight_layout()
    fig.legend(loc="upper right", fontsize = axis_font, frameon = false, ncol = 3, bbox_to_anchor=(0.6, 1.2))
    
    savefig(filepath)
    
    
end


"""
    info_timescales(step_msa::Array{Array{Int8,2},1}, nat_msa::Array{Int8,2}, 
                     h::Array{T,2}, J::Array{T,4}, steps::Array{Int,1}) where {T}

Analyzes timescales from the given multiple sequence alignments and returns relevant information.

# Arguments
- `step_msa`: An array of multiple sequence alignments at different evolutionary steps. *(non-optional)*
- `nat_msa`: The natural MSA used for comparison. *(non-optional)*
- `h`: A 2D array of external field parameters. *(non-optional)*
- `J`: A 4D array of interaction parameters, should be of shape `(q, L, q, L)`. *(non-optional)*
- `steps`: An array containing the Monte Carlo steps for the analysis. *(non-optional)*

# Returns
A tuple containing:
- `epis`: Indices of epistatic sites.
- `varr`: Indices of variable sites.
- `cons`: Indices of conserved sites.
- `cde_wt`: Conditional dependency entropy of the wild type.
- `entr`: Entropy values corresponding to the MSAs.

# Errors
- May throw errors related to mismatched sizes of input arrays.
"""


function info_timescales(step_msa::Array{Array{Int8,2},1}, nat_msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4}, steps::Array{Int,1}) where {T}
    
    L, M = size(step_msa[1])
    
    wt = step_msa[1][:,1]
    
    cie = CIE(nat_msa)
    
    cde_wt = cont_dep_entr(wt, h, J)
    
    freqs = [reshape(compute_weighted_frequencies(x,22,0.)[1],(21, L)) 
        for x in step_msa]
    entr = [get_entropy(f) for f in freqs];
    entr = hcat(entr...)';
    
    epis= [x in sortperm(cie .- cde_wt, rev=true)[1:10] for x in 1:length(cie)]
    d_from_neg_bisec = (cde_wt .+ cie .- 2) ./ sqrt(2)
    varr = [x in sortperm(d_from_neg_bisec, rev=true)[1:10] for x in 1:length(cie)]
    cons = [x in sortperm(d_from_neg_bisec, rev=false)[1:10] for x in 1:length(cie)]
    
    
    
    return epis, varr, cons, cde_wt, entr
    
end
    
   
"""
    pair_dist_freq(msa; n_seq::Int = 500)

Calculates the frequency of pairwise Hamming distances in a multiple sequence alignment.

# Arguments
- `msa`: The multiple sequence alignment matrix. *(non-optional)*
- `n_seq`: Number of sequences to consider for pairwise comparisons (default is `500`). *(optional)*

# Returns
- An array of frequencies for each pairwise Hamming distance.

# Errors
- Throws an error if the input MSA is not properly formatted.
"""


function pair_dist_freq(msa; n_seq::Int = 500)
    L = size(msa,1)
    d = pairwise_ham_dist(Int8.(msa), n_seq = n_seq, all = true); 
    count_d = [sum(d .== i) for i in 1:L] ./length(d)
    return count_d
end
    
  
"""
    check_pairwise(f_nat::Array{T,1}, f_sim::Array{T,1}, save_path::String) where T

Generates a plot comparing the frequency of pairwise Hamming distances for natural and simulated sequences.

# Arguments
- `f_nat`: Frequencies for natural sequences. *(non-optional)*
- `f_sim`: Frequencies for simulated sequences. *(non-optional)*
- `save_path`: The path where the generated plot will be saved. *(non-optional)*

# Returns
This function does not return a value; it saves the generated plot to the specified file path.

# Errors
- May throw errors related to plotting or file-saving issues, depending on the provided `save_path`.
"""

function check_pairwise(f_nat::Array{T,1}, f_sim::Array{T,1}, save_path::String) where T
    close("all")
    # Create a new figure
    plt.plot(f_nat, color="blue", label="nat", linewidth=2)
    # Plot the second histogram for 'sim' with no fill
    plt.plot(f_sim, color="red", label="sim", linewidth=2)

    plt.yscale("log")
    # Add labels and legend
    plt.xlabel("Pairwise Hamming")
    plt.ylabel("Frequency")
    plt.legend()

    # Save the plot to the specified path
    plt.savefig(save_path)
end


function check_energy(path::String, nat_msa::Array{Int8,2}, msa::Array{Int8,2}, h::Array{Float64,2}, J::Array{Float64,4})
    en_nat = energy(nat_msa, h, J)
    en_sim = energy(msa, h, J)
    
    close("all")
    # Plot the first histogram for 'nat' with no fill
    plt.hist(en_nat, bins=30, histtype="step", density = true, color="blue", label="nat", linewidth=2)

    # Plot the second histogram for 'sim' with no fill
    plt.hist(en_sim, bins=30, histtype="step", density = true, color="red", label="sim", linewidth=2)

    # Add labels and legend
    plt.xlabel("Energy")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.legend()

    # Save the plot to the specified path
    plt.savefig(path)
end


function check_energy(path::String, nat_msa::Array{Int8,2}, msa::Array{Int8,2}, h::Array{Float64,2}, J::Array{Float64,4}, w::Array{Float64,1})
    en_nat = energy(nat_msa, h, J)
    en_sim = energy(msa, h, J)
    
    close("all")
    # Plot the first histogram for 'nat' with no fill
    plt.hist(en_nat, weights = w, bins=30, histtype="step", density = true, color="blue", label="nat", linewidth=2)

    # Plot the second histogram for 'sim' with no fill
    plt.hist(en_sim, bins=30, histtype="step", color="red", density = true, label="sim", linewidth=2)

    # Add labels and legend
    plt.xlabel("Energy")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.legend()

    # Save the plot to the specified path
    plt.savefig(path)
end



function check_dist_from_wt(wt::Array{T,1}, msa::Array{T,2}, msa_nat::Array{T,2}, save_path::String) where {T}
    L = length(wt)
    d_nat = ham_dist(wt, msa_nat) ./L
    d_sim =  ham_dist(wt, msa) ./L
    
    close("all")
    
     # Plot the first histogram for 'nat' with no fill
    plt.hist(d_nat, bins=30, histtype="step", density = true, color="blue", label="nat", linewidth=2)

    # Plot the second histogram for 'sim' with no fill
    plt.hist(d_sim, bins=30, histtype="step", color="red", density = true, label="sim", linewidth=2)

    plt.yscale("log")
    # Add labels and legend
    plt.xlabel("Hamming from Wt")
    plt.ylabel("Frequency")
    plt.legend()

    # Save the plot to the specified path
    plt.savefig(save_path)
end