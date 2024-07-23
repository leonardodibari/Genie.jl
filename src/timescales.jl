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
    
    