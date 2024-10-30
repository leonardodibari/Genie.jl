

function conn_corr(msa, nat_msa, q, L)

    f1_nat_ql, f2_nat_ql = compute_weighted_frequencies(nat_msa, 22, 0.2); f1_nat = reshape(f1_nat_ql,q,L); f2_nat = reshape(f2_nat_ql,q,L,q,L);
    f1_ql, f2_ql = compute_weighted_frequencies(msa, 22, 0.); f1 = reshape(f1_ql,q,L); f2 = reshape(f2_ql,q,L,q,L);   

    c_ql = triu(f2_ql - f1_ql * f1_ql', 21); c_nat_ql = triu(f2_nat_ql - f1_nat_ql * f1_nat_ql', 21);

    c_nat = []; c = []; 
    for i in 1:L
        for j in 1:L
            if i !== j
                for a in 1:q
                    for b in 1:q
                        push!(c_nat, f2_nat[a,i,b,j] - f1_nat[a,i] * f1_nat[b,j])
                        push!(c, f2[a,i,b,j] - f1[a,i] * f1[b,j])
                    end
                end
            end
        end
    end


    println(cor(c_nat[:], c[:]))
    println(cor(c_nat_ql[:], c_ql[:]))
    
end









function check_equilibration(folder::String, nat_msa::Array{Int8,2}, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1})
    
    N_msa = length(steps)
    L, M = size(step_msa[end])
    
    if size(nat_msa,1) !== L
        error("Length of natural and artificial sequences is different \n
            remember that you need to have the msa in a (L, M) format of Int8")
    elseif N_msa !== length(step_msa)
        error("Length of step_msa and of steps is different")
    end
    
    f1_nat,f2_nat = compute_weighted_frequencies(nat_msa, 22, 0.2); conn_nat = triu(f2_nat - f1_nat * f1_nat', 21);
    w = compute_weights(nat_msa, 22, 0.2)[1];
    dist_nat = mean(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    err_dist_nat = std(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    cor1 = []; cor2 = []; corconn = []; dist = []; err_dist = [];
    a = zeros(21*L,21*L);
    for i in 1:N_msa
        f1,f2 = compute_weighted_frequencies(step_msa[i], 22, 0.); 
        conn = triu(f2 - f1 * f1',21);
        if i == N_msa
            a .= conn
        end
        push!(cor1, cor(f1[:], f1_nat[:]))
        push!(cor2, cor(f2[:], f2_nat[:]))
        push!(corconn, cor(conn[:], conn_nat[:]))
        dists = ham_dist(step_msa[1], step_msa[i])       
        push!(dist, mean(dists))
        push!(err_dist, var(dists))
    end 
    
    f1_end,f2_end = DCAUtils.compute_weighted_frequencies(step_msa[end], 22, 0.); 
    
    pc_nat, PC, expl_var = perform_pca(nat_msa, n_dim = 2);
    pc_sil = get_projection(PC, step_msa[end]);
            
    close("all")
    plt.plot(steps, cor1, label = "1-point"); plt.plot(steps, cor2, label = "2-point"); plt.plot(steps, corconn, label = "Conn cor"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Pearson correlation"); plt.xscale("log");savefig(joinpath(folder,"freq_corr_evol.png"));
    
    
    close("all")
    plt.plot(f1_nat[:], f1_end[:], "o", label = "1-point"); plt.legend(); plt.xlabel("Natural frequencies"); plt.ylabel("Artificial frequencies"); savefig(joinpath(folder,"final_freq_natvssim.png"));
    
    
    close("all")
    plt.plot(steps, dist./L, label = "Simulated"); plt.plot(steps, [dist_nat for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Mean hamming from wt (%)"); plt.xscale("log");savefig(joinpath(folder,"mean_hamming.png"))

    
    close("all")
    plt.plot(steps, (err_dist)./L, label = "Simulated"); plt.plot(steps, [err_dist_nat^2 for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Variance hamming from wt (%)");plt.xscale("log"); savefig(joinpath(folder,"var_hamming.png"))
    
    close("all")
    plot_density(pc_nat[:,1], pc_nat[:,2])
    plt.title("Natural sequences")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    savefig(joinpath(folder, "pca_nat.png"))
    
    
    close("all")
    plot_density(pc_sil[:,1], pc_sil[:,2])
    plt.title("Simulated sequences")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    savefig(joinpath(folder, "pca_sim.png"))
    
end
    

function check_equilibration(folder::String, nat_msa_path::String, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1})
    
    nat_msa = read_fasta_alignment(nat_msa_path, 0.9)
    check_equilibration(folder, nat_msa, step_msa, steps)
end